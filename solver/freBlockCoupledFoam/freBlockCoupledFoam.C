/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    freBlockCoupledFoam

Group
    AcousticSolvers

Description
    Block-coupled frequency-domain acoustic solver for a single, non-
    decomposed mesh. The real/imaginary pressure pair is solved as one
    2x2 block Helmholtz system:

        [ A  -(B1 + B2) ] [Pim] = [bPim]
        [ (B1 + B2)  A ] [Pre]   [bPre]

    where A is the main Helmholtz operator and B1/B2 are PML-induced
    coupling operators. Matrix assembly is local/serial in OpenFOAM and
    the PETSc solve uses a direct backend (default: MUMPS via LU).

    This solver intentionally aborts in decomposed MPI runs. Use it as the
    serial reference for equation form and damping behavior.
\*---------------------------------------------------------------------------*/

#include <petscksp.h>
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "isoAdvection.H"
#include "cutFaceAdvect.H"
#include "surfaceIteratorPLIC.H"
#include "reconstructionSchemes.H"
#include "upwind.H"
#include "processorPolyPatch.H"
#include "processorBC.H"

static inline scalar twoPi() { return constant::mathematical::twoPi; }

// ------------------------------------------------------------------------- //
// Matrix/source helpers
// ------------------------------------------------------------------------- //

// Add scalar boundary source contributions from boundary coefficients.
static void addBoundarySourceSimple
(
    const fvScalarMatrix& M,
    scalarField& source,
    const bool couples = true
)
{
    const auto& bpsi = M.psi().boundaryField();
    const auto& bcoeffs = M.boundaryCoeffs();
    const auto& addr = M.lduAddr();

    forAll(bpsi, patchi)
    {
        const fvPatchField<scalar>& ptf = bpsi[patchi];
        if (!ptf.size()) continue;

        const Field<scalar>& pbc = bcoeffs[patchi];
        const labelUList& paddr = addr.patchAddr(patchi);

        if (!ptf.coupled())
        {
            forAll(paddr, facei)
            {
                source[paddr[facei]] += pbc[facei];
            }
        }
        else if (couples)
        {
            const tmp<Field<scalar>> tpnf = ptf.patchNeighbourField();
            const Field<scalar>& pnf = tpnf();
            forAll(paddr, facei)
            {
                source[paddr[facei]] += pbc[facei]*pnf[facei];
            }
        }
    }
}

// Insert an fvScalarMatrix (ldu form) into one PETSc block position.
static void insertScalarOpIntoBlock
(
    Mat& M,
    const fvScalarMatrix& op,
    const globalIndex& globalCells,
    const PetscInt rowOffset,
    const PetscInt colOffset,
    const PetscScalar scale
)
{
    const lduMatrix& L = op;

    const scalarField& diag = L.diag();

    const lduAddressing& addr = L.lduAddr();
    const labelUList& owner = addr.lowerAddr();
    const labelUList& neigh = addr.upperAddr();

    // Diagonal
    forAll(diag, i)
    {
        PetscInt r = globalCells.toGlobal(i) + rowOffset;
        PetscInt c = globalCells.toGlobal(i) + colOffset;
        MatSetValue(M, r, c, scale*diag[i], INSERT_VALUES);
    }

    // Off-diagonal couplings (if present)
    if (L.hasUpper())
    {
        const scalarField& upper = L.upper();
        const scalarField* lowerPtr = L.hasLower() ? &L.lower() : nullptr;

        forAll(upper, f)
        {
            const label i = owner[f];
            const label j = neigh[f];

            const PetscInt ri = globalCells.toGlobal(i) + rowOffset;
            const PetscInt rj = globalCells.toGlobal(j) + rowOffset;

            const PetscInt ci = globalCells.toGlobal(i) + colOffset;
            const PetscInt cj = globalCells.toGlobal(j) + colOffset;

            const scalar upperVal = upper[f];
            const scalar lowerVal = lowerPtr ? (*lowerPtr)[f] : upperVal;

            MatSetValue(M, ri, cj, scale*upperVal, INSERT_VALUES); // owner->neigh
            MatSetValue(M, rj, ci, scale*lowerVal, INSERT_VALUES); // neigh->owner
        }
    }
}

static scalarField sourceWithBoundary
(
    const fvScalarMatrix& M
)
{
    scalarField source(M.source());
    addBoundarySourceSimple(M, source);
    return source;
}

static void assembleBlockSystem
(
    Mat& M,
    const globalIndex& globalCells,
    const PetscInt N,
    fvScalarMatrix& AopPim,
    fvScalarMatrix& AopPre,
    fvScalarMatrix& B1Pre,
    fvScalarMatrix& B2Pre,
    fvScalarMatrix& B1Pim,
    fvScalarMatrix& B2Pim
)
{
    // Apply boundary contributions to matrices before insertion.
    AopPre.boundaryManipulate(AopPre.psi().boundaryFieldRef());
    AopPim.boundaryManipulate(AopPim.psi().boundaryFieldRef());
    B1Pre.boundaryManipulate(B1Pre.psi().boundaryFieldRef());
    B2Pre.boundaryManipulate(B2Pre.psi().boundaryFieldRef());
    B1Pim.boundaryManipulate(B1Pim.psi().boundaryFieldRef());
    B2Pim.boundaryManipulate(B2Pim.psi().boundaryFieldRef());

    // Coupled system:
    // [ A  -(B1 + B2) ] [Pim] = [bPim]
    // [ (B1 + B2)  A ] [Pre]  [bPre]
    insertScalarOpIntoBlock(M, AopPim, globalCells, 0, 0, +1.0);
    insertScalarOpIntoBlock(M, B1Pre, globalCells, 0, N, -1.0);
    insertScalarOpIntoBlock(M, B2Pre, globalCells, 0, N, -1.0);

    insertScalarOpIntoBlock(M, B1Pim, globalCells, N, 0, +1.0);
    insertScalarOpIntoBlock(M, B2Pim, globalCells, N, 0, +1.0);
    insertScalarOpIntoBlock(M, AopPre, globalCells, N, N, +1.0);
}

static void setBlockRhs
(
    Vec& b,
    const globalIndex& globalCells,
    const PetscInt N,
    const scalarField& bPim,
    const scalarField& bPre
)
{
    forAll(bPim, i)
    {
        const PetscInt gi = globalCells.toGlobal(i);
        VecSetValue(b, gi + 0, (PetscScalar)bPim[i], INSERT_VALUES);
        VecSetValue(b, gi + N, (PetscScalar)bPre[i], INSERT_VALUES);
    }
}

// ------------------------------------------------------------------------- //
// PETSc lifecycle and diagnostics
// ------------------------------------------------------------------------- //

static void initializePetscSystem
(
    const PetscInt nLocal,
    const PetscInt N,
    Mat& M,
    Vec& x,
    Vec& b,
    KSP& ksp
)
{
    MatCreate(PETSC_COMM_WORLD, &M);
    MatSetSizes(M, 2*nLocal, 2*nLocal, 2*N, 2*N);
    MatSetType(M, (Pstream::nProcs() == 1) ? MATSEQAIJ : MATMPIAIJ);
    MatSetUp(M);

    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, 2*nLocal, 2*N);
    VecSetFromOptions(x);
    VecDuplicate(x, &b);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, M, M);

    // Helmholtz-friendly defaults (can be overridden by PETSc options).
    PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");
    PetscOptionsSetValue(nullptr, "-pc_type", "lu");
    PetscOptionsSetValue
    (
        nullptr,
        "-pc_factor_mat_solver_type",
        "mumps"
    );
    KSPSetFromOptions(ksp);
}

static void reportMatrixStats
(
    Mat& M,
    Vec& b
)
{
    PetscReal n1 = 0.0, ninf = 0.0, nfrob = 0.0;
    PetscBool isSym = PETSC_FALSE;
    MatNorm(M, NORM_1, &n1);
    MatNorm(M, NORM_INFINITY, &ninf);
    MatNorm(M, NORM_FROBENIUS, &nfrob);
    MatIsSymmetric(M, 1e-12, &isSym);

    Vec d;
    VecDuplicate(b, &d);
    MatGetDiagonal(M, d);
    PetscReal dmin = 0.0, dmax = 0.0;
    PetscInt imin = 0, imax = 0;
    VecMin(d, &imin, &dmin);
    VecMax(d, &imax, &dmax);
    VecDestroy(&d);

    PetscReal bNorm0 = 0.0;
    VecNorm(b, NORM_2, &bNorm0);

    Info<< "Matrix stats: ||A||1=" << n1
        << " ||A||inf=" << ninf
        << " ||A||F=" << nfrob
        << " diag[min,max]=(" << dmin << "," << dmax << ")"
        << " symmetric=" << (isSym ? "yes" : "no")
        << " rhsNorm=" << bNorm0 << nl;
}

static void reportKspStats
(
    KSP& ksp,
    Vec& b
)
{
    PetscInt kspIters = 0;
    PetscReal kspRes = 0.0;
    PetscReal bNorm = 0.0;
    KSPConvergedReason kspReason;
    KSPGetIterationNumber(ksp, &kspIters);
    KSPGetResidualNorm(ksp, &kspRes);
    KSPGetConvergedReason(ksp, &kspReason);
    VecNorm(b, NORM_2, &bNorm);
    Info<< "KSP iters: " << kspIters
        << " residual: " << kspRes
        << " rhsNorm: " << bNorm
        << " reason: " << kspReason << nl;
}

static void scatterBlockSolution
(
    Vec& x,
    const globalIndex& globalCells,
    const PetscInt N,
    volScalarField& Pim,
    volScalarField& Pre
)
{
    // Block ordering in x: [Pim block, Pre block]
    const PetscScalar* xArr = nullptr;
    VecGetArrayRead(x, &xArr);
    forAll(Pim, i)
    {
        const PetscInt gi = globalCells.toGlobal(i);
        Pim[i] = (scalar)PetscRealPart(xArr[gi + 0]);
        Pre[i] = (scalar)PetscRealPart(xArr[gi + N]);
    }
    VecRestoreArrayRead(x, &xArr);
}

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "computePMLCoefs.H"
    #include "computeAlphaf.H"

    simpleControl simple(mesh);

    PetscInitialize(&argc, &argv, nullptr, nullptr);

    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "freBlockCoupledFoam uses serial assembly. "
            << "Run on one MPI rank and use OMP_NUM_THREADS for multicore."
            << exit(FatalError);
    }

    // Global indexing for block system
    globalIndex globalCells(mesh.nCells());
    const PetscInt N      = (PetscInt)globalCells.size();
    const PetscInt nLocal = (PetscInt)mesh.nCells();

    Mat M;
    Vec x, b;
    KSP ksp;
    bool dumpedMatrixStats = false;

    initializePetscSystem(nLocal, N, M, x, b, ksp);

    Info<< "\nStarting time loop\n" << endl;

    const dimensionedScalar k2
    (
        "k2",
        dimless/dimLength/dimLength,
        sqr(twoPi()*f/cg).value()
    );

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            MatZeroEntries(M);
            VecSet(b, 0.0);

            fvScalarMatrix AopPre
            (
                fvm::laplacian(1 - (rhol - rhog)/rhol*alphaf, Pre)
              + fvm::Sp(k2*(1 + ((kl - kg)/kg)*alpha1) - SC0, Pre)
            );

            fvScalarMatrix AopPim
            (
                fvm::laplacian(1 - (rhol - rhog)/rhol*alphaf, Pim)
              + fvm::Sp(k2*(1 + ((kl - kg)/kg)*alpha1) - SC0, Pim)
            );

            // Coupling operators built from the opposite field:
            // Pim-row couples to Pre via (laplacian(TC1,Pre) + Sp(SC1,Pre))
            fvScalarMatrix couplingLaplPre(fvm::laplacian(TC1, Pre));
            fvScalarMatrix couplingMassPre(fvm::Sp(SC1, Pre));

            // Pre-row couples to Pim via (laplacian(TC1,Pim) + Sp(SC1,Pim))
            fvScalarMatrix couplingLaplPim(fvm::laplacian(TC1, Pim));
            fvScalarMatrix couplingMassPim(fvm::Sp(SC1, Pim));

            assembleBlockSystem
            (
                M, globalCells, N,
                AopPim, AopPre,
                couplingLaplPre, couplingMassPre,
                couplingLaplPim, couplingMassPim
            );

            // RHS from source terms (includes BC contributions)
            scalarField bPim(sourceWithBoundary(AopPim));
            scalarField bCouplingLaplPre(sourceWithBoundary(couplingLaplPre));
            scalarField bCouplingMassPre(sourceWithBoundary(couplingMassPre));
            bPim -= bCouplingLaplPre;
            bPim -= bCouplingMassPre;

            scalarField bPre(sourceWithBoundary(AopPre));
            scalarField bCouplingLaplPim(sourceWithBoundary(couplingLaplPim));
            scalarField bCouplingMassPim(sourceWithBoundary(couplingMassPim));
            bPre += bCouplingLaplPim;
            bPre += bCouplingMassPim;

            setBlockRhs(b, globalCells, N, bPim, bPre);

            MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
            VecAssemblyBegin(b);
            VecAssemblyEnd(b);

            if (!dumpedMatrixStats)
            {
                dumpedMatrixStats = true;
                reportMatrixStats(M, b);
            }

            KSPSolve(ksp, b, x);
            reportKspStats(ksp, b);
            scatterBlockSolution(x, globalCells, N, Pim, Pre);

            Pre.correctBoundaryConditions();
            Pim.correctBoundaryConditions();
        }

        Ure == 1/(2*constant::mathematical::pi*f*rho) * fvc::grad(Pim);
        Uim == -1/(2*constant::mathematical::pi*f*rho) * fvc::grad(Pre);
        pa == Foam::sqrt(Pim*Pim + Pre*Pre);
        pr == 0.25*(kl*alpha1 + kg*(1-alpha1))*(Pre*Pre + Pim*Pim)
            - 0.25*rho*((Ure&Ure) + (Uim&Uim));
        momFlux == 0.5*rho*(Ure*Ure + Uim*Uim);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    KSPDestroy(&ksp);
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&M);
    PetscFinalize();

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
