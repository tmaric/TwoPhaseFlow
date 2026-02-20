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
    grpAcousticSolvers

Description
    Block-coupled frequency-domain acoustic solver (Pre/Pim) using PETSc.
    Cylindrical PML with isotropic damping, scalar PML coefficients.
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

// Add boundary source contributions (scalar-only helper)
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

// Insert an fvScalarMatrix (ldu) into PETSc block matrix
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

    // Global indexing for block system
    globalIndex globalCells(mesh.nCells());
    const PetscInt N      = (PetscInt)globalCells.size();
    const PetscInt nLocal = (PetscInt)mesh.nCells();

    Mat M;
    Vec x, b;
    KSP ksp;
    bool dumpedMatrixStats = false;

    MatCreate(PETSC_COMM_WORLD, &M);
    MatSetSizes(M, 2*nLocal, 2*nLocal, 2*N, 2*N);
    if (Pstream::nProcs() == 1)
    {
        MatSetType(M, MATSEQAIJ);
    }
    else
    {
        MatSetType(M, MATMPIAIJ);
    }
    MatSetUp(M);

    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, 2*nLocal, 2*N);
    VecSetFromOptions(x);
    VecDuplicate(x, &b);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, M, M);
    // Helmholtz-friendly defaults (can be overridden by PETSc options)
    PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");
    PetscOptionsSetValue(nullptr, "-pc_type", "lu");
    if (Pstream::nProcs() == 1)
    {
        // Use PETSc's built-in sequential LU (no external solver required)
        PetscOptionsSetValue(nullptr, "-pc_factor_mat_solver_type", "petsc");
    }
    else
    {
        // Requires PETSc built with MUMPS for parallel LU
        PetscOptionsSetValue(nullptr, "-pc_factor_mat_solver_type", "mumps");
    }
    KSPSetFromOptions(ksp);

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

            // Coupling operators for Pre (used in Pim equation)
            fvScalarMatrix B1Pre(fvm::laplacian(TC1, Pre));
            fvScalarMatrix B2Pre(fvm::Sp(SC1, Pre));

            // Coupling operators for Pim (used in Pre equation)
            fvScalarMatrix B1Pim(fvm::laplacian(TC1, Pim));
            fvScalarMatrix B2Pim(fvm::Sp(SC1, Pim));

            // Apply boundary contributions to matrices
            AopPre.boundaryManipulate(Pre.boundaryFieldRef());
            AopPim.boundaryManipulate(Pim.boundaryFieldRef());
            B1Pre.boundaryManipulate(Pre.boundaryFieldRef());
            B2Pre.boundaryManipulate(Pre.boundaryFieldRef());
            B1Pim.boundaryManipulate(Pim.boundaryFieldRef());
            B2Pim.boundaryManipulate(Pim.boundaryFieldRef());

            // Assemble blocks:
            // [ A  -(B1 + B2) ] [Pim] = [bPim]
            // [ (B1 + B2)  A ] [Pre]  [bPre]
            insertScalarOpIntoBlock(M, AopPim, globalCells, 0, 0, +1.0);
            insertScalarOpIntoBlock(M, B1Pre, globalCells, 0, N, -1.0);
            insertScalarOpIntoBlock(M, B2Pre, globalCells, 0, N, -1.0);

            insertScalarOpIntoBlock(M, B1Pim, globalCells, N, 0, +1.0);
            insertScalarOpIntoBlock(M, B2Pim, globalCells, N, 0, +1.0);
            insertScalarOpIntoBlock(M, AopPre, globalCells, N, N, +1.0);

            // RHS from source terms (includes BC contributions)
            scalarField bPim(AopPim.source());
            addBoundarySourceSimple(AopPim, bPim);
            scalarField bB1Pre(B1Pre.source());
            addBoundarySourceSimple(B1Pre, bB1Pre);
            scalarField bB2Pre(B2Pre.source());
            addBoundarySourceSimple(B2Pre, bB2Pre);
            bPim -= bB1Pre;
            bPim -= bB2Pre;

            scalarField bPre(AopPre.source());
            addBoundarySourceSimple(AopPre, bPre);
            scalarField bB1Pim(B1Pim.source());
            addBoundarySourceSimple(B1Pim, bB1Pim);
            scalarField bB2Pim(B2Pim.source());
            addBoundarySourceSimple(B2Pim, bB2Pim);
            bPre += bB1Pim;
            bPre += bB2Pim;

            forAll(bPim, i)
            {
                const PetscInt gi = globalCells.toGlobal(i);
                VecSetValue(b, gi + 0, (PetscScalar)bPim[i], INSERT_VALUES);
                VecSetValue(b, gi + N, (PetscScalar)bPre[i], INSERT_VALUES);
            }

            MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
            VecAssemblyBegin(b);
            VecAssemblyEnd(b);

            if (!dumpedMatrixStats)
            {
                dumpedMatrixStats = true;

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

            KSPSolve(ksp, b, x);

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

            // Scatter back to fields
            const PetscScalar* xArr = nullptr;
            VecGetArrayRead(x, &xArr);
            forAll(Pim, i)
            {
                const PetscInt gi = globalCells.toGlobal(i);
                Pim[i] = (scalar)PetscRealPart(xArr[gi + 0]);
                Pre[i] = (scalar)PetscRealPart(xArr[gi + N]);
            }
            VecRestoreArrayRead(x, &xArr);

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
