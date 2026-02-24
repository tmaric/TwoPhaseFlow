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
    freBCMFoam

Group
    AcousticSolvers

Description
    Block-coupled frequency-domain acoustic solver (Pre/Pim) using PETSc.
    Rectangle/spherical PML with scalar damping coefficients.
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

static scalar sumField(const scalarField& f)
{
    scalar s = 0.0;
    forAll(f, i) { s += f[i]; }
    return s;
}

static scalar sumMagField(const scalarField& f)
{
    scalar s = 0.0;
    forAll(f, i) { s += mag(f[i]); }
    return s;
}

static scalar globalSum(const scalar local)
{
    scalar g = local;
    if (Pstream::parRun())
    {
        reduce(g, sumOp<scalar>());
    }
    return g;
}

static void reportFvMatrixBreakdown
(
    const word& name,
    const fvScalarMatrix& op
)
{
    const lduMatrix& L = op;
    const auto& bpsi = op.psi().boundaryField();
    const auto& bcoeffs = op.boundaryCoeffs();
    const auto& icoeffs = op.internalCoeffs();

    scalar sDiag = globalSum(sumField(L.diag()));
    scalar aDiag = globalSum(sumMagField(L.diag()));
    scalar sUpper = 0.0, aUpper = 0.0, sLower = 0.0, aLower = 0.0;

    if (L.hasUpper())
    {
        sUpper = globalSum(sumField(L.upper()));
        aUpper = globalSum(sumMagField(L.upper()));
    }
    if (L.hasLower())
    {
        sLower = globalSum(sumField(L.lower()));
        aLower = globalSum(sumMagField(L.lower()));
    }

    scalar bCoupled = 0.0, bUncoupled = 0.0, iCoupled = 0.0, iUncoupled = 0.0;
    label nCoupledFaces = 0, nUncoupledFaces = 0;
    forAll(bpsi, patchi)
    {
        const bool coupled = bpsi[patchi].coupled();
        const scalar sb = sumMagField(bcoeffs[patchi]);
        const scalar si = sumMagField(icoeffs[patchi]);
        if (coupled)
        {
            bCoupled += sb;
            iCoupled += si;
            nCoupledFaces += bcoeffs[patchi].size();
        }
        else
        {
            bUncoupled += sb;
            iUncoupled += si;
            nUncoupledFaces += bcoeffs[patchi].size();
        }
    }
    bCoupled = globalSum(bCoupled);
    bUncoupled = globalSum(bUncoupled);
    iCoupled = globalSum(iCoupled);
    iUncoupled = globalSum(iUncoupled);
    reduce(nCoupledFaces, sumOp<label>());
    reduce(nUncoupledFaces, sumOp<label>());

    if (Pstream::master())
    {
        Info<< "FV matrix [" << name << "] "
            << "diag(sum,abs)=(" << sDiag << "," << aDiag << ") "
            << "upper(sum,abs)=(" << sUpper << "," << aUpper << ") "
            << "lower(sum,abs)=(" << sLower << "," << aLower << ") "
            << "bCoeff|coupled=" << bCoupled
            << " bCoeff|uncoupled=" << bUncoupled
            << " iCoeff|coupled=" << iCoupled
            << " iCoeff|uncoupled=" << iUncoupled
            << " nCoupledFaces=" << nCoupledFaces
            << " nUncoupledFaces=" << nUncoupledFaces
            << nl;
    }
}

static List<labelList> gatherProcessorPatchNeighbourGlobals
(
    const fvScalarMatrix& op,
    const globalIndex& globalCells
)
{
    const auto& bpsi = op.psi().boundaryField();
    const lduAddressing& addr = op.lduAddr();
    const polyBoundaryMesh& patches = op.psi().mesh().boundaryMesh();

    List<labelList> neighbGlobals(bpsi.size());

    if (!Pstream::parRun())
    {
        return neighbGlobals;
    }

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Send local owner-cell global ids on each processor patch.
    forAll(bpsi, patchi)
    {
        if (!bpsi[patchi].coupled() || !isA<processorPolyPatch>(patches[patchi]))
        {
            continue;
        }

        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(patches[patchi]);
        const labelUList& paddr = addr.patchAddr(patchi);

        labelList sendIds(paddr.size());
        forAll(paddr, facei)
        {
            sendIds[facei] = globalCells.toGlobal(paddr[facei]);
        }

        UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
        toNbr << sendIds;
    }

    pBufs.finishedSends();

    // Receive neighbour owner-cell global ids (same face ordering).
    forAll(bpsi, patchi)
    {
        if (!bpsi[patchi].coupled() || !isA<processorPolyPatch>(patches[patchi]))
        {
            continue;
        }

        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(patches[patchi]);
        UIPstream fromNbr(procPatch.neighbProcNo(), pBufs);
        fromNbr >> neighbGlobals[patchi];
    }

    return neighbGlobals;
}

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
    scalarField diagWithProc(diag);

    const lduAddressing& addr = L.lduAddr();
    const labelUList& owner = addr.lowerAddr();
    const labelUList& neigh = addr.upperAddr();

    // Processor-patch couplings are not in lower/upper arrays in parallel.
    // Rebuild their implicit contributions explicitly.
    List<labelList> neighbGlobals;
    if (Pstream::parRun() && rowOffset == colOffset)
    {
        const auto& bpsi = op.psi().boundaryField();
        const auto& bcoeffs = op.boundaryCoeffs();
        const auto& icoeffs = op.internalCoeffs();
        const lduAddressing& baddr = op.lduAddr();
        const polyBoundaryMesh& patches = op.psi().mesh().boundaryMesh();
        neighbGlobals = gatherProcessorPatchNeighbourGlobals(op, globalCells);

        forAll(bpsi, patchi)
        {
            if (!bpsi[patchi].coupled() || !isA<processorPolyPatch>(patches[patchi]))
            {
                continue;
            }

            const labelUList& paddr = baddr.patchAddr(patchi);
            const Field<scalar>& pbc = bcoeffs[patchi];
            const labelList& nbrIds = neighbGlobals[patchi];

            if (nbrIds.size() != paddr.size())
            {
                FatalErrorInFunction
                    << "Processor patch exchange size mismatch on patch "
                    << patches[patchi].name()
                    << ": local faces " << paddr.size()
                    << " neighbour ids " << nbrIds.size()
                    << exit(FatalError);
            }

            const Field<scalar>& pic = icoeffs[patchi];

            forAll(paddr, facei)
            {
                const PetscInt cNbr = nbrIds[facei] + colOffset;
                const label ownCell = paddr[facei];
                const PetscInt r = globalCells.toGlobal(ownCell) + rowOffset;

                // Rebuild processor-face implicit stencil in matrix form:
                // diagonal contribution from internalCoeff and neighbour
                // coupling from boundaryCoeff.
                diagWithProc[ownCell] += pic[facei];
                MatSetValue(M, r, cNbr, -scale*pbc[facei], INSERT_VALUES);
            }
        }
    }

    // Diagonal
    forAll(diagWithProc, i)
    {
        PetscInt r = globalCells.toGlobal(i) + rowOffset;
        PetscInt c = globalCells.toGlobal(i) + colOffset;
        MatSetValue(M, r, c, scale*diagWithProc[i], INSERT_VALUES);
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

// Insert transpose(op) into one PETSc block position.
static void insertScalarOpTransposeIntoBlock
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

    forAll(diag, i)
    {
        const PetscInt r = globalCells.toGlobal(i) + rowOffset;
        const PetscInt c = globalCells.toGlobal(i) + colOffset;
        MatSetValue(M, r, c, scale*diag[i], INSERT_VALUES);
    }

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

            // Transpose swaps upper/lower contributions.
            MatSetValue(M, ri, cj, scale*lowerVal, INSERT_VALUES);
            MatSetValue(M, rj, ci, scale*upperVal, INSERT_VALUES);
        }
    }
}

static scalarField sourceWithBoundary
(
    const fvScalarMatrix& M,
    const bool includeCoupled
)
{
    scalarField source(M.source());
    addBoundarySourceSimple(M, source, includeCoupled);
    return source;
}

static void assembleBlockSystem
(
    Mat& M,
    const globalIndex& globalCells,
    const PetscInt N,
    fvScalarMatrix& AopPim,
    fvScalarMatrix& AopPre,
    fvScalarMatrix& BPre
)
{
    // Apply boundary contributions to matrices before insertion.
    AopPre.boundaryManipulate(AopPre.psi().boundaryFieldRef());
    AopPim.boundaryManipulate(AopPim.psi().boundaryFieldRef());
    BPre.boundaryManipulate(BPre.psi().boundaryFieldRef());

    // Coupled system:
    // [ A  -(B1 + B2) ] [Pim] = [bPim]
    // [ (B1 + B2)  A ] [Pre]  [bPre]
    insertScalarOpIntoBlock(M, AopPim, globalCells, 0, 0, +1.0);
    insertScalarOpIntoBlock(M, BPre, globalCells, 0, N, -1.0);

    // Enforce skew block structure from one coupling operator:
    // A10 = -A01^T, with A01 = -BPre.
    insertScalarOpTransposeIntoBlock(M, BPre, globalCells, N, 0, +1.0);
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
    MatSetType(M, MATAIJ);
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
    Vec& b,
    const PetscInt N
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

    PetscReal n01 = 0.0, n10 = 0.0, nSkew = 0.0, nSym = 0.0;
    IS is0 = nullptr, is1 = nullptr;
    Mat M01 = nullptr, M10 = nullptr, M01T = nullptr, Mtmp = nullptr;
    ISCreateStride(PETSC_COMM_WORLD, N, 0, 1, &is0);
    ISCreateStride(PETSC_COMM_WORLD, N, N, 1, &is1);
    MatCreateSubMatrix(M, is0, is1, MAT_INITIAL_MATRIX, &M01);
    MatCreateSubMatrix(M, is1, is0, MAT_INITIAL_MATRIX, &M10);
    MatNorm(M01, NORM_FROBENIUS, &n01);
    MatNorm(M10, NORM_FROBENIUS, &n10);
    MatTranspose(M01, MAT_INITIAL_MATRIX, &M01T);

    // Symmetry/skew diagnostics for off-diagonal blocks.
    MatDuplicate(M10, MAT_COPY_VALUES, &Mtmp);
    MatAXPY(Mtmp, -1.0, M01T, DIFFERENT_NONZERO_PATTERN); // M10 - M01^T
    MatNorm(Mtmp, NORM_FROBENIUS, &nSym);
    MatDestroy(&Mtmp);

    MatDuplicate(M10, MAT_COPY_VALUES, &Mtmp);
    MatAXPY(Mtmp, +1.0, M01T, DIFFERENT_NONZERO_PATTERN); // M10 + M01^T
    MatNorm(Mtmp, NORM_FROBENIUS, &nSkew);
    MatDestroy(&Mtmp);

    MatDestroy(&M01);
    MatDestroy(&M10);
    MatDestroy(&M01T);
    ISDestroy(&is0);
    ISDestroy(&is1);

    Info<< "Matrix stats: ||A||1=" << n1
        << " ||A||inf=" << ninf
        << " ||A||F=" << nfrob
        << " ||A01||F=" << n01
        << " ||A10||F=" << n10
        << " ||A10-A01T||F=" << nSym
        << " ||A10+A01T||F=" << nSkew
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
    Vec xAll = nullptr;
    VecScatter allScatter = nullptr;
    VecScatterCreateToAll(x, &allScatter, &xAll);
    VecScatterBegin(allScatter, x, xAll, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(allScatter, x, xAll, INSERT_VALUES, SCATTER_FORWARD);

    const PetscScalar* xArr = nullptr;
    VecGetArrayRead(xAll, &xArr);
    forAll(Pim, i)
    {
        const PetscInt gi = globalCells.toGlobal(i);
        Pim[i] = (scalar)PetscRealPart(xArr[gi + 0]);
        Pre[i] = (scalar)PetscRealPart(xArr[gi + N]);
    }
    VecRestoreArrayRead(xAll, &xArr);
    VecScatterDestroy(&allScatter);
    VecDestroy(&xAll);
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
    bool dumpedOpStats = false;

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
            fvScalarMatrix couplingPre(couplingLaplPre + couplingMassPre);

            // Pre-row couples to Pim via (laplacian(TC1,Pim) + Sp(SC1,Pim))
            fvScalarMatrix couplingLaplPim(fvm::laplacian(TC1, Pim));
            fvScalarMatrix couplingMassPim(fvm::Sp(SC1, Pim));
            fvScalarMatrix couplingPim(couplingLaplPim + couplingMassPim);

            if (!dumpedOpStats)
            {
                reportFvMatrixBreakdown("AopPre beforeBoundaryManipulate", AopPre);
                reportFvMatrixBreakdown("AopPim beforeBoundaryManipulate", AopPim);
                reportFvMatrixBreakdown("BPre beforeBoundaryManipulate", couplingPre);
            }

            assembleBlockSystem
            (
                M, globalCells, N,
                AopPim, AopPre,
                couplingPre
            );

            if (!dumpedOpStats)
            {
                reportFvMatrixBreakdown("AopPre afterBoundaryManipulate", AopPre);
                reportFvMatrixBreakdown("AopPim afterBoundaryManipulate", AopPim);
                reportFvMatrixBreakdown("BPre afterBoundaryManipulate", couplingPre);
                dumpedOpStats = true;
            }

            // RHS from source terms (includes BC contributions)
            // Keep processor-coupled interfaces implicit in matrix so the
            // decomposed equation matches the undecomposed serial form.
            scalarField bPim(sourceWithBoundary(AopPim, false));
            scalarField bCouplingPre(sourceWithBoundary(couplingPre, false));
            bPim -= bCouplingPre;

            scalarField bPre(sourceWithBoundary(AopPre, false));
            scalarField bCouplingPim(sourceWithBoundary(couplingPim, false));
            bPre += bCouplingPim;

            setBlockRhs(b, globalCells, N, bPim, bPre);

            MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
            VecAssemblyBegin(b);
            VecAssemblyEnd(b);

            if (!dumpedMatrixStats)
            {
                dumpedMatrixStats = true;
                reportMatrixStats(M, b, N);
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
