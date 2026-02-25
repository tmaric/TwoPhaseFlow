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
#include <cmath>
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
#include "syncTools.H"

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
    const fvMesh& mesh = op.psi().mesh();
    const auto& bpsi = op.psi().boundaryField();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    List<labelList> neighbGlobals(bpsi.size());

    if (!Pstream::parRun())
    {
        return neighbGlobals;
    }

    const label nInternalFaces = mesh.nInternalFaces();
    const label nBoundaryFaces = mesh.nBoundaryFaces();
    const labelList& faceOwner = mesh.faceOwner();

    // Boundary-face field containing owner-cell global ids.
    labelList neighbByBFace(nBoundaryFaces, -1);
    for (label bFacei = 0; bFacei < nBoundaryFaces; ++bFacei)
    {
        const label facei = nInternalFaces + bFacei;
        neighbByBFace[bFacei] = globalCells.toGlobal(faceOwner[facei]);
    }

    // Swap to neighbour side using OpenFOAM's patch-face mapping.
    syncTools::swapBoundaryFaceList(mesh, neighbByBFace);

    // Slice back per patch.
    forAll(patches, patchi)
    {
        if (!bpsi[patchi].coupled() || !isA<processorPolyPatch>(patches[patchi]))
        {
            continue;
        }

        const polyPatch& pp = patches[patchi];
        const label start = pp.start() - nInternalFaces;
        neighbGlobals[patchi].setSize(pp.size());
        forAll(neighbGlobals[patchi], facei)
        {
            neighbGlobals[patchi][facei] = neighbByBFace[start + facei];
        }
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

    if (Pstream::parRun())
    {
        const auto& bpsi = op.psi().boundaryField();
        const auto& bcoeffs = op.boundaryCoeffs();
        const auto& icoeffs = op.internalCoeffs();
        const lduAddressing& baddr = op.lduAddr();
        const polyBoundaryMesh& patches = op.psi().mesh().boundaryMesh();
        const List<labelList> neighbGlobals
        (
            gatherProcessorPatchNeighbourGlobals(op, globalCells)
        );

        forAll(bpsi, patchi)
        {
            if (!bpsi[patchi].coupled() || !isA<processorPolyPatch>(patches[patchi]))
            {
                continue;
            }

            const labelUList& paddr = baddr.patchAddr(patchi);
            const Field<scalar>& pbc = bcoeffs[patchi];
            const Field<scalar>& pic = icoeffs[patchi];
            const labelList& nbrIds = neighbGlobals[patchi];

            forAll(paddr, facei)
            {
                const label own = paddr[facei];
                const PetscInt r = globalCells.toGlobal(own) + rowOffset;
                const PetscInt cNbr = nbrIds[facei] + colOffset;

                // Convert processor-face contributions to fully implicit form.
                diagWithProc[own] += -pic[facei];
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

static void reportMatrixActionChecksum
(
    Mat& M,
    const globalIndex& globalCells,
    const PetscInt nLocal,
    const PetscInt N
)
{
    Vec x, y;
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, 2*nLocal, 2*N);
    VecSetFromOptions(x);
    VecDuplicate(x, &y);

    for (label i = 0; i < nLocal; ++i)
    {
        const PetscInt gi = globalCells.toGlobal(i);
        const PetscScalar v0 = std::sin(0.001*(gi + 1)) + 0.1*std::cos(0.013*(gi + 3));
        const PetscScalar v1 = std::cos(0.002*(gi + 2)) - 0.05*std::sin(0.017*(gi + 5));
        VecSetValue(x, gi + 0, v0, INSERT_VALUES);
        VecSetValue(x, gi + N, v1, INSERT_VALUES);
    }
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);

    MatMult(M, x, y);

    PetscReal yNorm2 = 0.0;
    VecNorm(y, NORM_2, &yNorm2);

    VecAbs(y);
    PetscReal yNorm1 = 0.0;
    VecNorm(y, NORM_1, &yNorm1);

    Info<< "Matrix action checksum: ||A*x||2=" << yNorm2
        << " ||A*x||1(abs)=" << yNorm1 << nl;

    VecDestroy(&x);
    VecDestroy(&y);
}

static void reportCouplingEntrySample
(
    Mat& M,
    const PetscInt N
)
{
    PetscInt rStart = 0, rEnd = 0;
    MatGetOwnershipRange(M, &rStart, &rEnd);

    bool printed = false;
    for (PetscInt r = rStart; r < rEnd && !printed; ++r)
    {
        if (r >= N) { break; }

        PetscInt ncols = 0;
        const PetscInt* cols = nullptr;
        const PetscScalar* vals = nullptr;
        MatGetRow(M, r, &ncols, &cols, &vals);
        for (PetscInt k = 0; k < ncols; ++k)
        {
            const PetscInt c = cols[k];
            if (c < N || c >= 2*N) { continue; }

            PetscScalar a01 = vals[k];
            if (mag(PetscRealPart(a01)) < 1e-18) { continue; }

            Pout<< "Coupling sample proc=" << Pstream::myProcNo()
                << " row=" << r
                << " col=" << c
                << " A01=" << PetscRealPart(a01) << nl;
            printed = true;
            break;
        }
        MatRestoreRow(M, r, &ncols, &cols, &vals);
    }

    if (!printed)
    {
        Pout<< "Coupling sample proc=" << Pstream::myProcNo()
            << " none in owned rows" << nl;
    }
}

static void reportRequestedCouplingEntries(Mat& M)
{
    const PetscInt rows[2] = {58717, 95250};
    const PetscInt cols[2] = {154072, 190494};
    PetscInt rStart = 0, rEnd = 0;
    MatGetOwnershipRange(M, &rStart, &rEnd);

    for (int i = 0; i < 2; ++i)
    {
        if (rows[i] < rStart || rows[i] >= rEnd) { continue; }
        PetscScalar a = 0.0;
        MatGetValues(M, 1, &rows[i], 1, &cols[i], &a);
        Pout<< "Requested entry proc=" << Pstream::myProcNo()
            << " A(" << rows[i] << "," << cols[i] << ")="
            << PetscRealPart(a) << nl;
    }
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

            // Pre-row couples to Pim via (laplacian(TC1,Pim) + Sp(SC1,Pim))
            fvScalarMatrix couplingLaplPim(fvm::laplacian(TC1, Pim));
            fvScalarMatrix couplingMassPim(fvm::Sp(SC1, Pim));

            if (!dumpedOpStats)
            {
                reportFvMatrixBreakdown("AopPre beforeBoundaryManipulate", AopPre);
                reportFvMatrixBreakdown("AopPim beforeBoundaryManipulate", AopPim);
                reportFvMatrixBreakdown("B1Pre beforeBoundaryManipulate", couplingLaplPre);
                reportFvMatrixBreakdown("B2Pre beforeBoundaryManipulate", couplingMassPre);
            }

            assembleBlockSystem
            (
                M, globalCells, N,
                AopPim, AopPre,
                couplingLaplPre, couplingMassPre,
                couplingLaplPim, couplingMassPim
            );

            if (!dumpedOpStats)
            {
                reportFvMatrixBreakdown("AopPre afterBoundaryManipulate", AopPre);
                reportFvMatrixBreakdown("AopPim afterBoundaryManipulate", AopPim);
                reportFvMatrixBreakdown("B1Pre afterBoundaryManipulate", couplingLaplPre);
                reportFvMatrixBreakdown("B2Pre afterBoundaryManipulate", couplingMassPre);
                dumpedOpStats = true;
            }

            // RHS from source terms (includes BC contributions)
            scalarField bPim(sourceWithBoundary(AopPim, false));
            scalarField bCouplingLaplPre(sourceWithBoundary(couplingLaplPre, false));
            scalarField bCouplingMassPre(sourceWithBoundary(couplingMassPre, false));
            bPim -= bCouplingLaplPre;
            bPim -= bCouplingMassPre;

            scalarField bPre(sourceWithBoundary(AopPre, false));
            scalarField bCouplingLaplPim(sourceWithBoundary(couplingLaplPim, false));
            scalarField bCouplingMassPim(sourceWithBoundary(couplingMassPim, false));
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
                reportMatrixActionChecksum(M, globalCells, nLocal, N);
                reportCouplingEntrySample(M, N);
                reportRequestedCouplingEntries(M);
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
