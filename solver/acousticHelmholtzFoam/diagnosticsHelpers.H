// Diagnostics helpers for acousticHelmholtzFoam.
//
// reportFvMatrixBreakdown(): operator-level coefficient balance check.
// reportMatrixStats(): global PETSc diagnostics + serial/MPI fingerprints.
// reportKspStats(): linear solve convergence report.

#include <cstdlib>
#include <fstream>

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
    if (!Pstream::parRun())
    {
        IS is0 = nullptr, is1 = nullptr;
        Mat M01 = nullptr, M10 = nullptr, M01T = nullptr, Mtmp = nullptr;
        ISCreateStride(PETSC_COMM_WORLD, N, 0, 1, &is0);
        ISCreateStride(PETSC_COMM_WORLD, N, N, 1, &is1);
        MatCreateSubMatrix(M, is0, is1, MAT_INITIAL_MATRIX, &M01);
        MatCreateSubMatrix(M, is1, is0, MAT_INITIAL_MATRIX, &M10);
        MatNorm(M01, NORM_FROBENIUS, &n01);
        MatNorm(M10, NORM_FROBENIUS, &n10);
        MatTranspose(M01, MAT_INITIAL_MATRIX, &M01T);

        MatDuplicate(M10, MAT_COPY_VALUES, &Mtmp);
        MatAXPY(Mtmp, -1.0, M01T, DIFFERENT_NONZERO_PATTERN);
        MatNorm(Mtmp, NORM_FROBENIUS, &nSym);
        MatDestroy(&Mtmp);

        MatDuplicate(M10, MAT_COPY_VALUES, &Mtmp);
        MatAXPY(Mtmp, +1.0, M01T, DIFFERENT_NONZERO_PATTERN);
        MatNorm(Mtmp, NORM_FROBENIUS, &nSkew);
        MatDestroy(&Mtmp);

        MatDestroy(&M01);
        MatDestroy(&M10);
        MatDestroy(&M01T);
        ISDestroy(&is0);
        ISDestroy(&is1);
    }

    // Entry-based fingerprints are robust in MPI and serial, unlike direct
    // off-block submatrix extraction on distributed matrices.
    PetscInt rStart = 0, rEnd = 0;
    MatGetOwnershipRange(M, &rStart, &rEnd);
    scalar sumVal = 0.0, sumAbs = 0.0, sumSq = 0.0;
    scalar h1 = 0.0, h2 = 0.0;
    scalar sumVal01 = 0.0, sumVal10 = 0.0;
    scalar sumAbs01 = 0.0, sumAbs10 = 0.0;
    scalar sumSq01 = 0.0, sumSq10 = 0.0;
    scalar h1_01 = 0.0, h1_10 = 0.0;

    for (PetscInt r = rStart; r < rEnd; ++r)
    {
        PetscInt ncols = 0;
        const PetscInt* cols = nullptr;
        const PetscScalar* vals = nullptr;
        MatGetRow(M, r, &ncols, &cols, &vals);
        for (PetscInt k = 0; k < ncols; ++k)
        {
            const PetscInt c = cols[k];
            const scalar v = PetscRealPart(vals[k]);
            const scalar av = mag(v);
            sumVal += v;
            sumAbs += av;
            sumSq += v*v;
            h1 += v*(1.0 + 1e-6*(r + 1) + 1e-12*(c + 1));
            h2 += v*(0.5*(r + 1) - 0.25*(c + 1));

            if (r < N && c >= N)
            {
                sumVal01 += v;
                sumAbs01 += av;
                sumSq01 += v*v;
                h1_01 += v*(1.0 + 1e-6*(r + 1) + 1e-12*(c - N + 1));
            }
            else if (r >= N && c < N)
            {
                sumVal10 += v;
                sumAbs10 += av;
                sumSq10 += v*v;
                h1_10 += v*(1.0 + 1e-6*(r - N + 1) + 1e-12*(c + 1));
            }
        }
        MatRestoreRow(M, r, &ncols, &cols, &vals);
    }

    if (Pstream::parRun())
    {
        reduce(sumVal, sumOp<scalar>());
        reduce(sumAbs, sumOp<scalar>());
        reduce(sumSq, sumOp<scalar>());
        reduce(h1, sumOp<scalar>());
        reduce(h2, sumOp<scalar>());
        reduce(sumVal01, sumOp<scalar>());
        reduce(sumAbs01, sumOp<scalar>());
        reduce(sumSq01, sumOp<scalar>());
        reduce(h1_01, sumOp<scalar>());
        reduce(sumVal10, sumOp<scalar>());
        reduce(sumAbs10, sumOp<scalar>());
        reduce(sumSq10, sumOp<scalar>());
        reduce(h1_10, sumOp<scalar>());
    }

    const scalar n01Direct = Foam::sqrt(max(sumSq01, scalar(0)));
    const scalar n10Direct = Foam::sqrt(max(sumSq10, scalar(0)));

    Info<< "Matrix stats: ||A||1=" << n1
        << " ||A||inf=" << ninf
        << " ||A||F=" << nfrob
        << " ||A01||F=" << n01
        << " ||A10||F=" << n10
        << " ||A10-A01T||F=" << nSym
        << " ||A10+A01T||F=" << nSkew
        << " diag[min,max]=( " << dmin << "," << dmax << ")"
        << " symmetric=" << (isSym ? "yes" : "no")
        << " rhsNorm=" << bNorm0
        << " fp(sum,sumAbs,sumSq)=(" << sumVal << "," << sumAbs << "," << sumSq << ")"
        << " fp(h1,h2)=(" << h1 << "," << h2 << ")"
        << " fp01(sum,abs,h1,frob)=(" << sumVal01 << "," << sumAbs01 << "," << h1_01 << "," << n01Direct << ")"
        << " fp10(sum,abs,h1,frob)=(" << sumVal10 << "," << sumAbs10 << "," << h1_10 << "," << n10Direct << ")"
        << " offDiagBlockSubmat=" << (Pstream::parRun() ? "disabledInParallel" : "enabledInSerial")
        << nl;
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

static bool dumpRowsEnabled()
{
    const char* env = std::getenv("FREBCM_DUMP_ROWS");
    return env && word(env) != "0";
}

static void dumpProcessorInterfaceRows
(
    Mat& M,
    const fvScalarMatrix& op,
    const globalIndex& globalCells,
    const PetscInt N
)
{
    if (!dumpRowsEnabled())
    {
        return;
    }

    const auto& bpsi = op.psi().boundaryField();
    const lduAddressing& addr = op.lduAddr();
    const polyBoundaryMesh& patches = op.psi().mesh().boundaryMesh();

    DynamicList<label> localRows;
    forAll(bpsi, patchi)
    {
        if (!bpsi[patchi].coupled() || !isA<processorPolyPatch>(patches[patchi]))
        {
            continue;
        }

        const labelUList& paddr = addr.patchAddr(patchi);
        forAll(paddr, facei)
        {
            const label celli = paddr[facei];
            if (!localRows.found(celli))
            {
                localRows.append(celli);
            }
            if (localRows.size() >= 6)
            {
                break;
            }
        }
        if (localRows.size() >= 6)
        {
            break;
        }
    }

    if (!localRows.size())
    {
        return;
    }

    const fileName outName
    (
        "matrixRows_rank" + Foam::name(Pstream::myProcNo()) + ".txt"
    );
    OFstream os(outName);
    os << "rank " << Pstream::myProcNo() << " selectedRows " << localRows.size() << nl;

    forAll(localRows, i)
    {
        const PetscInt g = globalCells.toGlobal(localRows[i]);
        for (label blk = 0; blk < 2; ++blk)
        {
            const PetscInt row = (blk == 0) ? g : (g + N);
            PetscInt ncols = 0;
            const PetscInt* cols = nullptr;
            const PetscScalar* vals = nullptr;
            MatGetRow(M, row, &ncols, &cols, &vals);
            os << "row " << row << " ncols " << ncols << nl;
            for (PetscInt k = 0; k < ncols; ++k)
            {
                os << "  " << cols[k] << " " << PetscRealPart(vals[k]) << nl;
            }
            MatRestoreRow(M, row, &ncols, &cols, &vals);
        }
    }
}
