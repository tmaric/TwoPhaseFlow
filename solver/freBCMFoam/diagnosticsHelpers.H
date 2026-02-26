// Diagnostics helpers for freBCMFoam.
//
// reportFvMatrixBreakdown(): local/operator-level sanity check before PETSc.
// reportMatrixStats(): global PETSc matrix check after assembly.
// reportKspStats(): linear solve convergence summary each iteration.

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

    Info<< "Matrix stats: ||A||1=" << n1
        << " ||A||inf=" << ninf
        << " ||A||F=" << nfrob
        << " ||A01||F=" << n01
        << " ||A10||F=" << n10
        << " ||A10-A01T||F=" << nSym
        << " ||A10+A01T||F=" << nSkew
        << " diag[min,max]=( " << dmin << "," << dmax << ")"
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
