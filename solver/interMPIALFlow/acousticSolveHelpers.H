// PETSc lifecycle and post-processing helpers for interMPIALFlow acoustic step.
//
// This follows freBCMFoam behaviour: distributed AIJ + direct factorization,
// gather solved block vector to all ranks, and construct RHS with processor
// coupled source terms excluded to avoid double counting.

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

    PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");
    PetscOptionsSetValue(nullptr, "-pc_type", "lu");
    PetscOptionsSetValue(nullptr, "-pc_factor_mat_solver_type", "mumps");
    KSPSetFromOptions(ksp);
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

static void buildRhs
(
    const fvScalarMatrix& AopPim,
    const fvScalarMatrix& AopPre,
    const fvScalarMatrix& B1Pre,
    const fvScalarMatrix& B2Pre,
    const fvScalarMatrix& B1Pim,
    const fvScalarMatrix& B2Pim,
    scalarField& bPim,
    scalarField& bPre
)
{
    bPim = sourceWithBoundary(AopPim, true, false);
    bPim -= sourceWithBoundary(B1Pre, true, false);
    bPim -= sourceWithBoundary(B2Pre, true, false);

    bPre = sourceWithBoundary(AopPre, true, false);
    bPre += sourceWithBoundary(B1Pim, true, false);
    bPre += sourceWithBoundary(B2Pim, true, false);
}

static void updateDerivedAcousticFields
(
    volScalarField& pr,
    volTensorField& momFlux,
    const volScalarField& Pim,
    const volScalarField& Pre,
    const volScalarField& alpha1,
    const volScalarField& rho,
    const dimensionedScalar& f,
    const dimensionedScalar& kl,
    const dimensionedScalar& kg
)
{
    const tmp<volVectorField> Ure = 1/(2*constant::mathematical::pi*f*rho) * fvc::grad(Pim);
    const tmp<volVectorField> Uim = -1/(2*constant::mathematical::pi*f*rho) * fvc::grad(Pre);

    pr == 0.25*(kl*alpha1 + kg*(1-alpha1))*(Pre*Pre + Pim*Pim)
        - 0.25*rho*((Ure()&Ure()) + (Uim()&Uim()));
    momFlux == 0.5*rho*(Ure()*Ure() + Uim()*Uim());
}
