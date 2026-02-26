// PETSc solve helpers for freBCMFoam.
//
// After distributed KSPSolve, scatterBlockSolution() gathers the full block
// vector on each rank to keep OpenFOAM field updates straightforward.
// buildRhs() and updateDerivedAcousticFields() mirror the serial reference
// formulation to ease one-to-one comparisons.

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
    const fvScalarMatrix& mainOpPim,
    const fvScalarMatrix& mainOpPre,
    const fvScalarMatrix& couplingOpPre,
    const fvScalarMatrix& couplingOpPim,
    scalarField& bPim,
    scalarField& bPre
)
{
    bPim = sourceWithBoundary(mainOpPim, false);
    bPim -= sourceWithBoundary(couplingOpPre, false);

    bPre = sourceWithBoundary(mainOpPre, false);
    bPre += sourceWithBoundary(couplingOpPim, false);
}

static void updateDerivedAcousticFields
(
    volVectorField& Ure,
    volVectorField& Uim,
    volScalarField& pa,
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
    Ure == 1/(2*constant::mathematical::pi*f*rho) * fvc::grad(Pim);
    Uim == -1/(2*constant::mathematical::pi*f*rho) * fvc::grad(Pre);
    pa == Foam::sqrt(Pim*Pim + Pre*Pre);
    pr == 0.25*(kl*alpha1 + kg*(1-alpha1))*(Pre*Pre + Pim*Pim)
        - 0.25*rho*((Ure&Ure) + (Uim&Uim));
    momFlux == 0.5*rho*(Ure*Ure + Uim*Uim);
}
