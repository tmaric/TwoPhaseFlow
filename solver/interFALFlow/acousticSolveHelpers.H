// PETSc lifecycle and post-processing helpers for interFALFlow acoustic step.
//
// This follows acousticHelmholtzFoam behaviour: distributed AIJ + direct
// factorization, with rank-local solution extraction after KSPSolve(). The
// assembled matrix, RHS, and block ordering stay unchanged.

struct AcousticPetscSystem
{
    Mat M = nullptr;
    Vec x = nullptr;
    Vec b = nullptr;
    KSP ksp = nullptr;
    IS scatterFromIs = nullptr;
    IS scatterToIs = nullptr;
    Vec xLocal = nullptr;
    VecScatter localScatter = nullptr;
    List<PetscInt> localBlockIds;
    PetscInt nLocal = -1;
    PetscInt N = -1;
    PetscInt rowNnzUpper = -1;
};

static void destroyPetscSystem
(
    AcousticPetscSystem& system
)
{
    if (system.localScatter) VecScatterDestroy(&system.localScatter);
    if (system.xLocal) VecDestroy(&system.xLocal);
    if (system.scatterToIs) ISDestroy(&system.scatterToIs);
    if (system.scatterFromIs) ISDestroy(&system.scatterFromIs);
    if (system.ksp) KSPDestroy(&system.ksp);
    if (system.x) VecDestroy(&system.x);
    if (system.b) VecDestroy(&system.b);
    if (system.M) MatDestroy(&system.M);

    system.localScatter = nullptr;
    system.xLocal = nullptr;
    system.scatterToIs = nullptr;
    system.scatterFromIs = nullptr;
    system.localBlockIds.clear();
    system.ksp = nullptr;
    system.x = nullptr;
    system.b = nullptr;
    system.M = nullptr;
    system.nLocal = -1;
    system.N = -1;
    system.rowNnzUpper = -1;
}

static void initializePetscSystem
(
    const PetscInt nLocal,
    const PetscInt N,
    const PetscInt rowNnzUpper,
    AcousticPetscSystem& system
)
{
    MatCreate(PETSC_COMM_WORLD, &system.M);
    MatSetSizes(system.M, 2*nLocal, 2*nLocal, 2*N, 2*N);
    MatSetType(system.M, MATAIJ);
    if (Pstream::parRun())
    {
        // Rows are inserted by global block index, while PETSc owns rows in a
        // contiguous global range. Use a conservative per-row upper bound so
        // preallocation stays valid without changing the assembled operator.
        MatMPIAIJSetPreallocation
        (
            system.M,
            rowNnzUpper,
            nullptr,
            rowNnzUpper,
            nullptr
        );
    }
    else
    {
        MatSeqAIJSetPreallocation(system.M, rowNnzUpper, nullptr);
    }
    MatSetUp(system.M);

    VecCreate(PETSC_COMM_WORLD, &system.x);
    VecSetSizes(system.x, 2*nLocal, 2*N);
    VecSetFromOptions(system.x);
    VecDuplicate(system.x, &system.b);

    KSPCreate(PETSC_COMM_WORLD, &system.ksp);
    KSPSetOperators(system.ksp, system.M, system.M);

    PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");
    PetscOptionsSetValue(nullptr, "-pc_type", "lu");
    PetscOptionsSetValue(nullptr, "-pc_factor_mat_solver_type", "mumps");
    KSPSetFromOptions(system.ksp);

    system.nLocal = nLocal;
    system.N = N;
    system.rowNnzUpper = rowNnzUpper;
}

static void ensurePetscSystem
(
    const PetscInt nLocal,
    const PetscInt N,
    const PetscInt rowNnzUpper,
    AcousticPetscSystem& system
)
{
    const bool needsRebuild =
        !system.M
     || !system.x
     || !system.b
     || !system.ksp
     || system.nLocal != nLocal
     || system.N != N
     || system.rowNnzUpper != rowNnzUpper;

    if (needsRebuild)
    {
        destroyPetscSystem(system);
        initializePetscSystem(nLocal, N, rowNnzUpper, system);
    }
}

static void resetPetscSystem
(
    AcousticPetscSystem& system
)
{
    MatZeroEntries(system.M);
    VecSet(system.x, 0.0);
    VecSet(system.b, 0.0);
}

static void ensureSolutionScatter
(
    const globalIndex& globalCells,
    AcousticPetscSystem& system
)
{
    List<PetscInt> currentIds(2*system.nLocal);
    for (PetscInt i = 0; i < system.nLocal; ++i)
    {
        const PetscInt gi = globalCells.toGlobal(i);
        currentIds[i] = gi;
        currentIds[i + system.nLocal] = gi + system.N;
    }

    bool needsRebuild =
        !system.localScatter
     || system.localBlockIds.size() != currentIds.size();

    if (!needsRebuild)
    {
        forAll(currentIds, i)
        {
            if (system.localBlockIds[i] != currentIds[i])
            {
                needsRebuild = true;
                break;
            }
        }
    }

    if (!needsRebuild)
    {
        return;
    }

    if (system.localScatter) VecScatterDestroy(&system.localScatter);
    if (system.xLocal) VecDestroy(&system.xLocal);
    if (system.scatterToIs) ISDestroy(&system.scatterToIs);
    if (system.scatterFromIs) ISDestroy(&system.scatterFromIs);

    system.localBlockIds.transfer(currentIds);

    ISCreateGeneral
    (
        PETSC_COMM_SELF,
        system.localBlockIds.size(),
        system.localBlockIds.begin(),
        PETSC_COPY_VALUES,
        &system.scatterFromIs
    );
    ISCreateStride
    (
        PETSC_COMM_SELF,
        2*system.nLocal,
        0,
        1,
        &system.scatterToIs
    );

    VecCreateSeq(PETSC_COMM_SELF, 2*system.nLocal, &system.xLocal);
    VecScatterCreate
    (
        system.x,
        system.scatterFromIs,
        system.xLocal,
        system.scatterToIs,
        &system.localScatter
    );
}

static void scatterBlockSolution
(
    AcousticPetscSystem& system,
    const globalIndex& globalCells,
    volScalarField& Pim,
    volScalarField& Pre
)
{
    ensureSolutionScatter(globalCells, system);

    VecScatterBegin
    (
        system.localScatter,
        system.x,
        system.xLocal,
        INSERT_VALUES,
        SCATTER_FORWARD
    );
    VecScatterEnd
    (
        system.localScatter,
        system.x,
        system.xLocal,
        INSERT_VALUES,
        SCATTER_FORWARD
    );

    const PetscScalar* xArr = nullptr;
    VecGetArrayRead(system.xLocal, &xArr);
    forAll(Pim, i)
    {
        Pim[i] = (scalar)PetscRealPart(xArr[i]);
        Pre[i] = (scalar)PetscRealPart(xArr[i + system.nLocal]);
    }
    VecRestoreArrayRead(system.xLocal, &xArr);
}

static PetscInt estimateBlockPreallocationUpper
(
    const fvScalarMatrix& op
)
{
    const label nCells = op.psi().mesh().nCells();
    labelList localNbrCount(nCells, 0);
    labelList remoteNbrCount(nCells, 0);

    const lduAddressing& addr = op.lduAddr();
    const labelUList& owner = addr.lowerAddr();
    const labelUList& neigh = addr.upperAddr();

    forAll(owner, facei)
    {
        ++localNbrCount[owner[facei]];
        ++localNbrCount[neigh[facei]];
    }

    const auto& bpsi = op.psi().boundaryField();
    const polyBoundaryMesh& patches = op.psi().mesh().boundaryMesh();

    forAll(bpsi, patchi)
    {
        if
        (
            bpsi[patchi].coupled()
         && isA<processorPolyPatch>(patches[patchi])
        )
        {
            const labelUList& paddr = addr.patchAddr(patchi);
            forAll(paddr, facei)
            {
                ++remoteNbrCount[paddr[facei]];
            }
        }
    }

    PetscInt rowNnzUpper = 0;
    for (label celli = 0; celli < nCells; ++celli)
    {
        const PetscInt rowNnz =
            static_cast<PetscInt>
            (
                2*(1 + localNbrCount[celli] + remoteNbrCount[celli])
            );
        rowNnzUpper = max(rowNnzUpper, rowNnz);
    }

    if (Pstream::parRun())
    {
        reduce(rowNnzUpper, maxOp<PetscInt>());
    }

    return rowNnzUpper;
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
    volScalarField& Pa,
    volScalarField& Ua,
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
    const volScalarField pSq("pSq", Pre*Pre + Pim*Pim);
    const volScalarField uSq("uSq", (Ure() & Ure()) + (Uim() & Uim()));

    Pa == Foam::sqrt(pSq);
    Ua == Foam::sqrt(uSq);

    pr == 0.25*(kl*alpha1 + kg*(1-alpha1))*pSq
        - 0.25*rho*uSq;
    momFlux == 0.5*rho*(Ure()*Ure() + Uim()*Uim());
}
