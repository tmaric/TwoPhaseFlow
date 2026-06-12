// PETSc solve helpers for acousticHelmholtzFoam.
//
// Interleaved cell ordering keeps each rank's OpenFOAM cells in its PETSc
// ownership range:
//   cell g -> Pim row 2*g, Pre row 2*g + 1.
// buildRhs() mirrors the serial reference formulation to ease one-to-one
// comparisons.

static void setDefaultPetscOption(const char* name, const char* value)
{
    PetscBool optionSet = PETSC_FALSE;
    PetscOptionsHasName(nullptr, nullptr, name, &optionSet);
    if (!optionSet)
    {
        PetscOptionsSetValue(nullptr, name, value);
    }
}

static void buildBlockPreallocation
(
    const fvMesh& mesh,
    List<PetscInt>& dNnz,
    List<PetscInt>& oNnz
)
{
    dNnz.setSize(2*mesh.nCells(), 2);
    oNnz.setSize(2*mesh.nCells(), 0);

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    forAll(neighbour, facei)
    {
        for (label component = 0; component < 2; ++component)
        {
            dNnz[2*owner[facei] + component] += 2;
            dNnz[2*neighbour[facei] + component] += 2;
        }
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    forAll(patches, patchi)
    {
        if (!isA<processorPolyPatch>(patches[patchi]))
        {
            continue;
        }

        const labelUList& faceCells = patches[patchi].faceCells();
        forAll(faceCells, facei)
        {
            for (label component = 0; component < 2; ++component)
            {
                oNnz[2*faceCells[facei] + component] += 2;
            }
        }
    }
}

static void initializePetscSystem
(
    const fvMesh& mesh,
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

    List<PetscInt> dNnz;
    List<PetscInt> oNnz;
    buildBlockPreallocation(mesh, dNnz, oNnz);
    if (Pstream::parRun())
    {
        MatMPIAIJSetPreallocation(M, 0, dNnz.begin(), 0, oNnz.begin());
    }
    else
    {
        MatSeqAIJSetPreallocation(M, 0, dNnz.begin());
    }
    MatSetUp(M);
    MatSetOption(M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, 2*nLocal, 2*N);
    VecSetFromOptions(x);
    VecDuplicate(x, &b);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, M, M);

    setDefaultPetscOption("-ksp_type", "preonly");
    setDefaultPetscOption("-pc_type", "lu");
    setDefaultPetscOption("-pc_factor_mat_solver_type", "mumps");
    KSPSetFromOptions(ksp);
}

static void scatterBlockSolution
(
    Vec& x,
    const globalIndex& globalCells,
    volScalarField& Pim,
    volScalarField& Pre
)
{
    PetscInt ownershipStart = 0;
    PetscInt ownershipEnd = 0;
    VecGetOwnershipRange(x, &ownershipStart, &ownershipEnd);

    const PetscInt expectedStart =
        Pim.size() ? 2*(PetscInt)globalCells.toGlobal(0) : ownershipStart;
    if
    (
        ownershipStart != expectedStart
     || ownershipEnd - ownershipStart != 2*(PetscInt)Pim.size()
    )
    {
        FatalErrorInFunction
            << "PETSc/OpenFOAM ownership mismatch: PETSc range ["
            << ownershipStart << ',' << ownershipEnd << "), expected start "
            << expectedStart << " and local size " << 2*Pim.size()
            << exit(FatalError);
    }

    const PetscScalar* xArr = nullptr;
    VecGetArrayRead(x, &xArr);
    forAll(Pim, i)
    {
        Pim[i] = (scalar)PetscRealPart(xArr[2*i]);
        Pre[i] = (scalar)PetscRealPart(xArr[2*i + 1]);
    }
    VecRestoreArrayRead(x, &xArr);
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
