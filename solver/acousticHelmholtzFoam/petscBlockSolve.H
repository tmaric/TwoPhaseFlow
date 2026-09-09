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
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsHasName(nullptr, nullptr, name, &optionSet));
    if (!optionSet)
    {
        PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsSetValue(nullptr, name, value));
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
    const globalIndex& globalCells,
    const std::vector<Foam::acoustic::FaceOperator>& interfaceOps,
    const PetscInt nLocal,
    const PetscInt N,
    Mat& M,
    Vec& x,
    Vec& b,
    KSP& ksp
)
{
    PetscCallAbort(PETSC_COMM_WORLD, MatCreate(PETSC_COMM_WORLD, &M));
    PetscCallAbort(PETSC_COMM_WORLD, MatSetSizes(M, 2*nLocal, 2*nLocal, 2*N, 2*N));
    PetscCallAbort(PETSC_COMM_WORLD, MatSetType(M, MATAIJ));

    List<PetscInt> dNnz;
    List<PetscInt> oNnz;
    buildBlockPreallocation(mesh, dNnz, oNnz);
    // Additional columns occur only in the same pressure-component block.
    // Count a superset of the union; retaining unused slots is harmless.
    List<labelHashSet> extra(mesh.nCells());
    for (const auto& op:interfaceOps)
    {
        for (const label id:op.cells)
        {
            extra[op.owner].insert(id);
            if(op.neighbour>=0) extra[op.neighbour].insert(id);
        }
    }
    forAll(extra,i) for(const label id:extra[i])
        for(label c=0;c<2;++c)
        {
            if(globalCells.isLocal(id)) ++dNnz[2*i+c];
            else ++oNnz[2*i+c];
        }
    if (Pstream::parRun())
    {
        PetscCallAbort(PETSC_COMM_WORLD, MatMPIAIJSetPreallocation(M, 0, dNnz.begin(), 0, oNnz.begin()));
    }
    else
    {
        PetscCallAbort(PETSC_COMM_WORLD, MatSeqAIJSetPreallocation(M, 0, dNnz.begin()));
    }
    PetscCallAbort(PETSC_COMM_WORLD, MatSetUp(M));
    PetscCallAbort(PETSC_COMM_WORLD, MatSetOption(M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));

    PetscCallAbort(PETSC_COMM_WORLD, VecCreate(PETSC_COMM_WORLD, &x));
    PetscCallAbort(PETSC_COMM_WORLD, VecSetSizes(x, 2*nLocal, 2*N));
    PetscCallAbort(PETSC_COMM_WORLD, VecSetFromOptions(x));
    PetscCallAbort(PETSC_COMM_WORLD, VecDuplicate(x, &b));

    PetscCallAbort(PETSC_COMM_WORLD, KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetOperators(ksp, M, M));

    setDefaultPetscOption("-ksp_type", "preonly");
    setDefaultPetscOption("-pc_type", "lu");
    setDefaultPetscOption("-pc_factor_mat_solver_type", "mumps");
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetFromOptions(ksp));
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
    PetscCallAbort(PETSC_COMM_WORLD, VecGetOwnershipRange(x, &ownershipStart, &ownershipEnd));

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
    PetscCallAbort(PETSC_COMM_WORLD, VecGetArrayRead(x, &xArr));
    forAll(Pim, i)
    {
        Pim[i] = (scalar)PetscRealPart(xArr[2*i]);
        Pre[i] = (scalar)PetscRealPart(xArr[2*i + 1]);
    }
    PetscCallAbort(PETSC_COMM_WORLD, VecRestoreArrayRead(x, &xArr));
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

static void verifyPetscSolution(Mat M,Vec x,Vec b,KSP ksp)
{
    KSPConvergedReason reason;
    PetscCallAbort(PETSC_COMM_WORLD, KSPGetConvergedReason(ksp,&reason));
    if(reason<=0)
        FatalErrorInFunction << "PETSc solve did not converge, reason=" << reason << exit(FatalError);
    Vec residual;
    PetscCallAbort(PETSC_COMM_WORLD, VecDuplicate(b,&residual));
    PetscCallAbort(PETSC_COMM_WORLD, MatMult(M,x,residual));
    PetscCallAbort(PETSC_COMM_WORLD, VecAXPY(residual,-1.0,b));
    PetscReal r2,b2,ri,bi,xi,mi;
    PetscCallAbort(PETSC_COMM_WORLD, VecNorm(residual,NORM_2,&r2));
    PetscCallAbort(PETSC_COMM_WORLD, VecNorm(b,NORM_2,&b2));
    PetscCallAbort(PETSC_COMM_WORLD, VecNorm(residual,NORM_INFINITY,&ri));
    PetscCallAbort(PETSC_COMM_WORLD, VecNorm(b,NORM_INFINITY,&bi));
    PetscCallAbort(PETSC_COMM_WORLD, VecNorm(x,NORM_INFINITY,&xi));
    PetscCallAbort(PETSC_COMM_WORLD, MatNorm(M,NORM_INFINITY,&mi));
    PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&residual));
    if(!std::isfinite(r2) || !std::isfinite(b2) || !std::isfinite(xi)
        || !std::isfinite(mi))
        FatalErrorInFunction << "Nonfinite PETSc solution or residual" << exit(FatalError);
    Info<< "Acoustic solve: reason=" << reason
        << " relativeResidual=" << r2/Foam::max(b2,PetscReal(VSMALL))
        << " backwardError=" << ri/Foam::max(mi*xi+bi,PetscReal(VSMALL)) << nl;
}
