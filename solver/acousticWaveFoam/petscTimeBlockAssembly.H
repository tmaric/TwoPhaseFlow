// PETSc assembly helpers for the monolithic time-domain PML system.

static constexpr PetscInt timeBlockSize = 6;
static constexpr PetscInt pComponent = 0;
static constexpr PetscInt qComponent = 1;
static constexpr PetscInt vComponent = 2;
static constexpr PetscInt uxComponent = 3;
static constexpr PetscInt uyComponent = 4;
static constexpr PetscInt uzComponent = 5;

static bool sameTimeBlockPatchLayout(const labelList& a, const labelList& b)
{
    if (a.size() != b.size())
    {
        return false;
    }

    forAll(a, i)
    {
        if (a[i] != b[i])
        {
            return false;
        }
    }

    return true;
}

static void setDefaultTimePetscOption(const char* name, const char* value)
{
    PetscBool optionSet = PETSC_FALSE;
    PetscOptionsHasName(nullptr, nullptr, name, &optionSet);
    if (!optionSet)
    {
        PetscOptionsSetValue(nullptr, name, value);
    }
}

static const labelList& timeBlockLocalGlobalIds
(
    const fvMesh& mesh,
    const globalIndex& globalCells
)
{
    static const fvMesh* cachedMesh = nullptr;
    static label cachedNLocal = -1;
    static label cachedNGlobal = -1;
    static labelList localGlobalIds;

    if
    (
        cachedMesh != &mesh
     || cachedNLocal != mesh.nCells()
     || cachedNGlobal != globalCells.size()
    )
    {
        localGlobalIds.setSize(mesh.nCells());
        forAll(localGlobalIds, celli)
        {
            localGlobalIds[celli] = globalCells.toGlobal(celli);
        }

        cachedMesh = &mesh;
        cachedNLocal = mesh.nCells();
        cachedNGlobal = globalCells.size();
    }

    return localGlobalIds;
}

static const List<labelList>& timeBlockProcessorNeighbourGlobals
(
    const fvScalarMatrix& op,
    const globalIndex& globalCells
)
{
    const auto& bpsi = op.psi().boundaryField();
    const lduInterfaceFieldPtrsList interfaces = bpsi.scalarInterfaces();

    static const fvMesh* cachedMesh = nullptr;
    static label cachedNLocal = -1;
    static label cachedNGlobal = -1;
    static labelList cachedPatchFaces;
    static List<labelList> neighbourGlobals;

    if (!Pstream::parRun())
    {
        neighbourGlobals.setSize(bpsi.size());
        return neighbourGlobals;
    }

    labelList currentPatchFaces(bpsi.size(), -1);
    forAll(bpsi, patchi)
    {
        if (!interfaces.set(patchi))
        {
            continue;
        }

        const lduInterface& interface = interfaces[patchi].interface();
        if (isA<processorLduInterface>(interface))
        {
            currentPatchFaces[patchi] = interface.faceCells().size();
        }
    }

    if
    (
        cachedMesh == &op.psi().mesh()
     && cachedNLocal == op.psi().mesh().nCells()
     && cachedNGlobal == globalCells.size()
     && sameTimeBlockPatchLayout(cachedPatchFaces, currentPatchFaces)
    )
    {
        return neighbourGlobals;
    }

    neighbourGlobals.setSize(bpsi.size());
    const labelList& localGlobalIds =
        timeBlockLocalGlobalIds(op.psi().mesh(), globalCells);

    const label startOfRequests = UPstream::nRequests();
    forAll(bpsi, patchi)
    {
        if (!interfaces.set(patchi))
        {
            continue;
        }

        const lduInterface& interface = interfaces[patchi].interface();
        if (!isA<processorLduInterface>(interface))
        {
            continue;
        }

        interface.initInternalFieldTransfer
        (
            Pstream::commsTypes::nonBlocking,
            localGlobalIds,
            interface.faceCells()
        );
    }

    UPstream::waitRequests(startOfRequests);

    forAll(bpsi, patchi)
    {
        if (!interfaces.set(patchi))
        {
            continue;
        }

        const lduInterface& interface = interfaces[patchi].interface();
        if (!isA<processorLduInterface>(interface))
        {
            continue;
        }

        neighbourGlobals[patchi] = interface.internalFieldTransfer
        (
            Pstream::commsTypes::nonBlocking,
            localGlobalIds
        );
    }

    cachedMesh = &op.psi().mesh();
    cachedNLocal = op.psi().mesh().nCells();
    cachedNGlobal = globalCells.size();
    cachedPatchFaces.transfer(currentPatchFaces);
    return neighbourGlobals;
}

static void buildTimeBlockPreallocation
(
    const fvMesh& mesh,
    List<PetscInt>& dNnz,
    List<PetscInt>& oNnz
)
{
    labelList localNeighbours(mesh.nCells(), 0);
    labelList remoteNeighbours(mesh.nCells(), 0);

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    forAll(neighbour, facei)
    {
        ++localNeighbours[owner[facei]];
        ++localNeighbours[neighbour[facei]];
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
            ++remoteNeighbours[faceCells[facei]];
        }
    }

    dNnz.setSize(timeBlockSize*mesh.nCells(), 2);
    oNnz.setSize(timeBlockSize*mesh.nCells(), 0);

    forAll(localNeighbours, celli)
    {
        const label nLocal = localNeighbours[celli];
        const label nRemote = remoteNeighbours[celli];

        dNnz[timeBlockSize*celli + pComponent] = 2;
        dNnz[timeBlockSize*celli + qComponent] = 6 + 4*nLocal;
        dNnz[timeBlockSize*celli + vComponent] = 2;

        oNnz[timeBlockSize*celli + qComponent] = 4*nRemote;

        for (PetscInt cmpt = uxComponent; cmpt <= uzComponent; ++cmpt)
        {
            dNnz[timeBlockSize*celli + cmpt] = 3 + 2*nLocal;
            oNnz[timeBlockSize*celli + cmpt] = 2*nRemote;
        }
    }
}

static void initializeTimeBlockPetscSystem
(
    const fvMesh& mesh,
    const PetscInt nLocal,
    const PetscInt nGlobal,
    Mat& matrix,
    Vec& solution,
    Vec& rhs,
    KSP& ksp,
    const word& linearSolver
)
{
    MatCreate(PETSC_COMM_WORLD, &matrix);
    MatSetSizes
    (
        matrix,
        timeBlockSize*nLocal,
        timeBlockSize*nLocal,
        timeBlockSize*nGlobal,
        timeBlockSize*nGlobal
    );
    MatSetType(matrix, MATAIJ);

    List<PetscInt> dNnz;
    List<PetscInt> oNnz;
    buildTimeBlockPreallocation(mesh, dNnz, oNnz);
    if (Pstream::parRun())
    {
        MatMPIAIJSetPreallocation
        (
            matrix,
            0,
            dNnz.begin(),
            0,
            oNnz.begin()
        );
    }
    else
    {
        MatSeqAIJSetPreallocation(matrix, 0, dNnz.begin());
    }

    MatSetUp(matrix);
    MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

    VecCreate(PETSC_COMM_WORLD, &solution);
    VecSetSizes
    (
        solution,
        timeBlockSize*nLocal,
        timeBlockSize*nGlobal
    );
    VecSetFromOptions(solution);
    VecSet(solution, 0.0);
    VecDuplicate(solution, &rhs);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, matrix, matrix);

    if (linearSolver == "iterative")
    {
        setDefaultTimePetscOption("-ksp_type", "fgmres");
        setDefaultTimePetscOption("-ksp_gmres_restart", "50");
        setDefaultTimePetscOption("-ksp_rtol", "1e-8");
        setDefaultTimePetscOption("-ksp_max_it", "300");
        setDefaultTimePetscOption
        (
            "-ksp_norm_type",
            "unpreconditioned"
        );
        setDefaultTimePetscOption("-ksp_initial_guess_nonzero", "true");
        setDefaultTimePetscOption("-ksp_diagonal_scale", "true");
        setDefaultTimePetscOption("-ksp_diagonal_scale_fix", "true");
        setDefaultTimePetscOption("-pc_type", "bjacobi");
        setDefaultTimePetscOption("-sub_ksp_type", "preonly");
        setDefaultTimePetscOption("-sub_pc_type", "ilu");
        setDefaultTimePetscOption("-sub_pc_factor_levels", "1");
        setDefaultTimePetscOption
        (
            "-sub_pc_factor_shift_type",
            "nonzero"
        );
    }
    else
    {
        setDefaultTimePetscOption("-ksp_type", "preonly");
        setDefaultTimePetscOption("-pc_type", "lu");
        setDefaultTimePetscOption
        (
            "-pc_factor_mat_solver_type",
            "mumps"
        );
    }

    KSPSetFromOptions(ksp);
}

static scalarField timeBlockSourceWithBoundary(const fvScalarMatrix& op)
{
    scalarField source(op.source());
    const auto& bpsi = op.psi().boundaryField();
    const auto& boundaryCoeffs = op.boundaryCoeffs();
    const auto& addressing = op.lduAddr();
    const polyBoundaryMesh& patches = op.psi().mesh().boundaryMesh();

    forAll(bpsi, patchi)
    {
        const fvPatchField<scalar>& patchField = bpsi[patchi];
        if (!patchField.size())
        {
            continue;
        }

        if
        (
            patchField.coupled()
         && isA<processorPolyPatch>(patches[patchi])
        )
        {
            continue;
        }

        const Field<scalar>& patchCoeffs = boundaryCoeffs[patchi];
        const labelUList& patchAddressing = addressing.patchAddr(patchi);

        if (!patchField.coupled())
        {
            forAll(patchAddressing, facei)
            {
                source[patchAddressing[facei]] += patchCoeffs[facei];
            }
        }
        else
        {
            const tmp<Field<scalar>> neighbour =
                patchField.patchNeighbourField();
            forAll(patchAddressing, facei)
            {
                source[patchAddressing[facei]] +=
                    patchCoeffs[facei]*neighbour()[facei];
            }
        }
    }

    return source;
}

static void insertTimeBlockOperator
(
    Mat& matrix,
    const fvScalarMatrix& op,
    const globalIndex& globalCells,
    const PetscInt rowComponent,
    const PetscInt colComponent,
    const scalar scale,
    const scalarField* rowScale = nullptr
)
{
    const lduMatrix& ldu = op;
    scalarField diagonal(ldu.diag());
    const labelList& localGlobalIds =
        timeBlockLocalGlobalIds(op.psi().mesh(), globalCells);
    const List<labelList>& neighbourGlobals =
        timeBlockProcessorNeighbourGlobals(op, globalCells);

    const auto rowMultiplier = [&](const label celli)
    {
        return scale*(rowScale ? (*rowScale)[celli] : scalar(1));
    };

    const auto& bpsi = op.psi().boundaryField();
    const auto& internalCoeffs = op.internalCoeffs();
    const auto& boundaryCoeffs = op.boundaryCoeffs();
    const lduAddressing& addressing = op.lduAddr();

    forAll(bpsi, patchi)
    {
        if (bpsi[patchi].coupled())
        {
            continue;
        }

        const labelUList& patchAddressing = addressing.patchAddr(patchi);
        const Field<scalar>& patchInternalCoeffs = internalCoeffs[patchi];
        forAll(patchAddressing, facei)
        {
            diagonal[patchAddressing[facei]] += patchInternalCoeffs[facei];
        }
    }

    if (Pstream::parRun())
    {
        const polyBoundaryMesh& patches = op.psi().mesh().boundaryMesh();
        forAll(bpsi, patchi)
        {
            if
            (
                !bpsi[patchi].coupled()
             || !isA<processorPolyPatch>(patches[patchi])
            )
            {
                continue;
            }

            const labelUList& patchAddressing = addressing.patchAddr(patchi);
            const Field<scalar>& patchInternalCoeffs = internalCoeffs[patchi];
            const Field<scalar>& patchBoundaryCoeffs = boundaryCoeffs[patchi];
            const labelList& neighbourIds = neighbourGlobals[patchi];

            forAll(patchAddressing, facei)
            {
                const label celli = patchAddressing[facei];
                diagonal[celli] += patchInternalCoeffs[facei];

                const PetscInt row =
                    timeBlockSize*localGlobalIds[celli] + rowComponent;
                const PetscInt column =
                    timeBlockSize*neighbourIds[facei] + colComponent;
                const PetscScalar value =
                    -rowMultiplier(celli)*patchBoundaryCoeffs[facei];
                MatSetValue(matrix, row, column, value, ADD_VALUES);
            }
        }
    }

    const labelUList& owner = addressing.lowerAddr();
    const labelUList& neighbour = addressing.upperAddr();
    if (ldu.hasUpper())
    {
        const scalarField& upper = ldu.upper();
        const scalarField* lower = ldu.hasLower() ? &ldu.lower() : nullptr;

        forAll(upper, facei)
        {
            const label own = owner[facei];
            const label nei = neighbour[facei];

            const PetscInt ownRow =
                timeBlockSize*localGlobalIds[own] + rowComponent;
            const PetscInt neiRow =
                timeBlockSize*localGlobalIds[nei] + rowComponent;
            const PetscInt ownColumn =
                timeBlockSize*localGlobalIds[own] + colComponent;
            const PetscInt neiColumn =
                timeBlockSize*localGlobalIds[nei] + colComponent;

            const PetscScalar ownValue =
                rowMultiplier(own)*upper[facei];
            const PetscScalar neiValue =
                rowMultiplier(nei)*(lower ? (*lower)[facei] : upper[facei]);

            MatSetValue(matrix, ownRow, neiColumn, ownValue, ADD_VALUES);
            MatSetValue(matrix, neiRow, ownColumn, neiValue, ADD_VALUES);
        }
    }

    forAll(diagonal, celli)
    {
        const PetscInt row =
            timeBlockSize*localGlobalIds[celli] + rowComponent;
        const PetscInt column =
            timeBlockSize*localGlobalIds[celli] + colComponent;
        const PetscScalar value = rowMultiplier(celli)*diagonal[celli];
        MatSetValue(matrix, row, column, value, ADD_VALUES);
    }
}

static void addTimeBlockCellDiagonal
(
    Mat& matrix,
    const globalIndex& globalCells,
    const fvMesh& mesh,
    const PetscInt rowComponent,
    const PetscInt colComponent,
    const scalarField& coefficients
)
{
    const labelList& localGlobalIds =
        timeBlockLocalGlobalIds(mesh, globalCells);

    forAll(coefficients, celli)
    {
        const PetscInt row =
            timeBlockSize*localGlobalIds[celli] + rowComponent;
        const PetscInt column =
            timeBlockSize*localGlobalIds[celli] + colComponent;
        const PetscScalar value = coefficients[celli];
        MatSetValue(matrix, row, column, value, ADD_VALUES);
    }
}

static void addTimeBlockOperatorSource
(
    scalarField& rhs,
    const fvScalarMatrix& op,
    const scalar scale,
    const scalarField* rowScale = nullptr
)
{
    const scalarField source(timeBlockSourceWithBoundary(op));
    forAll(rhs, celli)
    {
        rhs[celli] +=
            scale*(rowScale ? (*rowScale)[celli] : scalar(1))*source[celli];
    }
}

static void setTimeBlockRhs
(
    Vec& rhs,
    const scalarField& pRhs,
    const scalarField& qRhs,
    const scalarField& vRhs,
    const scalarField& uxRhs,
    const scalarField& uyRhs,
    const scalarField& uzRhs
)
{
    PetscScalar* values = nullptr;
    VecGetArray(rhs, &values);

    forAll(pRhs, celli)
    {
        values[timeBlockSize*celli + pComponent] = pRhs[celli];
        values[timeBlockSize*celli + qComponent] = qRhs[celli];
        values[timeBlockSize*celli + vComponent] = vRhs[celli];
        values[timeBlockSize*celli + uxComponent] = uxRhs[celli];
        values[timeBlockSize*celli + uyComponent] = uyRhs[celli];
        values[timeBlockSize*celli + uzComponent] = uzRhs[celli];
    }

    VecRestoreArray(rhs, &values);
}

static void scatterTimeBlockSolution
(
    Vec& solution,
    const globalIndex& globalCells,
    volScalarField& p,
    volScalarField& q,
    volScalarField& v,
    volScalarField& ux,
    volScalarField& uy,
    volScalarField& uz,
    volVectorField& u
)
{
    PetscInt ownershipStart = 0;
    PetscInt ownershipEnd = 0;
    VecGetOwnershipRange(solution, &ownershipStart, &ownershipEnd);

    const PetscInt expectedStart = p.size()
      ? timeBlockSize*static_cast<PetscInt>(globalCells.toGlobal(0))
      : ownershipStart;

    if
    (
        ownershipStart != expectedStart
     || ownershipEnd - ownershipStart
        != timeBlockSize*static_cast<PetscInt>(p.size())
    )
    {
        FatalErrorInFunction
            << "PETSc/OpenFOAM ownership mismatch for time block: ["
            << ownershipStart << ',' << ownershipEnd << "), expected start "
            << expectedStart << " and local size " << timeBlockSize*p.size()
            << exit(FatalError);
    }

    const PetscScalar* values = nullptr;
    VecGetArrayRead(solution, &values);
    forAll(p, celli)
    {
        p[celli] = PetscRealPart(values[timeBlockSize*celli + pComponent]);
        q[celli] = PetscRealPart(values[timeBlockSize*celli + qComponent]);
        v[celli] = PetscRealPart(values[timeBlockSize*celli + vComponent]);
        ux[celli] = PetscRealPart(values[timeBlockSize*celli + uxComponent]);
        uy[celli] = PetscRealPart(values[timeBlockSize*celli + uyComponent]);
        uz[celli] = PetscRealPart(values[timeBlockSize*celli + uzComponent]);
    }
    VecRestoreArrayRead(solution, &values);

    p.correctBoundaryConditions();
    q.correctBoundaryConditions();
    v.correctBoundaryConditions();
    ux.correctBoundaryConditions();
    uy.correctBoundaryConditions();
    uz.correctBoundaryConditions();

    u.replace(vector::X, ux);
    u.replace(vector::Y, uy);
    u.replace(vector::Z, uz);
    u.correctBoundaryConditions();
}
