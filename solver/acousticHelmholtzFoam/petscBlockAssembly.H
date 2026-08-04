// Block assembly helpers for acousticHelmholtzFoam.
//
// Design intent:
// - Keep the algebra identical to the serial reference block system.
// - Assemble matrix entries explicitly into PETSc, including processor links.
// - Treat B1/B2 in off-diagonal blocks consistently with serial formulation:
//   the off-diagonal PETSc blocks should contain the full B1 + B2 operator.

static bool samePatchFaceLayout
(
    const labelList& a,
    const labelList& b
)
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

static const labelList& cachedLocalGlobalIds
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

static const List<labelList>& cachedProcessorPatchNeighbourGlobals
(
    const fvScalarMatrix& op,
    const globalIndex& globalCells
)
{
    const auto& bpsi = op.psi().boundaryField();
    const lduInterfaceFieldPtrsList interfaces =
        op.psi().boundaryField().scalarInterfaces();

    static const fvMesh* cachedMesh = nullptr;
    static label cachedNLocal = -1;
    static label cachedNGlobal = -1;
    static labelList cachedPatchFaces;
    static List<labelList> neighbGlobals;

    if (!Pstream::parRun())
    {
        neighbGlobals.setSize(bpsi.size());
        return neighbGlobals;
    }

    labelList currentPatchFaces(bpsi.size(), -1);
    forAll(bpsi, patchi)
    {
        if (!interfaces.set(patchi))
        {
            continue;
        }
        const lduInterface& lduInt = interfaces[patchi].interface();
        if (isA<processorLduInterface>(lduInt))
        {
            currentPatchFaces[patchi] = lduInt.faceCells().size();
        }
    }

    if
    (
        cachedMesh == &op.psi().mesh()
     && cachedNLocal == op.psi().mesh().nCells()
     && cachedNGlobal == globalCells.size()
     && samePatchFaceLayout(cachedPatchFaces, currentPatchFaces)
    )
    {
        return neighbGlobals;
    }

    neighbGlobals.setSize(bpsi.size());
    const labelList& localGlobalIds =
        cachedLocalGlobalIds(op.psi().mesh(), globalCells);

    // Processor-neighbour IDs depend only on mesh topology, so exchange once.
    // Use nonblocking transfers: posting scheduled/blocking sends on both sides
    // can deadlock when a processor-patch message exceeds MPI's eager limit.
    const label startOfRequests = UPstream::nRequests();
    forAll(bpsi, patchi)
    {
        if (!interfaces.set(patchi))
        {
            continue;
        }

        const lduInterface& lduInt = interfaces[patchi].interface();
        if (!isA<processorLduInterface>(lduInt))
        {
            continue;
        }

        const labelUList& faceCells = lduInt.faceCells();

        lduInt.initInternalFieldTransfer
        (
            Pstream::commsTypes::nonBlocking,
            localGlobalIds,
            faceCells
        );
    }

    UPstream::waitRequests(startOfRequests);

    forAll(bpsi, patchi)
    {
        if (!interfaces.set(patchi))
        {
            continue;
        }

        const lduInterface& lduInt = interfaces[patchi].interface();
        if (!isA<processorLduInterface>(lduInt))
        {
            continue;
        }

        const labelUList& faceCells = lduInt.faceCells();
        neighbGlobals[patchi] =
            lduInt.internalFieldTransfer
            (
                Pstream::commsTypes::nonBlocking,
                localGlobalIds
            );

        if (neighbGlobals[patchi].size() != faceCells.size())
        {
            FatalErrorInFunction
                << "Processor exchange size mismatch on patch index " << patchi
                << ": expected " << faceCells.size()
                << " got ids " << neighbGlobals[patchi].size()
                << exit(FatalError);
        }
    }

    cachedMesh = &op.psi().mesh();
    cachedNLocal = op.psi().mesh().nCells();
    cachedNGlobal = globalCells.size();
    cachedPatchFaces.transfer(currentPatchFaces);
    return neighbGlobals;
}

static void addBoundarySourceSimple
(
    const fvScalarMatrix& M,
    scalarField& source,
    const bool couples = true,
    const bool includeProcessorCoupled = true
)
{
    const auto& bpsi = M.psi().boundaryField();
    const auto& bcoeffs = M.boundaryCoeffs();
    const auto& addr = M.lduAddr();
    const polyBoundaryMesh& patches = M.psi().mesh().boundaryMesh();

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
            if
            (
                !includeProcessorCoupled
             && isA<processorPolyPatch>(patches[patchi])
            )
            {
                continue;
            }
            const tmp<Field<scalar>> tpnf = ptf.patchNeighbourField();
            const Field<scalar>& pnf = tpnf();
            forAll(paddr, facei)
            {
                source[paddr[facei]] += pbc[facei]*pnf[facei];
            }
        }
    }
}

static scalarField sourceWithBoundary
(
    const fvScalarMatrix& M,
    const bool includeCoupled,
    const bool includeProcessorCoupled = true
)
{
    scalarField source(M.source());
    addBoundarySourceSimple(M, source, includeCoupled, includeProcessorCoupled);
    return source;
}

static void insertScalarOpIntoBlock
(
    Mat& M,
    const fvScalarMatrix& op,
    const globalIndex& globalCells,
    const PetscInt rowComponent,
    const PetscInt colComponent,
    const PetscScalar scale,
    const bool includeDiag = true
)
{
    const lduMatrix& L = op;
    const scalarField& diag = L.diag();
    scalarField diagWithProc(diag);
    const labelList& localGlobalIds =
        cachedLocalGlobalIds(op.psi().mesh(), globalCells);

    const lduAddressing& addr = L.lduAddr();
    const labelUList& owner = addr.lowerAddr();
    const labelUList& neigh = addr.upperAddr();

    List<DynamicList<PetscInt>> rowCols(op.psi().mesh().nCells());
    List<DynamicList<PetscScalar>> rowVals(op.psi().mesh().nCells());
    const List<labelList>& neighbGlobals =
        cachedProcessorPatchNeighbourGlobals(op, globalCells);

    // Physical boundary conditions contribute both a source coefficient and
    // an internal coefficient. OpenFOAM normally adds the latter to the
    // diagonal in fvMatrix::solve(); PETSc assembly must do this explicitly.
    if (includeDiag)
    {
        const auto& bpsi = op.psi().boundaryField();
        const auto& icoeffs = op.internalCoeffs();
        const lduAddressing& baddr = op.lduAddr();
        forAll(bpsi, patchi)
        {
            if (bpsi[patchi].coupled())
            {
                continue;
            }

            const labelUList& paddr = baddr.patchAddr(patchi);
            const Field<scalar>& pic = icoeffs[patchi];
            forAll(paddr, facei)
            {
                diagWithProc[paddr[facei]] += pic[facei];
            }
        }
    }

    // Processor coupling is explicit here so PETSc sees the same operator that
    // OpenFOAM would otherwise apply implicitly through ldu interfaces.
    if (Pstream::parRun())
    {
        const auto& bpsi = op.psi().boundaryField();
        const auto& bcoeffs = op.boundaryCoeffs();
        const auto& icoeffs = op.internalCoeffs();
        const lduAddressing& baddr = op.lduAddr();
        forAll(bpsi, patchi)
        {
            if
            (
                !bpsi[patchi].coupled()
             || !isA<processorPolyPatch>(op.psi().mesh().boundaryMesh()[patchi])
            )
            {
                continue;
            }

            const labelUList& paddr = baddr.patchAddr(patchi);
            const Field<scalar>& pbc = bcoeffs[patchi];
            const Field<scalar>& pic = icoeffs[patchi];
            const labelList& nbrIds = neighbGlobals[patchi];

            if (nbrIds.size() != paddr.size())
            {
                FatalErrorInFunction
                    << "Processor patch exchange size mismatch on patch "
                    << op.psi().mesh().boundaryMesh()[patchi].name()
                    << ": local faces " << paddr.size()
                    << " neighbour ids " << nbrIds.size()
                    << exit(FatalError);
            }

            forAll(paddr, facei)
            {
                const label ownCell = paddr[facei];
                const PetscInt cNbr = 2*nbrIds[facei] + colComponent;

                if (nbrIds[facei] < 0 || nbrIds[facei] >= globalCells.size())
                {
                    FatalErrorInFunction
                        << "Invalid neighbour global cell id on patch "
                        << op.psi().mesh().boundaryMesh()[patchi].name()
                        << " face " << facei
                        << ": " << nbrIds[facei]
                        << " not in [0," << (globalCells.size() - 1) << "]"
                        << exit(FatalError);
                }

                if (includeDiag)
                {
                    diagWithProc[ownCell] += pic[facei];
                }
                rowCols[ownCell].append(cNbr);
                rowVals[ownCell].append(-scale*pbc[facei]);
            }
        }
    }

    if (L.hasUpper())
    {
        const scalarField& upper = L.upper();
        const scalarField* lowerPtr = L.hasLower() ? &L.lower() : nullptr;

        forAll(upper, f)
        {
            const label i = owner[f];
            const label j = neigh[f];

            const PetscInt ci = 2*localGlobalIds[i] + colComponent;
            const PetscInt cj = 2*localGlobalIds[j] + colComponent;

            const scalar upperVal = upper[f];
            const scalar lowerVal = lowerPtr ? (*lowerPtr)[f] : upperVal;

            rowCols[i].append(cj);
            rowVals[i].append(scale*upperVal);
            rowCols[j].append(ci);
            rowVals[j].append(scale*lowerVal);
        }
    }

    forAll(rowCols, i)
    {
        if (includeDiag)
        {
            rowCols[i].append(2*localGlobalIds[i] + colComponent);
            rowVals[i].append(scale*diagWithProc[i]);
        }

        const PetscInt row = 2*localGlobalIds[i] + rowComponent;
        MatSetValues
        (
            M,
            1,
            &row,
            rowCols[i].size(),
            rowCols[i].begin(),
            rowVals[i].begin(),
            ADD_VALUES
        );
    }
}

static void assembleBlockSystem
(
    Mat& M,
    const globalIndex& globalCells,
    fvScalarMatrix& AopPim,
    fvScalarMatrix& AopPre,
    fvScalarMatrix& B1Pre,
    fvScalarMatrix& B2Pre,
    fvScalarMatrix& B1Pim,
    fvScalarMatrix& B2Pim
)
{
    AopPre.boundaryManipulate(AopPre.psi().boundaryFieldRef());
    AopPim.boundaryManipulate(AopPim.psi().boundaryFieldRef());
    B1Pre.boundaryManipulate(B1Pre.psi().boundaryFieldRef());
    B2Pre.boundaryManipulate(B2Pre.psi().boundaryFieldRef());
    B1Pim.boundaryManipulate(B1Pim.psi().boundaryFieldRef());
    B2Pim.boundaryManipulate(B2Pim.psi().boundaryFieldRef());

    insertScalarOpIntoBlock(M, AopPim, globalCells, 0, 0, +1.0);
    insertScalarOpIntoBlock(M, B1Pre, globalCells, 0, 1, -1.0, true);
    insertScalarOpIntoBlock(M, B2Pre, globalCells, 0, 1, -1.0, true);

    insertScalarOpIntoBlock(M, B1Pim, globalCells, 1, 0, +1.0, true);
    insertScalarOpIntoBlock(M, B2Pim, globalCells, 1, 0, +1.0, true);
    insertScalarOpIntoBlock(M, AopPre, globalCells, 1, 1, +1.0);
}

static void setBlockRhs
(
    Vec& b,
    const scalarField& bPim,
    const scalarField& bPre
)
{
    PetscScalar* bArr = nullptr;
    VecGetArray(b, &bArr);
    forAll(bPim, i)
    {
        bArr[2*i] = (PetscScalar)bPim[i];
        bArr[2*i + 1] = (PetscScalar)bPre[i];
    }
    VecRestoreArray(b, &bArr);
}
