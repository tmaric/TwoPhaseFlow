// Acoustic block assembly helpers for interFALFlow.
//
// Design intent:
// - Keep the tensor-PML algebra identical to acousticHelmholtzFoam and the
//   serial block reference system.
// - Assemble matrix entries explicitly into PETSc, including processor links.
// - Treat B1/B2 in off-diagonal blocks consistently with the current
//   acousticHelmholtzFoam formulation: the off-diagonal PETSc blocks contain
//   the full B1 + B2 operator.

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
    static label cachedNLocal = -1;
    static label cachedNGlobal = -1;
    static labelList localGlobalIds;

    if
    (
        cachedNLocal != mesh.nCells()
     || cachedNGlobal != globalCells.size()
    )
    {
        localGlobalIds.setSize(mesh.nCells());
        forAll(localGlobalIds, celli)
        {
            localGlobalIds[celli] = globalCells.toGlobal(celli);
        }

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
        if (!isA<processorLduInterface>(lduInt))
        {
            continue;
        }

        currentPatchFaces[patchi] = lduInt.faceCells().size();
    }

    const bool needsRebuild =
        cachedNLocal != op.psi().mesh().nCells()
     || cachedNGlobal != globalCells.size()
     || !samePatchFaceLayout(cachedPatchFaces, currentPatchFaces);

    if (!needsRebuild)
    {
        return neighbGlobals;
    }

    neighbGlobals.setSize(bpsi.size());
    const labelList& localGlobalIds =
        cachedLocalGlobalIds(op.psi().mesh(), globalCells);

    // Processor-neighbour IDs depend only on mesh topology, so exchange once.
    // Post both sides before waiting to avoid deadlock above MPI's eager limit.
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
        if
        (
            neighbGlobals[patchi].size() != faceCells.size()
        )
        {
            FatalErrorInFunction
                << "Processor exchange size mismatch on patch index " << patchi
                << ": expected " << faceCells.size()
                << " got ids " << neighbGlobals[patchi].size()
                << exit(FatalError);
        }
    }

    cachedNLocal = op.psi().mesh().nCells();
    cachedNGlobal = globalCells.size();
    cachedPatchFaces.transfer(currentPatchFaces);

    return neighbGlobals;
}

static void exchangeProcessorPatchBoundaryCoeffs
(
    const fvScalarMatrix& op,
    List<scalarField>& neighbBcoeffs
)
{
    const auto& bpsi = op.psi().boundaryField();
    const auto& bcoeffs = op.boundaryCoeffs();
    const lduInterfaceFieldPtrsList interfaces =
        op.psi().boundaryField().scalarInterfaces();

    neighbBcoeffs.setSize(bpsi.size());

    if (!Pstream::parRun())
    {
        return;
    }

    // nonBlocking send posts both the receive and send into interface-owned
    // buffers. Wait before copying the received values into output fields.
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

        const processorLduInterface& procInt =
            refCast<const processorLduInterface>(lduInt);
        procInt.send(Pstream::commsTypes::nonBlocking, bcoeffs[patchi]);
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

        const processorLduInterface& procInt =
            refCast<const processorLduInterface>(lduInt);
        const label nFaces = lduInt.faceCells().size();

        scalarField recvCoeff(nFaces);
        procInt.receive(Pstream::commsTypes::nonBlocking, recvCoeff);
        neighbBcoeffs[patchi].transfer(recvCoeff);
    }
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
    const PetscInt rowOffset,
    const PetscInt colOffset,
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

    const List<labelList>& neighbGlobals =
        cachedProcessorPatchNeighbourGlobals(op, globalCells);
    List<scalarField> neighbBcoeffs;
    // Processor coupling is inserted explicitly so PETSc sees the same
    // operator that OpenFOAM would otherwise apply implicitly via interfaces.
    if (Pstream::parRun())
    {
        const auto& bpsi = op.psi().boundaryField();
        const auto& bcoeffs = op.boundaryCoeffs();
        const auto& icoeffs = op.internalCoeffs();
        const lduAddressing& baddr = op.lduAddr();
        exchangeProcessorPatchBoundaryCoeffs(op, neighbBcoeffs);

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
            const scalarField& nbrBcoeffs = neighbBcoeffs[patchi];

            if (nbrIds.size() != paddr.size() || nbrBcoeffs.size() != paddr.size())
            {
                FatalErrorInFunction
                    << "Processor patch exchange size mismatch on patch "
                    << op.psi().mesh().boundaryMesh()[patchi].name()
                    << ": local faces " << paddr.size()
                    << " neighbour ids " << nbrIds.size()
                    << " neighbour coeffs " << nbrBcoeffs.size()
                    << exit(FatalError);
            }

            forAll(paddr, facei)
            {
                const label ownCell = paddr[facei];
                const PetscInt r = localGlobalIds[ownCell] + rowOffset;
                const PetscInt cNbr = nbrIds[facei] + colOffset;

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

                MatSetValue(M, r, cNbr, -scale*pbc[facei], ADD_VALUES);
            }
        }
    }

    if (includeDiag)
    {
        forAll(diagWithProc, i)
        {
            PetscInt r = localGlobalIds[i] + rowOffset;
            PetscInt c = localGlobalIds[i] + colOffset;
            MatSetValue(M, r, c, scale*diagWithProc[i], ADD_VALUES);
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

            const PetscInt ri = localGlobalIds[i] + rowOffset;
            const PetscInt rj = localGlobalIds[j] + rowOffset;

            const PetscInt ci = localGlobalIds[i] + colOffset;
            const PetscInt cj = localGlobalIds[j] + colOffset;

            const scalar upperVal = upper[f];
            const scalar lowerVal = lowerPtr ? (*lowerPtr)[f] : upperVal;

            MatSetValue(M, ri, cj, scale*upperVal, ADD_VALUES);
            MatSetValue(M, rj, ci, scale*lowerVal, ADD_VALUES);
        }
    }
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
    AopPre.boundaryManipulate(AopPre.psi().boundaryFieldRef());
    AopPim.boundaryManipulate(AopPim.psi().boundaryFieldRef());
    B1Pre.boundaryManipulate(B1Pre.psi().boundaryFieldRef());
    B2Pre.boundaryManipulate(B2Pre.psi().boundaryFieldRef());
    B1Pim.boundaryManipulate(B1Pim.psi().boundaryFieldRef());
    B2Pim.boundaryManipulate(B2Pim.psi().boundaryFieldRef());

    insertScalarOpIntoBlock(M, AopPim, globalCells, 0, 0, +1.0);
    insertScalarOpIntoBlock(M, B1Pre, globalCells, 0, N, -1.0, true);
    insertScalarOpIntoBlock(M, B2Pre, globalCells, 0, N, -1.0, true);

    insertScalarOpIntoBlock(M, B1Pim, globalCells, N, 0, +1.0, true);
    insertScalarOpIntoBlock(M, B2Pim, globalCells, N, 0, +1.0, true);
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
