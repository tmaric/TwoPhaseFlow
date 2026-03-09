// Block assembly helpers for acousticHelmholtzFoam.
//
// Design intent:
// - Keep the algebra identical to the serial reference block system.
// - Assemble matrix entries explicitly into PETSc, including processor links.
// - Treat B1/B2 in off-diagonal blocks consistently with serial formulation:
//   B1 contributes face-coupling only, B2 carries diagonal mass term.

static List<labelList> gatherProcessorPatchNeighbourGlobals
(
    const fvScalarMatrix& op,
    const globalIndex& globalCells,
    List<scalarField>& neighbBcoeffs
)
{
    const auto& bpsi = op.psi().boundaryField();
    const auto& bcoeffs = op.boundaryCoeffs();
    const lduInterfaceFieldPtrsList interfaces =
        op.psi().boundaryField().scalarInterfaces();

    List<labelList> neighbGlobals(bpsi.size());
    neighbBcoeffs.setSize(bpsi.size());

    if (!Pstream::parRun())
    {
        return neighbGlobals;
    }

    labelList localGlobalIds(op.psi().mesh().nCells());
    forAll(localGlobalIds, celli)
    {
        localGlobalIds[celli] = globalCells.toGlobal(celli);
    }

    // Stage 1: send local global-cell ids and patch coefficients to neighbours.
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
        const labelUList& faceCells = lduInt.faceCells();

        lduInt.initInternalFieldTransfer
        (
            Pstream::commsTypes::scheduled,
            localGlobalIds,
            faceCells
        );
        procInt.send(Pstream::commsTypes::scheduled, bcoeffs[patchi]);
    }

    // Stage 2: receive neighbour ids/coefficients in interface-face order.
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
        const labelUList& faceCells = lduInt.faceCells();
        neighbGlobals[patchi] =
            lduInt.internalFieldTransfer
            (
                Pstream::commsTypes::scheduled,
                localGlobalIds
            );

        scalarField recvCoeff(faceCells.size());
        procInt.receive(Pstream::commsTypes::scheduled, recvCoeff);
        if
        (
            neighbGlobals[patchi].size() != faceCells.size()
         || recvCoeff.size() != faceCells.size()
        )
        {
            FatalErrorInFunction
                << "Processor exchange size mismatch on patch index " << patchi
                << ": expected " << faceCells.size()
                << " got ids " << neighbGlobals[patchi].size()
                << " coeffs " << recvCoeff.size()
                << exit(FatalError);
        }

        neighbBcoeffs[patchi].transfer(recvCoeff);
    }

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
    const PetscInt rowOffset,
    const PetscInt colOffset,
    const PetscScalar scale,
    const bool includeDiag = true
)
{
    const lduMatrix& L = op;
    const scalarField& diag = L.diag();
    scalarField diagWithProc(diag);

    const lduAddressing& addr = L.lduAddr();
    const labelUList& owner = addr.lowerAddr();
    const labelUList& neigh = addr.upperAddr();

    List<labelList> neighbGlobals;
    List<scalarField> neighbBcoeffs;
    scalar procDiagAbs = 0.0;
    scalar procOffDiagAbs = 0.0;
    label procFaceCount = 0;
    // Processor coupling is explicit here so PETSc sees the same operator that
    // OpenFOAM would otherwise apply implicitly through ldu interfaces.
    if (Pstream::parRun())
    {
        const auto& bpsi = op.psi().boundaryField();
        const auto& bcoeffs = op.boundaryCoeffs();
        const auto& icoeffs = op.internalCoeffs();
        const lduAddressing& baddr = op.lduAddr();
        neighbGlobals =
            gatherProcessorPatchNeighbourGlobals(op, globalCells, neighbBcoeffs);
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
                const PetscInt r = globalCells.toGlobal(ownCell) + rowOffset;
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

                // For B1 off-block assembly, includeDiag is false: keep only the
                // face-coupling term and avoid injecting B1 diagonal entries.
                if (includeDiag)
                {
                    diagWithProc[ownCell] += pic[facei];
                }
                MatSetValue(M, r, cNbr, -scale*pbc[facei], INSERT_VALUES);
                if (includeDiag)
                {
                    procDiagAbs += mag(pic[facei]);
                }
                procOffDiagAbs += mag(pbc[facei]);
                ++procFaceCount;
            }
        }
    }

    if (Pstream::parRun())
    {
        procDiagAbs = globalSum(procDiagAbs);
        procOffDiagAbs = globalSum(procOffDiagAbs);
        reduce(procFaceCount, sumOp<label>());
        static label printedProcContribCount = 0;
        if (printedProcContribCount < 4)
        {
            if (Pstream::master())
            {
                Info<< "Proc-coupling contribution (first block): rowOffset="
                    << rowOffset << " colOffset=" << colOffset
                    << " nFaces=" << procFaceCount
                    << " |diagAdd|=" << procDiagAbs
                    << " |offDiag|=" << procOffDiagAbs
                    << nl;
            }
            ++printedProcContribCount;
        }
    }

    if (includeDiag)
    {
        forAll(diagWithProc, i)
        {
            PetscInt r = globalCells.toGlobal(i) + rowOffset;
            PetscInt c = globalCells.toGlobal(i) + colOffset;
            MatSetValue(M, r, c, scale*diagWithProc[i], INSERT_VALUES);
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

            const PetscInt ri = globalCells.toGlobal(i) + rowOffset;
            const PetscInt rj = globalCells.toGlobal(j) + rowOffset;

            const PetscInt ci = globalCells.toGlobal(i) + colOffset;
            const PetscInt cj = globalCells.toGlobal(j) + colOffset;

            const scalar upperVal = upper[f];
            const scalar lowerVal = lowerPtr ? (*lowerPtr)[f] : upperVal;

            MatSetValue(M, ri, cj, scale*upperVal, INSERT_VALUES);
            MatSetValue(M, rj, ci, scale*lowerVal, INSERT_VALUES);
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
    insertScalarOpIntoBlock(M, B1Pre, globalCells, 0, N, -1.0, false);
    insertScalarOpIntoBlock(M, B2Pre, globalCells, 0, N, -1.0, true);

    insertScalarOpIntoBlock(M, B1Pim, globalCells, N, 0, +1.0, false);
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
