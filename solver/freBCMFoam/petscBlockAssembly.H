// Distributed block-assembly helpers for freBCMFoam.
//
// This file maps OpenFOAM scalar operators into PETSc block rows/columns and
// explicitly handles processor-patch coupling when assembling diagonal blocks.
// Goal: preserve the same block equation structure as freBlockCoupledFoam
// while supporting decomposed MPI meshes.

static List<labelList> gatherProcessorPatchNeighbourGlobals
(
    const fvScalarMatrix& op,
    const globalIndex& globalCells
)
{
    const auto& bpsi = op.psi().boundaryField();
    const lduAddressing& addr = op.lduAddr();
    const polyBoundaryMesh& patches = op.psi().mesh().boundaryMesh();

    List<labelList> neighbGlobals(bpsi.size());

    if (!Pstream::parRun())
    {
        return neighbGlobals;
    }

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    forAll(bpsi, patchi)
    {
        if (!bpsi[patchi].coupled() || !isA<processorPolyPatch>(patches[patchi]))
        {
            continue;
        }

        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(patches[patchi]);
        const labelUList& paddr = addr.patchAddr(patchi);

        labelList sendIds(paddr.size());
        forAll(paddr, facei)
        {
            sendIds[facei] = globalCells.toGlobal(paddr[facei]);
        }

        UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
        toNbr << sendIds;
    }

    pBufs.finishedSends();

    forAll(bpsi, patchi)
    {
        if (!bpsi[patchi].coupled() || !isA<processorPolyPatch>(patches[patchi]))
        {
            continue;
        }

        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(patches[patchi]);
        UIPstream fromNbr(procPatch.neighbProcNo(), pBufs);
        fromNbr >> neighbGlobals[patchi];
    }

    return neighbGlobals;
}

static void addBoundarySourceSimple
(
    const fvScalarMatrix& M,
    scalarField& source,
    const bool couples = true
)
{
    const auto& bpsi = M.psi().boundaryField();
    const auto& bcoeffs = M.boundaryCoeffs();
    const auto& addr = M.lduAddr();

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
    const bool includeCoupled
)
{
    scalarField source(M.source());
    addBoundarySourceSimple(M, source, includeCoupled);
    return source;
}

static void insertScalarOpIntoBlock
(
    Mat& M,
    const fvScalarMatrix& op,
    const globalIndex& globalCells,
    const PetscInt rowOffset,
    const PetscInt colOffset,
    const PetscScalar scale
)
{
    const lduMatrix& L = op;

    const scalarField& diag = L.diag();
    scalarField diagWithProc(diag);

    const lduAddressing& addr = L.lduAddr();
    const labelUList& owner = addr.lowerAddr();
    const labelUList& neigh = addr.upperAddr();

    List<labelList> neighbGlobals;
    if (Pstream::parRun() && rowOffset == colOffset)
    {
        const auto& bpsi = op.psi().boundaryField();
        const auto& bcoeffs = op.boundaryCoeffs();
        const auto& icoeffs = op.internalCoeffs();
        const lduAddressing& baddr = op.lduAddr();
        const polyBoundaryMesh& patches = op.psi().mesh().boundaryMesh();
        neighbGlobals = gatherProcessorPatchNeighbourGlobals(op, globalCells);

        forAll(bpsi, patchi)
        {
            if (!bpsi[patchi].coupled() || !isA<processorPolyPatch>(patches[patchi]))
            {
                continue;
            }

            const labelUList& paddr = baddr.patchAddr(patchi);
            const Field<scalar>& pbc = bcoeffs[patchi];
            const labelList& nbrIds = neighbGlobals[patchi];

            if (nbrIds.size() != paddr.size())
            {
                FatalErrorInFunction
                    << "Processor patch exchange size mismatch on patch "
                    << patches[patchi].name()
                    << ": local faces " << paddr.size()
                    << " neighbour ids " << nbrIds.size()
                    << exit(FatalError);
            }

            const Field<scalar>& pic = icoeffs[patchi];

            forAll(paddr, facei)
            {
                const PetscInt cNbr = nbrIds[facei] + colOffset;
                const label ownCell = paddr[facei];
                const PetscInt r = globalCells.toGlobal(ownCell) + rowOffset;

                diagWithProc[ownCell] += pic[facei];
                MatSetValue(M, r, cNbr, -scale*pbc[facei], INSERT_VALUES);
            }
        }
    }

    forAll(diagWithProc, i)
    {
        PetscInt r = globalCells.toGlobal(i) + rowOffset;
        PetscInt c = globalCells.toGlobal(i) + colOffset;
        MatSetValue(M, r, c, scale*diagWithProc[i], INSERT_VALUES);
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

static void insertScalarOpTransposeIntoBlock
(
    Mat& M,
    const fvScalarMatrix& op,
    const globalIndex& globalCells,
    const PetscInt rowOffset,
    const PetscInt colOffset,
    const PetscScalar scale
)
{
    const lduMatrix& L = op;
    const scalarField& diag = L.diag();
    const lduAddressing& addr = L.lduAddr();
    const labelUList& owner = addr.lowerAddr();
    const labelUList& neigh = addr.upperAddr();

    forAll(diag, i)
    {
        const PetscInt r = globalCells.toGlobal(i) + rowOffset;
        const PetscInt c = globalCells.toGlobal(i) + colOffset;
        MatSetValue(M, r, c, scale*diag[i], INSERT_VALUES);
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

            MatSetValue(M, ri, cj, scale*lowerVal, INSERT_VALUES);
            MatSetValue(M, rj, ci, scale*upperVal, INSERT_VALUES);
        }
    }
}

static void assembleBlockSystem
(
    Mat& M,
    const globalIndex& globalCells,
    const PetscInt N,
    fvScalarMatrix& mainOpPim,
    fvScalarMatrix& mainOpPre,
    fvScalarMatrix& couplingOpPre
)
{
    mainOpPre.boundaryManipulate(mainOpPre.psi().boundaryFieldRef());
    mainOpPim.boundaryManipulate(mainOpPim.psi().boundaryFieldRef());
    couplingOpPre.boundaryManipulate(couplingOpPre.psi().boundaryFieldRef());

    insertScalarOpIntoBlock(M, mainOpPim, globalCells, 0, 0, +1.0);
    insertScalarOpIntoBlock(M, couplingOpPre, globalCells, 0, N, -1.0);

    insertScalarOpTransposeIntoBlock(M, couplingOpPre, globalCells, N, 0, +1.0);
    insertScalarOpIntoBlock(M, mainOpPre, globalCells, N, N, +1.0);
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
