static inline scalar twoPi() { return constant::mathematical::twoPi; }

static void addBoundarySourceSimple
(
    const fvScalarMatrix& M,
    scalarField& source,
    const bool includeProcessorCoupled
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
        else
        {
            const bool isProcessor = isA<processorPolyPatch>(patches[patchi]);
            if (includeProcessorCoupled || !isProcessor)
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
}

static scalarField sourceWithBoundary
(
    const fvScalarMatrix& M,
    const bool includeProcessorCoupled
)
{
    scalarField source(M.source());
    addBoundarySourceSimple(M, source, includeProcessorCoupled);
    return source;
}

static List<labelList> gatherProcessorPatchNeighbourGlobals
(
    const fvScalarMatrix& op,
    const globalIndex& globalCells
)
{
    const fvMesh& mesh = op.psi().mesh();
    const auto& bpsi = op.psi().boundaryField();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    List<labelList> neighbGlobals(bpsi.size());

    if (!Pstream::parRun())
    {
        return neighbGlobals;
    }

    const label nInternalFaces = mesh.nInternalFaces();
    const label nBoundaryFaces = mesh.nBoundaryFaces();
    const labelList& faceOwner = mesh.faceOwner();

    labelList neighbByBFace(nBoundaryFaces, -1);
    for (label bFacei = 0; bFacei < nBoundaryFaces; ++bFacei)
    {
        const label facei = nInternalFaces + bFacei;
        neighbByBFace[bFacei] = globalCells.toGlobal(faceOwner[facei]);
    }

    syncTools::swapBoundaryFaceList(mesh, neighbByBFace);

    forAll(patches, patchi)
    {
        if (!bpsi[patchi].coupled() || !isA<processorPolyPatch>(patches[patchi]))
        {
            continue;
        }

        const polyPatch& pp = patches[patchi];
        const label start = pp.start() - nInternalFaces;
        neighbGlobals[patchi].setSize(pp.size());
        forAll(neighbGlobals[patchi], facei)
        {
            neighbGlobals[patchi][facei] = neighbByBFace[start + facei];
        }
    }

    return neighbGlobals;
}

static void appendScalarOpTriplets
(
    const fvScalarMatrix& op,
    const globalIndex& globalCells,
    const PetscInt rowOffset,
    const PetscInt colOffset,
    const scalar scale,
    labelList& rows,
    labelList& cols,
    scalarField& vals
)
{
    const lduMatrix& L = op;
    const scalarField& diag = L.diag();
    scalarField diagWithProc(diag);

    const lduAddressing& addr = L.lduAddr();
    const labelUList& owner = addr.lowerAddr();
    const labelUList& neigh = addr.upperAddr();

    if (Pstream::parRun())
    {
        const auto& bpsi = op.psi().boundaryField();
        const auto& bcoeffs = op.boundaryCoeffs();
        const auto& icoeffs = op.internalCoeffs();
        const lduAddressing& baddr = op.lduAddr();
        const polyBoundaryMesh& patches = op.psi().mesh().boundaryMesh();
        const List<labelList> neighbGlobals
        (
            gatherProcessorPatchNeighbourGlobals(op, globalCells)
        );

        forAll(bpsi, patchi)
        {
            if (!bpsi[patchi].coupled() || !isA<processorPolyPatch>(patches[patchi]))
            {
                continue;
            }

            const labelUList& paddr = baddr.patchAddr(patchi);
            const Field<scalar>& pbc = bcoeffs[patchi];
            const Field<scalar>& pic = icoeffs[patchi];
            const labelList& nbrIds = neighbGlobals[patchi];

            forAll(paddr, facei)
            {
                const label own = paddr[facei];
                diagWithProc[own] += pic[facei];
                rows.append(globalCells.toGlobal(own) + rowOffset);
                cols.append(nbrIds[facei] + colOffset);
                vals.append(-scale*pbc[facei]);
            }
        }
    }

    forAll(diagWithProc, i)
    {
        rows.append(globalCells.toGlobal(i) + rowOffset);
        cols.append(globalCells.toGlobal(i) + colOffset);
        vals.append(scale*diagWithProc[i]);
    }

    if (L.hasUpper())
    {
        const scalarField& upper = L.upper();
        const scalarField* lowerPtr = L.hasLower() ? &L.lower() : nullptr;

        forAll(upper, f)
        {
            const label i = owner[f];
            const label j = neigh[f];
            const scalar upperVal = upper[f];
            const scalar lowerVal = lowerPtr ? (*lowerPtr)[f] : upperVal;

            rows.append(globalCells.toGlobal(i) + rowOffset);
            cols.append(globalCells.toGlobal(j) + colOffset);
            vals.append(scale*upperVal);

            rows.append(globalCells.toGlobal(j) + rowOffset);
            cols.append(globalCells.toGlobal(i) + colOffset);
            vals.append(scale*lowerVal);
        }
    }
}

static void appendBlockTriplets
(
    const globalIndex& globalCells,
    const PetscInt N,
    fvScalarMatrix& AopPim,
    fvScalarMatrix& AopPre,
    fvScalarMatrix& B1Pre,
    fvScalarMatrix& B2Pre,
    fvScalarMatrix& B1Pim,
    fvScalarMatrix& B2Pim,
    labelList& rows,
    labelList& cols,
    scalarField& vals
)
{
    AopPre.boundaryManipulate(AopPre.psi().boundaryFieldRef());
    AopPim.boundaryManipulate(AopPim.psi().boundaryFieldRef());
    B1Pre.boundaryManipulate(B1Pre.psi().boundaryFieldRef());
    B2Pre.boundaryManipulate(B2Pre.psi().boundaryFieldRef());
    B1Pim.boundaryManipulate(B1Pim.psi().boundaryFieldRef());
    B2Pim.boundaryManipulate(B2Pim.psi().boundaryFieldRef());

    appendScalarOpTriplets(AopPim, globalCells, 0, 0, +1.0, rows, cols, vals);
    appendScalarOpTriplets(B1Pre, globalCells, 0, N, -1.0, rows, cols, vals);
    appendScalarOpTriplets(B2Pre, globalCells, 0, N, -1.0, rows, cols, vals);

    appendScalarOpTriplets(B1Pim, globalCells, N, 0, +1.0, rows, cols, vals);
    appendScalarOpTriplets(B2Pim, globalCells, N, 0, +1.0, rows, cols, vals);
    appendScalarOpTriplets(AopPre, globalCells, N, N, +1.0, rows, cols, vals);
}
