// Block-assembly helpers for tensorLaplacianFoam.

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

    const lduAddressing& addr = L.lduAddr();
    const labelUList& owner = addr.lowerAddr();
    const labelUList& neigh = addr.upperAddr();

    forAll(diag, i)
    {
        PetscInt r = globalCells.toGlobal(i) + rowOffset;
        PetscInt c = globalCells.toGlobal(i) + colOffset;
        MatSetValue(M, r, c, scale*diag[i], ADD_VALUES);
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

            MatSetValue(M, ri, cj, scale*upperVal, ADD_VALUES);
            MatSetValue(M, rj, ci, scale*lowerVal, ADD_VALUES);
        }
    }
}

static scalarField sourceWithBoundary
(
    const fvScalarMatrix& M
)
{
    scalarField source(M.source());
    addBoundarySourceSimple(M, source);
    return source;
}

static scalar tensorFaceCoeff
(
    const vector& Sf,
    const tensor& Tf,
    const vector& d
)
{
    return (Sf & (Tf & d))/(magSqr(d) + VSMALL);
}

static void insertTensorLaplacianIntoBlock
(
    Mat& M,
    const volTensorField& T1,
    const volScalarField& psi,
    const globalIndex& globalCells,
    const PetscInt rowOffset,
    const PetscInt colOffset,
    const PetscScalar scale
)
{
    const fvMesh& mesh = psi.mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const vectorField& C = mesh.C();
    const vectorField& Sf = mesh.Sf();
    const tensorField& T1Cells = T1.primitiveField();

    forAll(neighbour, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        const vector d = C[nei] - C[own];
        const tensor Tf = 0.5*(T1Cells[own] + T1Cells[nei]);
        const scalar coeff = tensorFaceCoeff(Sf[facei], Tf, d);

        const PetscInt rowOwn = globalCells.toGlobal(own) + rowOffset;
        const PetscInt rowNei = globalCells.toGlobal(nei) + rowOffset;
        const PetscInt colOwn = globalCells.toGlobal(own) + colOffset;
        const PetscInt colNei = globalCells.toGlobal(nei) + colOffset;

        MatSetValue(M, rowOwn, colOwn, scale*(-coeff), ADD_VALUES);
        MatSetValue(M, rowOwn, colNei, scale*(+coeff), ADD_VALUES);
        MatSetValue(M, rowNei, colNei, scale*(-coeff), ADD_VALUES);
        MatSetValue(M, rowNei, colOwn, scale*(+coeff), ADD_VALUES);
    }

    const auto& psiB = psi.boundaryField();
    const auto& T1B = T1.boundaryField();
    const auto& CfB = mesh.Cf().boundaryField();
    const auto& SfB = mesh.Sf().boundaryField();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        const fvPatchField<scalar>& ptf = psiB[patchi];
        if (!ptf.size() || ptf.coupled() || !ptf.fixesValue())
        {
            continue;
        }

        const labelUList& faceCells = mesh.boundary()[patchi].faceCells();
        forAll(faceCells, facei)
        {
            const label own = faceCells[facei];
            const vector d = CfB[patchi][facei] - C[own];
            const scalar coeff = tensorFaceCoeff(SfB[patchi][facei], T1B[patchi][facei], d);
            const PetscInt rowOwn = globalCells.toGlobal(own) + rowOffset;
            const PetscInt colOwn = globalCells.toGlobal(own) + colOffset;
            MatSetValue(M, rowOwn, colOwn, scale*(-coeff), ADD_VALUES);
        }
    }
}

static scalarField tensorLaplacianBoundarySource
(
    const volTensorField& T1,
    const volScalarField& psi
)
{
    const fvMesh& mesh = psi.mesh();
    scalarField source(mesh.nCells(), Zero);

    const auto& psiB = psi.boundaryField();
    const auto& T1B = T1.boundaryField();
    const auto& CfB = mesh.Cf().boundaryField();
    const auto& SfB = mesh.Sf().boundaryField();
    const vectorField& C = mesh.C();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        const fvPatchField<scalar>& ptf = psiB[patchi];
        if (!ptf.size())
        {
            continue;
        }

        const labelUList& faceCells = mesh.boundary()[patchi].faceCells();
        const tmp<Field<scalar>> tpnf = ptf.coupled() ? ptf.patchNeighbourField() : tmp<Field<scalar>>();

        forAll(faceCells, facei)
        {
            if (!ptf.coupled() && !ptf.fixesValue())
            {
                continue;
            }

            const label own = faceCells[facei];
            const vector d = CfB[patchi][facei] - C[own];
            const scalar coeff = tensorFaceCoeff(SfB[patchi][facei], T1B[patchi][facei], d);
            const scalar boundaryValue = ptf.coupled() ? tpnf()[facei] : ptf[facei];
            source[own] += coeff*boundaryValue;
        }
    }

    return source;
}

static void assembleTensorBlockSystem
(
    Mat& M,
    const globalIndex& globalCells,
    const PetscInt N,
    fvScalarMatrix& AopPim,
    fvScalarMatrix& AopPre,
    fvScalarMatrix& B2Pre,
    fvScalarMatrix& B2Pim,
    const volTensorField& T1,
    const volScalarField& Pim,
    const volScalarField& Pre
)
{
    AopPre.boundaryManipulate(AopPre.psi().boundaryFieldRef());
    AopPim.boundaryManipulate(AopPim.psi().boundaryFieldRef());
    B2Pre.boundaryManipulate(B2Pre.psi().boundaryFieldRef());
    B2Pim.boundaryManipulate(B2Pim.psi().boundaryFieldRef());

    insertScalarOpIntoBlock(M, AopPim, globalCells, 0, 0, +1.0);
    insertTensorLaplacianIntoBlock(M, T1, Pre, globalCells, 0, N, -1.0);
    insertScalarOpIntoBlock(M, B2Pre, globalCells, 0, N, -1.0);

    insertTensorLaplacianIntoBlock(M, T1, Pim, globalCells, N, 0, +1.0);
    insertScalarOpIntoBlock(M, B2Pim, globalCells, N, 0, +1.0);
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
