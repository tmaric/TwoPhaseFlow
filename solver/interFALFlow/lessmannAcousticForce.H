// Build a Lessmann-style shifted gas-side acoustic interface force:
//   fac = (pShift*nOut + momFluxShift·nOut) |grad(alpha1)|
//
// Stage 1 implementation:
// - cell-centred normals from grad(alpha1)
// - donor-cell sampling at x + delta*nOut
// - fallback to a forward neighbour with low alpha1

static label selectLessmannDonorCell
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const scalar gasAlphaMax,
    const label celli,
    const point& samplePoint,
    const vector& nOut
)
{
    label donor = mesh.findCell(samplePoint);

    if (donor >= 0 && donor != celli && alpha1[donor] <= gasAlphaMax)
    {
        return donor;
    }

    const labelList& nbrs = mesh.cellCells()[celli];
    label bestNbr = -1;
    scalar bestAlpha = GREAT;
    scalar bestProj = -GREAT;

    forAll(nbrs, nbrI)
    {
        const label nei = nbrs[nbrI];
        const scalar proj = (mesh.C()[nei] - mesh.C()[celli]) & nOut;

        if (proj <= 0)
        {
            continue;
        }

        const scalar alphaNei = alpha1[nei];
        if
        (
            alphaNei <= gasAlphaMax
         && (alphaNei < bestAlpha || (alphaNei == bestAlpha && proj > bestProj))
        )
        {
            bestNbr = nei;
            bestAlpha = alphaNei;
            bestProj = proj;
        }
    }

    return bestNbr;
}

static void updateLessmannAcousticForce
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& pr,
    const volTensorField& momFlux,
    const scalar shiftCoeff,
    const scalar interfaceAlphaMin,
    const scalar interfaceAlphaMax,
    const scalar gasAlphaMax,
    const scalar minGradAlpha,
    volScalarField& prShift,
    volVectorField& acousticTraction,
    volVectorField& acousticForce
)
{
    const volVectorField gradAlpha(fvc::grad(alpha1));
    const scalarField& V = mesh.V();
    const vectorField& C = mesh.C();

    prShift = dimensionedScalar("zero", prShift.dimensions(), 0);
    acousticTraction = dimensionedVector("zero", acousticTraction.dimensions(), vector::zero);
    acousticForce = dimensionedVector("zero", acousticForce.dimensions(), vector::zero);

    forAll(alpha1, celli)
    {
        const scalar alphaCell = alpha1[celli];
        const vector gradAlphaCell = gradAlpha[celli];
        const scalar magGradAlpha = mag(gradAlphaCell);

        if
        (
            alphaCell <= interfaceAlphaMin
         || alphaCell >= interfaceAlphaMax
         || magGradAlpha <= minGradAlpha
        )
        {
            continue;
        }

        // alpha1=1 in liquid and alpha1=0 in gas, so -grad(alpha1) points gas-ward.
        const vector nOut = -gradAlphaCell/(magGradAlpha + VSMALL);
        const scalar localDelta =
            shiftCoeff*Foam::cbrt(Foam::max(V[celli], scalar(VSMALL)));
        const point samplePoint = C[celli] + localDelta*nOut;

        const label donor = selectLessmannDonorCell
        (
            mesh,
            alpha1,
            gasAlphaMax,
            celli,
            samplePoint,
            nOut
        );

        if (donor < 0)
        {
            continue;
        }

        const scalar sampledPr = pr[donor];
        const vector traction = sampledPr*nOut + (momFlux[donor] & nOut);

        prShift[celli] = sampledPr;
        acousticTraction[celli] = traction;
        acousticForce[celli] = traction*magGradAlpha;
    }

    prShift.correctBoundaryConditions();
    acousticTraction.correctBoundaryConditions();
    acousticForce.correctBoundaryConditions();
}
