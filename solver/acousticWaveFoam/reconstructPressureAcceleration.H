// Cell value of (1/rho) grad(p) reconstructed from the face fluxes of the
// pressure Laplacian, invRhof*snGrad(p)*magSf, as in fvc::reconstruct.
//
// Pure cells weight their faces equally. A mixed cell contains a kink of the
// velocity, which is far steeper in the gas than in the liquid, so an equal
// weighting overestimates its value. There each face is weighted by its
// density 1/invRhof, and the cell value follows the phase that carries the
// momentum of the cell.

//
// A non-orthogonal correction applied by the solver (nonOrthogonalCorrection)
// is not part of the uncorrected snGrad; its face flux is then passed as
// correctionFlux and added, so that the face fluxes stay consistent.

Foam::tmp<Foam::volVectorField> reconstructPressureAcceleration
(
    const Foam::volScalarField& p,
    const Foam::surfaceScalarField& invRhof,
    const Foam::volScalarField& alpha1,
    const Foam::scalarField* correctionFlux = nullptr
)
{
    using namespace Foam;

    const fvMesh& mesh = p.mesh();
    const scalar mixedTol = 1e-12;

    surfaceScalarField faceFlux
    (
        invRhof*fvc::snGrad(p)*mesh.magSf()
    );
    if (correctionFlux)
    {
        faceFlux.primitiveFieldRef() += *correctionFlux;
    }

    const auto weight = [&](const label celli, const scalar faceInvRho)
    {
        const scalar alpha = alpha1[celli];
        return (alpha > mixedTol && alpha < 1 - mixedTol)
            ? 1/faceInvRho
            : scalar(1);
    };

    tensorField normalSum(mesh.nCells(), Zero);
    vectorField fluxSum(mesh.nCells(), Zero);

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    forAll(owner, facei)
    {
        const vector nHat = Sf[facei]/magSf[facei];
        const tensor nS = nHat*Sf[facei];
        const vector nFlux = nHat*faceFlux[facei];

        const scalar wOwn = weight(owner[facei], invRhof[facei]);
        const scalar wNei = weight(neighbour[facei], invRhof[facei]);

        normalSum[owner[facei]] += wOwn*nS;
        fluxSum[owner[facei]] += wOwn*nFlux;
        normalSum[neighbour[facei]] += wNei*nS;
        fluxSum[neighbour[facei]] += wNei*nFlux;
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& faceCells = mesh.boundary()[patchi].faceCells();
        const vectorField& pSf = Sf.boundaryField()[patchi];
        const scalarField& pMagSf = magSf.boundaryField()[patchi];
        const scalarField& pFlux = faceFlux.boundaryField()[patchi];
        const scalarField& pInvRho = invRhof.boundaryField()[patchi];

        forAll(faceCells, facei)
        {
            const label celli = faceCells[facei];
            const vector nHat = pSf[facei]/pMagSf[facei];
            const scalar w = weight(celli, pInvRho[facei]);

            normalSum[celli] += w*(nHat*pSf[facei]);
            fluxSum[celli] += w*(nHat*pFlux[facei]);
        }
    }

    tmp<volVectorField> tacceleration
    (
        new volVectorField
        (
            IOobject
            (
                "reconstructedPressureAcceleration",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedVector(dimVelocity/dimTime, Zero),
            calculatedFvPatchVectorField::typeName
        )
    );

    // inv() of a tensor field also handles the empty directions of 1D and
    // 2D meshes, as in fvc::reconstruct.
    tacceleration.ref().primitiveFieldRef() = inv(normalSum) & fluxSum;

    return tacceleration;
}
