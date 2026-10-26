// Cell value of (1/rho) grad(p) reconstructed from the face fluxes of the
// pressure Laplacian, invRhof*snGrad(p)*magSf, as in fvc::reconstruct.
//
// Pure cells weight their faces equally. A mixed cell contains a kink of the
// velocity, which is far steeper in the gas than in the liquid, so an equal
// weighting overestimates its value. There each face is weighted by its
// density 1/invRhof, and the cell value follows the phase that carries the
// momentum of the cell.
//
// The face flux is passed in so that it is the flux of the discrete operator
// that was solved, including its non-orthogonal correction; see
// laplacianFaceFlux. Same cell algorithm as
// solver/acousticWaveFoam/reconstructPressureAcceleration.H.
//
// Boundary values: coupled patches are evaluated by their patch type. On the
// other patches the cell value is extrapolated and its normal component is
// replaced by the face flux, as fvc::grad does with the face gradient; at a
// sound-hard wall the normal component is therefore zero.

#ifndef reconstructPressureAcceleration_H
#define reconstructPressureAcceleration_H

#include "fvCFD.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "surfaceInterpolationScheme.H"
#include "snGradScheme.H"

namespace Foam
{
namespace acoustic
{

// Face flux gamma*snGrad(p)*|S| of fvm::laplacian(gamma, p), with the snGrad
// scheme of its laplacian scheme entry (Gauss <interpolation> <snGrad>), so
// that the explicit non-orthogonal correction of the solved operator is
// included. Equal to fvMatrix::flux() of that operator, which needs a
// flux-required field and then cannot be scaled by a cell field.
inline tmp<surfaceScalarField> laplacianFaceFlux
(
    const surfaceScalarField& gamma,
    const volScalarField& p
)
{
    const fvMesh& mesh = p.mesh();
    const word name("laplacian(" + gamma.name() + ',' + p.name() + ')');
    ITstream& is = mesh.laplacianScheme(name);

    const word scheme(is);
    if (scheme != "Gauss")
    {
        FatalErrorInFunction
            << "Face-flux velocity recovery needs a Gauss laplacian scheme for "
            << name << ", not " << scheme << exit(FatalError);
    }

    // The interpolation entry is not used for a face coefficient; it is read
    // to reach the snGrad entry.
    surfaceInterpolationScheme<scalar>::New(mesh, is);
    tmp<fv::snGradScheme<scalar>> tsnGrad(fv::snGradScheme<scalar>::New(mesh, is));

    return gamma*tsnGrad().snGrad(p)*mesh.magSf();
}

inline tmp<volVectorField> reconstructPressureAcceleration
(
    const surfaceScalarField& faceFlux,
    const surfaceScalarField& invRhof,
    const volScalarField& alpha1,
    const scalar mixedTol
)
{
    const fvMesh& mesh = faceFlux.mesh();

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
            dimensionedVector(faceFlux.dimensions()/dimArea, Zero),
            extrapolatedCalculatedFvPatchVectorField::typeName
        )
    );
    volVectorField& acceleration = tacceleration.ref();

    // inv() of a tensor field also handles the empty directions of 1D and
    // 2D meshes, as in fvc::reconstruct.
    acceleration.primitiveFieldRef() = inv(normalSum) & fluxSum;
    acceleration.correctBoundaryConditions();

    forAll(acceleration.boundaryField(), patchi)
    {
        fvPatchVectorField& pa = acceleration.boundaryFieldRef()[patchi];
        if (pa.coupled() || pa.empty())
        {
            continue;
        }
        const vectorField nHat(Sf.boundaryField()[patchi]/magSf.boundaryField()[patchi]);
        const scalarField faceNormal
        (
            faceFlux.boundaryField()[patchi]/magSf.boundaryField()[patchi]
        );
        pa == pa + nHat*(faceNormal - (nHat & pa));
    }

    return tacceleration;
}

// Cells in which the cell-centred Gauss gradient divided by the cell density
// is inconsistent with the face fluxes: mixed cells and cells with a face
// whose volume fraction differs from their own (pure cells beside the
// interface).
inline boolList interfaceVelocityCells
(
    const volScalarField& alpha1,
    const surfaceScalarField& alphaf,
    const scalar tol
)
{
    const fvMesh& mesh = alpha1.mesh();
    boolList cells(mesh.nCells(), false);

    forAll(alpha1, celli)
    {
        cells[celli] = alpha1[celli] > tol && alpha1[celli] < 1 - tol;
    }

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    forAll(owner, facei)
    {
        if (mag(alphaf[facei] - alpha1[owner[facei]]) > tol)
        {
            cells[owner[facei]] = true;
        }
        if (mag(alphaf[facei] - alpha1[neighbour[facei]]) > tol)
        {
            cells[neighbour[facei]] = true;
        }
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& faceCells = mesh.boundary()[patchi].faceCells();
        const scalarField& pAlphaf = alphaf.boundaryField()[patchi];
        forAll(faceCells, facei)
        {
            if (mag(pAlphaf[facei] - alpha1[faceCells[facei]]) > tol)
            {
                cells[faceCells[facei]] = true;
            }
        }
    }

    return cells;
}

} // End namespace acoustic
} // End namespace Foam

#endif
