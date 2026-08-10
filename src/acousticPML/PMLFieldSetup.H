// Shared geometry and damping-profile setup for acoustic PML fields.

#include "fvCFD.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "processorPolyPatch.H"

namespace Foam
{
namespace acousticPML
{

struct Config
{
    scalar L;
    scalar sigmaMax;
    vector sigmaMaxVector;
    scalar po;
    vector PMLMax;
    vector PMLMin;
};

inline Config readConfig(const IOdictionary& transportProperties)
{
    Config cfg;

    const dimensionedScalar L("L", dimless, transportProperties);
    const dimensionedScalar sigmaMax
    (
        "sigmaMax",
        dimless/dimTime,
        transportProperties.lookupOrDefault<scalar>("sigmaMax", 0.0)
    );
    const dimensionedScalar po
    (
        "po",
        dimless,
        transportProperties.lookupOrDefault<scalar>
        (
            "po",
            transportProperties.lookupOrDefault<scalar>("p", 2.0)
        )
    );

    cfg.L = L.value();
    cfg.sigmaMax = -sigmaMax.value();
    cfg.sigmaMaxVector = transportProperties.lookupOrDefault<vector>
    (
        "sigmaMaxVector",
        vector(cfg.sigmaMax, cfg.sigmaMax, cfg.sigmaMax)
    );
    cfg.po = po.value();
    cfg.PMLMax = Zero;
    cfg.PMLMin = Zero;

    const word PMLTypeWord = transportProperties.get<word>("PMLType");

    if (PMLTypeWord != "rectangle")
    {
        FatalErrorInFunction
            << "Unknown PMLType '" << PMLTypeWord << "' in transportProperties."
            << nl << "The supported type is: rectangle."
            << exit(FatalError);
    }

    dimensionedVector PMLMax;
    dimensionedVector PMLMin;
    transportProperties.lookup("PMLMax") >> PMLMax;
    transportProperties.lookup("PMLMin") >> PMLMin;
    cfg.PMLMax = PMLMax.value();
    cfg.PMLMin = PMLMin.value();

    return cfg;
}

inline tensor evaluateSigma(const point& x, const Config& cfg)
{
    vector damping(Zero);

    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        const scalar distToMax =
            cfg.PMLMax.component(cmpt) - x.component(cmpt);
        const scalar distToMin =
            x.component(cmpt) - cfg.PMLMin.component(cmpt);

        if (distToMax < cfg.L || distToMin < cfg.L)
        {
            const scalar dist =
                (distToMax < cfg.L) ? distToMax : distToMin;
            const scalar s = min
            (
                max((cfg.L - dist)/cfg.L, scalar(0)),
                scalar(1)
            );
            damping.component(cmpt) =
                cfg.sigmaMaxVector.component(cmpt)*Foam::pow(s, cfg.po);
        }
    }

    return tensor
    (
        damping.x(), 0, 0,
        0, damping.y(), 0,
        0, 0, damping.z()
    );
}

inline void updateSigma
(
    const fvMesh& mesh,
    const Config& cfg,
    volTensorField& sigma
)
{
    forAll(mesh.C(), celli)
    {
        sigma[celli] = evaluateSigma(mesh.C()[celli], cfg);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const surfaceVectorField& Cf = mesh.Cf();

    forAll(patches, patchi)
    {
        sigma.boundaryFieldRef()[patchi] = tensor::zero;

        const polyPatch& patch = patches[patchi];
        if
        (
            isA<emptyPolyPatch>(patch)
         || isA<wedgePolyPatch>(patch)
         || isA<processorPolyPatch>(patch)
        )
        {
            continue;
        }

        forAll(patch, facei)
        {
            sigma.boundaryFieldRef()[patchi][facei] =
                evaluateSigma(Cf.boundaryField()[patchi][facei], cfg);
        }
    }

    sigma.correctBoundaryConditions();
}

} // End namespace acousticPML
} // End namespace Foam
