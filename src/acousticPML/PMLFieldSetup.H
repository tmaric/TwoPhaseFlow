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

    // Axisymmetric wedge: azimuthal stretch of the cylindrical PML
    bool wedgeStretch;
    direction axialDir;
    direction radialDir;
    direction azimuthalDir;
};


// Coordinate direction of a wedge vector; the wedge must be aligned with the
// coordinate axes.
inline direction wedgeCoordinateDirection(const vector& v, const word& name)
{
    direction dir = 0;
    for (direction cmpt = 1; cmpt < vector::nComponents; ++cmpt)
    {
        if (mag(v.component(cmpt)) > mag(v.component(dir)))
        {
            dir = cmpt;
        }
    }

    if (mag(v.component(dir)) < 1 - 1e-6*mag(v))
    {
        FatalErrorInFunction
            << "wedgeAzimuthalStretch requires the wedge " << name << ' ' << v
            << " to be aligned with a coordinate axis."
            << exit(FatalError);
    }

    return dir;
}

inline Config readConfig
(
    const fvMesh& mesh,
    const IOdictionary& transportProperties
)
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

    // In an axisymmetric wedge the Cartesian tensor reproduces the exact
    // radial stretch r~ = r + (1/(i omega)) int sigma_r dr only together with
    // the azimuthal stretch s_phi = r~/r, i.e.
    // sigma_phi(r) = (1/r) int_0^r sigma_r(r') dr'. The wedge axis must pass
    // through the origin.
    cfg.wedgeStretch =
        transportProperties.getOrDefault<bool>("wedgeAzimuthalStretch", false);
    cfg.axialDir = 0;
    cfg.radialDir = 0;
    cfg.azimuthalDir = 0;

    if (cfg.wedgeStretch)
    {
        const wedgePolyPatch* wedgePtr = nullptr;
        for (const polyPatch& patch : mesh.boundaryMesh())
        {
            if (isA<wedgePolyPatch>(patch))
            {
                wedgePtr = &refCast<const wedgePolyPatch>(patch);
                break;
            }
        }

        if (!wedgePtr)
        {
            FatalErrorInFunction
                << "wedgeAzimuthalStretch is enabled but the mesh has no "
                << "wedge patch."
                << exit(FatalError);
        }

        cfg.axialDir = wedgeCoordinateDirection(wedgePtr->axis(), "axis");
        cfg.azimuthalDir =
            wedgeCoordinateDirection(wedgePtr->centreNormal(), "centre normal");
        cfg.radialDir = 3 - cfg.axialDir - cfg.azimuthalDir;

        if
        (
            cfg.PMLMin.component(cfg.radialDir) + cfg.L > 0
         && mag(cfg.sigmaMaxVector.component(cfg.radialDir)) > 0
        )
        {
            FatalErrorInFunction
                << "wedgeAzimuthalStretch does not support radial damping "
                << "next to the axis: set PMLMin component "
                << label(cfg.radialDir) << " below -L."
                << exit(FatalError);
        }

        const scalar span = mag(mesh.bounds().span());
        if (gMin(mesh.points().component(cfg.radialDir)) < -1e-9*span)
        {
            FatalErrorInFunction
                << "wedgeAzimuthalStretch assumes that the wedge axis passes "
                << "through the origin and the mesh lies on its positive "
                << "radial side."
                << exit(FatalError);
        }
    }

    return cfg;
}


// Primitive of the radial damping profile, F(r) = int^r sigma_r(r') dr'.
inline scalar radialDampingIntegral(const scalar r, const Config& cfg)
{
    const scalar rMax = cfg.PMLMax.component(cfg.radialDir);
    const scalar s = max((r - (rMax - cfg.L))/cfg.L, scalar(0));

    return cfg.sigmaMaxVector.component(cfg.radialDir)
       *(
            cfg.L*Foam::pow(min(s, scalar(1)), cfg.po + 1)/(cfg.po + 1)
          + max(r - rMax, scalar(0))
        );
}

inline tensor evaluateSigma(const point& pt, const Config& cfg)
{
    vector damping(Zero);

    // Wedge: the radial profile is evaluated at the distance from the axis
    point x(pt);
    scalar r = 0;
    if (cfg.wedgeStretch)
    {
        r = Foam::sqrt
        (
            sqr(pt.component(cfg.radialDir))
          + sqr(pt.component(cfg.azimuthalDir))
        );
        x.component(cfg.radialDir) = r;
    }

    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        if (cfg.wedgeStretch && cmpt == cfg.azimuthalDir)
        {
            continue;
        }

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

    if (cfg.wedgeStretch && r > VSMALL)
    {
        damping.component(cfg.azimuthalDir) =
            (radialDampingIntegral(r, cfg) - radialDampingIntegral(0, cfg))/r;
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
