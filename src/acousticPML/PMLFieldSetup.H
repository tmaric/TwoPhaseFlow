// Shared PML coefficient setup for the frequency-domain acoustic solvers.

#include "fvCFD.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "processorPolyPatch.H"

namespace Foam
{
namespace acousticPML
{

enum PMLType
{
    RECTANGLE,
    SPHERICAL,
    UNDEFINED
};

struct Config
{
    scalar cg;
    scalar cl;
    scalar omega;
    scalar L;
    scalar sigmaMax;
    vector sigmaMaxVector;
    scalar po;
    PMLType type;
    vector PMLMax;
    vector PMLMin;
    scalar Rmin;
    scalar Rmax;
    vector origin;
};

inline Config readConfig(const IOdictionary& transportProperties)
{
    Config cfg;

    const dimensionedScalar cg("cg", dimVelocity, transportProperties);
    const dimensionedScalar cl("cl", dimVelocity, transportProperties);
    const dimensionedScalar f("f", dimless/dimTime, transportProperties);
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

    cfg.cg = cg.value();
    cfg.cl = cl.value();
    cfg.omega = constant::mathematical::twoPi*f.value();
    cfg.L = L.value();
    cfg.sigmaMax = sigmaMax.value();
    cfg.sigmaMaxVector = transportProperties.lookupOrDefault<vector>
    (
        "sigmaMaxVector",
        vector(cfg.sigmaMax, cfg.sigmaMax, cfg.sigmaMax)
    );
    cfg.po = po.value();
    cfg.type = UNDEFINED;
    cfg.PMLMax = Zero;
    cfg.PMLMin = Zero;
    cfg.Rmin = 0.0;
    cfg.Rmax = 0.0;
    cfg.origin = Zero;

    const word PMLTypeWord = transportProperties.get<word>("PMLType");

    if (PMLTypeWord == "rectangle")
    {
        cfg.type = RECTANGLE;
        dimensionedVector PMLMax;
        dimensionedVector PMLMin;
        transportProperties.lookup("PMLMax") >> PMLMax;
        transportProperties.lookup("PMLMin") >> PMLMin;
        cfg.PMLMax = PMLMax.value();
        cfg.PMLMin = PMLMin.value();
    }
    else if (PMLTypeWord == "spherical")
    {
        cfg.type = SPHERICAL;
        dimensionedScalar Rmin("Rmin", dimLength, transportProperties);
        dimensionedScalar Rmax("Rmax", dimLength, transportProperties);
        dimensionedVector origin;
        transportProperties.lookup("origin") >> origin;
        cfg.Rmin = Rmin.value();
        cfg.Rmax = Rmax.value();
        cfg.origin = origin.value();
    }
    else
    {
        FatalErrorInFunction
            << "Unknown PMLType '" << PMLTypeWord << "' in transportProperties."
            << nl << "Valid types are: rectangle, spherical."
            << exit(FatalError);
    }

    return cfg;
}

inline void evaluatePoint
(
    const point& x,
    const Config& cfg,
    const scalar invSoundSpeedSqr,
    scalar& sigmaLocal,
    scalar& sc0,
    scalar& sc1,
    scalar& tc1,
    tensor& t0,
    tensor& t1
)
{
    sigmaLocal = 0.0;
    sc0 = 0.0;
    sc1 = 0.0;
    tc1 = 0.0;
    t0 = Zero;
    t1 = Zero;

    switch (cfg.type)
    {
        case RECTANGLE:
        {
            vector DpFunc(Zero);

            for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
            {
                const scalar distToMax = cfg.PMLMax.component(cmpt) - x.component(cmpt);
                const scalar distToMin = x.component(cmpt) - cfg.PMLMin.component(cmpt);

                if (distToMax < cfg.L || distToMin < cfg.L)
                {
                    const scalar distCmpt = (distToMax < cfg.L) ? distToMax : distToMin;
                    scalar s = (cfg.L - distCmpt)/cfg.L;
                    s = min(max(s, scalar(0)), scalar(1));
                    DpFunc.component(cmpt) =
                        cfg.sigmaMaxVector.component(cmpt)*Foam::pow(s, cfg.po);
                }
            }

            const scalar dx = DpFunc.x();
            const scalar dy = DpFunc.y();
            const scalar dz = DpFunc.z();
            const scalar dxyzOverOmega =
                (mag(cfg.omega) > VSMALL) ? dx*dy*dz/cfg.omega : 0.0;

            // Real and imaginary parts of diag(etaY*etaZ/etaX, ...),
            // excluding the identity contribution.
            const scalar t0x =
                (dx*(dy + dz - dx) - dy*dz)
               /(sqr(cfg.omega) + sqr(dx) + VSMALL);
            const scalar t0y =
                (dy*(dx + dz - dy) - dx*dz)
               /(sqr(cfg.omega) + sqr(dy) + VSMALL);
            const scalar t0z =
                (dz*(dx + dy - dz) - dx*dy)
               /(sqr(cfg.omega) + sqr(dz) + VSMALL);

            const scalar t1x =
                (dxyzOverOmega + cfg.omega*(dy + dz - dx))
               /(sqr(cfg.omega) + sqr(dx) + VSMALL);
            const scalar t1y =
                (dxyzOverOmega + cfg.omega*(dx + dz - dy))
               /(sqr(cfg.omega) + sqr(dy) + VSMALL);
            const scalar t1z =
                (dxyzOverOmega + cfg.omega*(dx + dy - dz))
               /(sqr(cfg.omega) + sqr(dz) + VSMALL);

            t0 = tensor(t0x, 0, 0, 0, t0y, 0, 0, 0, t0z);
            t1 = tensor(t1x, 0, 0, 0, t1y, 0, 0, 0, t1z);

            // Retain TC1 for legacy consumers; directional solvers use T1.
            tc1 = (t1x + t1y + t1z)/3.0;

            sigmaLocal = dx;
            if (mag(dy) > mag(sigmaLocal))
            {
                sigmaLocal = dy;
            }
            if (mag(dz) > mag(sigmaLocal))
            {
                sigmaLocal = dz;
            }
            sc0 = (dx*dy + dx*dz + dy*dz)*invSoundSpeedSqr;
            sc1 =
                (cfg.omega*(dx + dy + dz) - dxyzOverOmega)
               *invSoundSpeedSqr;
            break;
        }

        case SPHERICAL:
        {
            const scalar dist = mag(x - cfg.origin);

            if ((cfg.Rmax - dist) < cfg.L)
            {
                scalar s = (cfg.L - (cfg.Rmax - dist))/cfg.L;
                s = min(max(s, scalar(0)), scalar(1));
                sigmaLocal = cfg.sigmaMax*Foam::pow(s, cfg.po);
            }

            const scalar denom = sqr(sigmaLocal) + sqr(cfg.omega);
            const scalar sigmaCubeOverOmega =
                (mag(cfg.omega) > VSMALL)
              ? (sigmaLocal*sigmaLocal*sigmaLocal/cfg.omega)
              : 0.0;

            tc1 =
                (sigmaLocal > 0 || cfg.omega > 0)
              ? (sigmaCubeOverOmega + cfg.omega*sigmaLocal)/(denom + VSMALL)
              : 0.0;
            sc0 = (3*sigmaLocal*sigmaLocal)*invSoundSpeedSqr;
            sc1 =
                (cfg.omega*(3*sigmaLocal) - sigmaCubeOverOmega)
               *invSoundSpeedSqr;
            t1 = tensor(tc1, 0, 0, 0, tc1, 0, 0, 0, tc1);
            break;
        }

        case UNDEFINED:
        {
            FatalErrorInFunction
                << "Internal error: PML type is undefined."
                << exit(FatalError);
            break;
        }
    }
}

inline void updateFields
(
    const fvMesh& mesh,
    const Config& cfg,
    const volScalarField& alpha,
    volScalarField& SC0,
    volScalarField& SC1,
    volScalarField& TC1,
    volTensorField& T0,
    volTensorField& T1,
    volScalarField& sigma
)
{
    const auto& Ccs = mesh.C();

    forAll(Ccs, celli)
    {
        const scalar alphaCell = min(max(alpha[celli], scalar(0)), scalar(1));
        const scalar invSoundSpeedSqr =
            alphaCell/sqr(cfg.cl) + (1 - alphaCell)/sqr(cfg.cg);

        scalar sigmaLocal = 0.0;
        scalar sc0 = 0.0;
        scalar sc1 = 0.0;
        scalar tc1 = 0.0;
        tensor t0(Zero);
        tensor t1(Zero);
        evaluatePoint
        (
            Ccs[celli],
            cfg,
            invSoundSpeedSqr,
            sigmaLocal,
            sc0,
            sc1,
            tc1,
            t0,
            t1
        );

        sigma[celli] = sigmaLocal;
        SC0[celli] = sc0;
        SC1[celli] = sc1;
        TC1[celli] = tc1;
        T0[celli] = t0;
        T1[celli] = t1;
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const surfaceVectorField& Cf = mesh.Cf();

    forAll(patches, patchi)
    {
        sigma.boundaryFieldRef()[patchi] = 0.0;
        SC0.boundaryFieldRef()[patchi] = 0.0;
        SC1.boundaryFieldRef()[patchi] = 0.0;
        TC1.boundaryFieldRef()[patchi] = 0.0;
        T0.boundaryFieldRef()[patchi] = tensor::zero;
        T1.boundaryFieldRef()[patchi] = tensor::zero;

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
            const scalar alphaFace = min
            (
                max(alpha.boundaryField()[patchi][facei], scalar(0)),
                scalar(1)
            );
            const scalar invSoundSpeedSqr =
                alphaFace/sqr(cfg.cl) + (1 - alphaFace)/sqr(cfg.cg);

            scalar sigmaLocal = 0.0;
            scalar sc0 = 0.0;
            scalar sc1 = 0.0;
            scalar tc1 = 0.0;
            tensor t0(Zero);
            tensor t1(Zero);
            evaluatePoint
            (
                Cf.boundaryField()[patchi][facei],
                cfg,
                invSoundSpeedSqr,
                sigmaLocal,
                sc0,
                sc1,
                tc1,
                t0,
                t1
            );

            sigma.boundaryFieldRef()[patchi][facei] = sigmaLocal;
            SC0.boundaryFieldRef()[patchi][facei] = sc0;
            SC1.boundaryFieldRef()[patchi][facei] = sc1;
            TC1.boundaryFieldRef()[patchi][facei] = tc1;
            T0.boundaryFieldRef()[patchi][facei] = t0;
            T1.boundaryFieldRef()[patchi][facei] = t1;
        }
    }

    sigma.correctBoundaryConditions();
    SC0.correctBoundaryConditions();
    SC1.correctBoundaryConditions();
    TC1.correctBoundaryConditions();
    T0.correctBoundaryConditions();
    T1.correctBoundaryConditions();
}

} // End namespace acousticPML
} // End namespace Foam
