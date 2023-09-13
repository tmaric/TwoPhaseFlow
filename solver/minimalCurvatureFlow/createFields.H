#include "createRDeltaT.H"

Info<< "Reading field psi\n" << endl;
volScalarField psi 
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        U.time().constant(),
        U.db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Minimal curvature flow
autoPtr<surfaceTensionForceModel> surfTensionForce = surfaceTensionForceModel::New
(
    transportProperties.subDict("surfaceForces"),
    alpha1,
    phi,
    U
);

volScalarField Ktilde 
(
    IOobject
    (
        "Ktilde",
        runTime.timeName(),
        mesh, 
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    surfTensionForce->K()
);

// These initial values are not used.
surfaceScalarField gradKtildeFlux ("gradKtildeFlux", fvc::interpolate(fvc::grad(Ktilde)) & mesh.Sf());

dimensionedScalar lambdaf ("lambdaf", dimVol / dimTime, 1.);

surfaceScalarField alphaCo(
    IOobject
    (
        "alphaCo",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(
        "alphaCo", 
        dimless, 
        runTime.controlDict().get<scalar>("maxAlphaCo")
    )
);

volScalarField divError
(
    IOobject
    (
        "divError",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(phi)
);

autoPtr<advectionSchemes> advector
(
    advectionSchemes::New(alpha1,phi,U)
);





