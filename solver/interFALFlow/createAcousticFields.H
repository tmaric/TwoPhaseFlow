#include "readTransportProperties.H"

    Info<< "Reading field Pre and Pim\n" << endl;
    volScalarField Pre
    (
        IOobject
        (
            "Pre",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField Pim
    (
        IOobject
        (
            "Pim",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	surfaceScalarField alphaf
	(
	    IOobject
	    (
	        "alphaf",
	        runTime.timeName(),
	        mesh,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	     ), 
	     fvc::interpolate(alpha1)
	);

volScalarField rhoAcoustic
(
    IOobject
    (
        "rhoAcoustic",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    1/(alpha1/rhol + (1 - alpha1)/rhog)
);

surfaceScalarField invRhof
(
    IOobject
    (
        "invRhof",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alphaf/rhol + (1 - alphaf)/rhog
);

// radiation force pressure
volScalarField pr
     (
       IOobject
        (
            "pr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
         dimensionedScalar("pr", dimPressure, 0) 
    );

// Face mask: beta=1 on faces touching an interface cell, beta=0 elsewhere.
surfaceScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("beta", dimless, 0.0)
);

volScalarField Pa
(
    IOobject
    (
        "Pa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Pa", dimPressure, 0)
);

volScalarField Ua
(
    IOobject
    (
        "Ua",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Ua", dimVelocity, 0)
);

volTensorField momFlux
(
      IOobject
      (
          "momFlux",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
       ),
       mesh,
       dimensionedTensor("momFlux", dimDensity*dimVelocity*dimVelocity, tensor::zero) 
 );

volScalarField pMomFlux
(
    IOobject
    (
        "pMomFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("pMomFlux", dimPressure, 0)
);

volScalarField C0
(
    IOobject
    (
        "C0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("C0", dimless/sqr(dimLength), Zero)
);

volScalarField C1
(
    IOobject
    (
        "C1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("C1", dimless/sqr(dimLength), Zero)
);

volTensorField T0
(
    IOobject
    (
        "T0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("T0", dimless, tensor::zero)
);

volTensorField T1
(
    IOobject
    (
        "T1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("T1", dimless, tensor::zero)
);

volTensorField sigma
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
