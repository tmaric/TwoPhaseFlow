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

volScalarField SC0
(
    IOobject
    (
        "SC0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField SC1
(
    IOobject
    (
        "SC1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField TC1
(
    IOobject
    (
        "TC1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField sigma
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
