// #include "readThermodynamicProperties.H"
#include "readTransportProperties.H"
    Info<< "Reading field p\n" << endl;
    volScalarField pa
    (
        IOobject
        (
            "pa",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

pa == pa*(1-alpha1);

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
volVectorField Ua
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
    dimensionedVector("Ua", dimVelocity, vector::zero)
);

/*
volScalarField AcousticForcePotential
    (
       IOobject
        (
            "AcousticForcePotential",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
         dimensionedScalar("AcousticForcePotential", dimForce*dimLength,
         0) 
    );*/

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

volScalarField prAver
     (
       IOobject
        (
            "prAver",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	pr
    );

volTensorField momFluxAver
(
      IOobject
      (
          "momFluxAver",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
       ),
       momFlux
 );

mesh.setFluxRequired(pa.name()); 
