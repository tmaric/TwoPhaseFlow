// #include "readThermodynamicProperties.H"
#include "readTransportProperties.H"
    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField pDot
    (
        IOobject
        (
            "pDot",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("pDot", dimPressure/dimTime, Zero)
    );

	Info<< "Reading alpha field pr\n" << endl;
	volScalarField alpha1
	(
	    IOobject
	    (
	        "alpha.water",
	        runTime.timeName(),
	        mesh,
	        IOobject::READ_IF_PRESENT,
	        IOobject::NO_WRITE
	    ),
	    mesh,
	    dimensionedScalar(dimless, Zero)
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
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("U", dimVelocity, vector::zero)
);

// volScalarField rho
// (
//     IOobject
//     (
//         "rho",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//    rhol*alpha1 + rhog*(1-alpha1)
// ); 

volScalarField invc2
(
    IOobject
    (
        "invc2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1/sqr(cl) + (1 - alpha1)/sqr(cg)
);

volScalarField rho
(
    IOobject
    (
        "rho",
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
    );
*/
// radiation force pressure
volScalarField prT
     (
       IOobject
        (
            "prT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
         dimensionedScalar("prT", dimPressure,
         0) 
    );

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
         dimensionedScalar("pr", dimPressure,
         0) 
    );

volTensorField momFluxT
(
      IOobject
      (
          "momFluxT",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
       ),
       mesh,
       dimensionedTensor("momFluxT", dimDensity*dimVelocity*dimVelocity, tensor::zero) 
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

#   include "createPhi.H"

    mesh.setFluxRequired(p.name()); //mesh.schemesDict().setFluxRequired(p.name());

#include "createFvOptions.H"

Foam::advection::isoAdvection advector(alpha1, phi, U);

cutFaceAdvect cutFace(alpha1.mesh(), alpha1);
