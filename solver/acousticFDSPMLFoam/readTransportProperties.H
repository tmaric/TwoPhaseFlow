    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

     dimensionedScalar rhog
     (
 //        transportProperties.lookup("nu")
      "rhog",
      dimDensity,
      transportProperties
     );

     dimensionedScalar rhol
     (
 //        transportProperties.lookup("nu")
      "rhol",
      dimDensity,
      transportProperties
     );

     dimensionedScalar kg
     (
     //    thermodynamicProperties.lookup("psi")
      "kg",
      dimCompressibility/dimDensity,
      transportProperties
     );
     dimensionedScalar kl
     (
     //    thermodynamicProperties.lookup("psi")
      "kl",
      dimCompressibility/dimDensity,
      transportProperties
     );

     dimensionedScalar cg
     (
     //    thermodynamicProperties.lookup("psi")
      "cg",
      dimVelocity,
      transportProperties
     );

     dimensionedScalar f
     (
     //    thermodynamicProperties.lookup("psi")
      "f",
      dimless/dimTime,
      transportProperties
     );

Info<< "Creating PML coefficients fields \n" << endl;
//PML thickness
    dimensionedScalar L
    (
      "L", dimless, transportProperties
    );
//PML max x, y ,z
    dimensionedVector PMLMax
    (
      "PMLMax", dimless, transportProperties
    );
//PML min x, y, z
    dimensionedVector PMLMin
    (
      "PMLMin", dimless, transportProperties
    );
