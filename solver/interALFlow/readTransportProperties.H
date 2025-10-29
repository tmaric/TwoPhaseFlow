/*    Info<< "Reading transportProperties\n" << endl;

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
*/
/*    dimensionedScalar mu
    (
//        transportProperties.lookup("nu")
     "mu",
     dimViscosity*dimDensity,
     transportProperties
    );*/

/*     dimensionedScalar rhog
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
*/
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

/*     dimensionedScalar f
     (
     //    thermodynamicProperties.lookup("psi")
      "f",
      dimless/dimTime,
      transportProperties
     );*/

/*
dimensionedScalar p0
(
    "p0",
    dimPressure,
    transportProperties
);
*/

/*dimensionedScalar R
(
      //    thermodynamicProperties.lookup("psi")
        "R",
       dimLength,
       transportProperties
);*/
