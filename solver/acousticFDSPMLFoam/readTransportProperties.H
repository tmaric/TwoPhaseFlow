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

enum PMLType
{
    RECTANGLE,
    // SPHERICAL,
    CYLINDRICAL,
    UNDEFINED // A state for "not set" or "error"
};

const word PMLTypeWord = transportProperties.get<word>("PMLType");

PMLType selectedModel = UNDEFINED; // Default to undefined

if (PMLTypeWord == "rectangle")
{
    selectedModel = RECTANGLE;
}
/*else if (PMLTypeWord == "spherical")
{
    selectedModel = SPHERICAL;
}*/
else if (PMLTypeWord == "cylindrical")
{
    selectedModel = CYLINDRICAL;
}
else
{
    FatalErrorInFunction
        << "Unknown PMLType '" << PMLTypeWord
        << "' in transportProperties." << nl
        << "Valid types are: rectangle, spherical, cylindrical."
        << exit(FatalError);
}

// PML parameters (resolved by model)
dimensionedVector PMLMax;
dimensionedVector PMLMin;
dimensionedScalar Rmin;
dimensionedScalar Rmax;
dimensionedVector origin;

switch (selectedModel)
{
    case RECTANGLE:
    {
        transportProperties.lookup("PMLMax") >> PMLMax;
        transportProperties.lookup("PMLMin") >> PMLMin;
        break;
    }
    case CYLINDRICAL:
    {
        transportProperties.lookup("Rmin") >> Rmin;
        transportProperties.lookup("Rmax") >> Rmax;
        transportProperties.lookup("origin") >> origin;
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
