dimensionedScalar rhog
(
    "rhog",
    dimDensity,
    transportProperties
);

dimensionedScalar rhol
(
    "rhol",
    dimDensity,
    transportProperties
);

dimensionedScalar kg
(
    "kg",
    dimCompressibility/dimDensity,
    transportProperties
);

dimensionedScalar kl
(
    "kl",
    dimCompressibility/dimDensity,
    transportProperties
);

dimensionedScalar cg
(
    "cg",
    dimVelocity,
    transportProperties
);

dimensionedScalar f
(
    "f",
    dimless/dimTime,
    transportProperties
);

Info<< "Reading PML setup parameters\n" << endl;

dimensionedScalar L
(
    "L",
    dimless,
    transportProperties
);

dimensionedScalar sigmaMax
(
    "sigmaMax",
    dimless/dimTime,
    transportProperties.lookupOrDefault<scalar>("sigmaMax", 0.0)
);

dimensionedScalar po
(
    "po",
    dimless,
    transportProperties.lookupOrDefault<scalar>
    (
        "po",
        transportProperties.lookupOrDefault<scalar>("p", 2.0)
    )
);

enum PMLType
{
    RECTANGLE,
    SPHERICAL,
    UNDEFINED
};

const word PMLTypeWord = transportProperties.get<word>("PMLType");
PMLType selectedModel = UNDEFINED;

if (PMLTypeWord == "rectangle")
{
    selectedModel = RECTANGLE;
}
else if (PMLTypeWord == "spherical")
{
    selectedModel = SPHERICAL;
}
else
{
    FatalErrorInFunction
        << "Unknown PMLType '" << PMLTypeWord
        << "' in transportProperties. Valid types: rectangle, spherical."
        << exit(FatalError);
}

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
    case SPHERICAL:
    {
        transportProperties.lookup("Rmin") >> Rmin;
        transportProperties.lookup("Rmax") >> Rmax;
        transportProperties.lookup("origin") >> origin;
        break;
    }
    default:
    {
        FatalErrorInFunction
            << "Internal error: undefined PML type."
            << exit(FatalError);
        break;
    }
}
