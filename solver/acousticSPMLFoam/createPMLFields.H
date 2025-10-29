Info<< "Reading PML related fields \n" << endl;

volVectorField u
(
    IOobject
    (
        "u",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("u", dimPressure/dimLength, vector::zero)
);

    volScalarField v
    (
        IOobject
        (
            "v",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("v", dimPressure*dimTime, Zero)
    );

// volVectorField sgm
// (
//     IOobject
//     (
//         "sgm",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedVector("sgm", dimless/dimTime, vector::zero)
// );

// alpha coeficient in eq.15 
    volScalarField SC0
    (
        IOobject
        (
            "SC0",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("SC0", dimTime/(dimLength*dimLength), Zero)
    );
// beta coeficient in eq.15
    volScalarField SC1
    (
        IOobject
        (
            "SC1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("SC1", dimless/(dimLength*dimLength), Zero)
    );
// gamma coeficient in eq.15
    volScalarField SC2
    (
        IOobject
        (
            "SC2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("SC2", dimless/(dimLength*dimLength*dimTime), Zero)
    );
// tensor A in eq.16
volTensorField TC0
(
    IOobject
    (
        "TC0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("TC0", dimless/dimTime, tensor::zero)
);
//tensor B in eq.17
volTensorField TC1
(
    IOobject
    (
        "TC1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("TC1", dimless/dimTime, tensor::zero)
);
// tensor C in eq.16
volTensorField TC2
(
    IOobject
    (
        "TC2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("TC2", dimless/(dimTime*dimTime), tensor::zero)
);
// frequency
    //  dimensionedScalar f
    //  (
    //  //    thermodynamicProperties.lookup("psi")
    //   "f",
    //   dimless/dimTime,
    //   transportProperties
    //  );

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
/*else if (waveTypeWord == "spherical")
{
    selectedModel = SPHERICAL;
}*/
else if (PMLTypeWord == "cylindrical")
{
    selectedModel = CYLINDRICAL;
}
else
{
    // Error handling is now centralized and clear
    FatalErrorInFunction
        << "Unknown PMLType '" << PMLTypeWord
        << "' in transportProperties." << nl
        << "Valid types are: rectangle, spherical, cylindrical."
        << exit(FatalError);
}

//PML max x, y ,z
dimensionedVector PMLMax;
/*    (
      "PMLMax", dimless, transportProperties
    );*/
//PML min x, y, z
    dimensionedVector PMLMin;
/*    (
      "PMLMin", dimless, transportProperties
    );*/

dimensionedScalar Rmin;
/*(
      //    thermodynamicProperties.lookup("psi")
        "Rmin",
       dimLength,
       transportProperties
);*/
dimensionedScalar Rmax;
/*(
      //    thermodynamicProperties.lookup("psi")
        "Rmax",
       dimLength,
       transportProperties
);*/
    dimensionedVector origin;
/*    (
      "origin", dimless, transportProperties
    );*/


switch (selectedModel)
{
    case RECTANGLE:
    {
        //Info<< "  Wave type: linear" << endl;
        transportProperties.lookup("PMLMax") >> PMLMax;
        transportProperties.lookup("PMLMin") >> PMLMin;
        
       /* Info<< "  Read amplitude: " << amplitude << nl
            << "  Read period: " << period << endl;*/
        break; // <-- Don't forget the 'break'!
    }

    case CYLINDRICAL:
    {
        //Info<< "  Wave type: stokesFifth" << endl;
        transportProperties.lookup("Rmin") >> Rmin;
        transportProperties.lookup("Rmax") >> Rmax;
        transportProperties.lookup("origin") >> origin;

       /* Info<< "  Read amplitude: " << amplitude << nl
            << "  Read waveLength: " << waveLength << nl
            << "  Read stokesOrder: " << stokesOrder << endl;*/
        break;
    }
    case UNDEFINED:
    {
        // This case should not be reachable because of the 
        // FatalError, but it's good practice for completeness.
        FatalErrorInFunction
            << "Internal error: PML type is undefined."
            << exit(FatalError);
        break;
    }
}
