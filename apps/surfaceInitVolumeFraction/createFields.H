// Volume fraction field
volScalarField alpha
(
    IOobject
    (
        fieldName, 
        runTime.timeName(), 
        mesh, 
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

