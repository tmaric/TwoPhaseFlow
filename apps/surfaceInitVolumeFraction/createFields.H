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

surfaceScalarField alphaSurface
(
    IOobject
    (
        surfaceFieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(alpha)
);

