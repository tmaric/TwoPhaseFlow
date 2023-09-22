// Volume fraction field
#include "dimensionedVector.H"
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

volVectorField interfaceNormal 
(
    IOobject
    (
        IOobject::groupName(normalFieldName, alpha.group()),
        runTime.timeName(), 
        mesh, 
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimArea, Zero)
);

