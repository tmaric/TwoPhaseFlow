
// Reconstruct interface at time 0
advector->surf().reconstruct();
const volScalarField& RDF = mesh.lookupObject<volScalarField>("reconstructedDistanceFunction");
RDF.write();

#include "computeDelta.H"
delta.write();

volScalarField cSigma
(
    IOobject
    (
	"cSigma",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField cDelta
(
    IOobject
    (
	"cDelta",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    cSigma * delta, 
    "zeroGradient"
);
cDelta.write();
