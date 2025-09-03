// Calculate minimal discretization length.
// Assumes interface is resolved properly.
const scalar hMin = gMin(pow(mesh.deltaCoeffs(), -1)());

// Calculate the interface dirac
const scalar sigma = 4*hMin;      // 1 stdev within 2*h_min
const scalar m = 3.0;           // truncate at 3 sigma
const scalar support = m * sigma;

volScalarField delta
(
    IOobject
    (
	"delta",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("delta", dimless/dimLength, 0.0)
);

