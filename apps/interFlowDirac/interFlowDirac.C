#include "advectionSchemes.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "surfaceForces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Reconstructs a PLIC interface and initializes surface fields associated"
	" with PLIC centroids."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"
    #include "createDeltaFields.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    volScalarField RDF 
    (
        IOobject
        (
            "reconstructedDistanceFunction",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField cSigma
    (
	IOobject
	(
	    "cSigma", 
	    runTime.timeName(), 
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::NO_WRITE
	),
	mesh
    );


    // Compute interfaceDirac from reconstructed distance 
    #include "computeDelta.H"

    // Compute interfacial excess concentration 
    volScalarField cDelta ("cDelta", cSigma * delta);

    delta.write();
    cDelta.write();

    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
