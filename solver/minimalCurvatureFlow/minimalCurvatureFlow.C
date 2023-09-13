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
        "Solver for minimal curvature flow "
        " using isoAdvector phase-fraction-based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing.\n"
        "The solver is derived from interFlow"
    );

    Foam::argList::addBoolOption
    (
        "overwrite",
        "Update and overwrite the existing mesh useful for adaptive mesh refinement"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"

    const bool overwrite = args.found("overwrite");

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        ++runTime;

        if(overwrite)
        {
            runTime.setTime(runTime.value() - runTime.deltaTValue(), 1);
            runTime.writeAndEnd();
        }

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Store old flux.
        phi.oldTime(); 

        while (pimple.loop())
        {
            // Update surface tension force.
            surfTensionForce->correct();

            //  - Compute the mean curvature. 
            Ktilde = surfTensionForce->K(); 
            //  - Smooth the mean curvature. 
            for(label Ki= 0; Ki < 3; ++Ki)
            {
                Ktilde = fvc::average(Ktilde);
            }

            // - Compute the initial minimal curvature flux using lambdaf = 1 m^3 / s. 
            phi == -lambdaf * fvc::snGrad(Ktilde) * mesh.magSf();  

            // - Compute the CFL number
            alphaCo = mag(phi) * runTime.deltaT() * mesh.deltaCoeffs() / mesh.magSf(); 
            scalar gMaxAlphaCo = gMax(alphaCo); 
            Info << "Max alphaCo = " << gMaxAlphaCo << endl;
            Info << "Min alphaCo = " << gMin(alphaCo) << endl;

            // Scale lambdaf based on maximal Courant number
            scalar maxAlphaCo = runTime.controlDict().get<scalar>("maxAlphaCo");
            lambdaf *= (maxAlphaCo / gMaxAlphaCo);

            // Recompute the flux with the scaled lambdaf that ensure the Courant condition is upheld.
            phi == -lambdaf * fvc::snGrad(Ktilde) * mesh.magSf();  

            // - Compute the CFL number
            alphaCo = mag(phi) * runTime.deltaT() * mesh.deltaCoeffs() / mesh.magSf(); 
            gMaxAlphaCo = gMax(alphaCo); 
            Info << "Max alphaCo from scaled lambdaf = " << gMaxAlphaCo << endl;
            Info << "Min alphaCo from scaled lambdaf = " << gMin(alphaCo) << endl;

            // TODO(TM): apply boundary conditions on the scaled flux using U boundary conditions.

            // - Perform Helmholz decomposition to ensure \sum_f \phi_f = 0.
            fvScalarMatrix psiEqn
            (
                fvm::laplacian(psi) == fvc::div(phi)
            );
            psiEqn.solve(); 
            phi = phi - psiEqn.flux();

            // Report discrete volume conservation 
            Info << "Helmholz decomposition, gMax(fvc::div(phi)) = " << gMax(divError) << endl; 
            Info << "Helmholz decomposition, gAverage(fvc::div(phi)) = " << gAverage(divError) << endl; 
            divError = fvc::div(phi);

            // - Compute the CFL number
            alphaCo = mag(phi) * runTime.deltaT() * mesh.deltaCoeffs() / mesh.magSf(); 
            gMaxAlphaCo = gMax(alphaCo); 
            Info << "Max alphaCo from divfree phi = " << gMaxAlphaCo << endl;
            Info << "Min alphaCo from divfree phi = " << gMin(alphaCo) << endl;

            // Solve for volume fractions using CFL-abiding divergence-free minimal-curvature flux.
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"
        }

        // Reconstruct U for visualization.
        U = fvc::reconstruct(phi);

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
