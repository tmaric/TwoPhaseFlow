/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
   interMPIALFlow

Description
    Two-phase flow solver coupled with acoustic wave fields.

Author
    Jun Liu, MMA, TU Darmstadt
    Email: liu@mma.tu-darmstadt.de
SourceFiles
    interMPIALFlow.C
\*---------------------------------------------------------------------------*/

#include "advectionSchemes.H"
#include <petscksp.h>
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
#include "cutFaceAdvect.H"
#include "surfaceIteratorPLIC.H"
#include "reconstructionSchemes.H"
#include "upwind.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "processorPolyPatch.H"
#include "processorBC.H"
#include "syncTools.H"
#include "acousticBlockPetsc.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using isoAdvector phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing.\n"
        "Acoustic-flow coupling skeleton for interMPIALFlow."
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
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    const bool overwrite = args.found("overwrite");

    turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    PetscInitialize(&argc, &argv, nullptr, nullptr);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        if(overwrite)
        {
            runTime.setTime(runTime.value() - runTime.deltaTValue(), 1);
            runTime.writeAndEnd();
        }

        Info<< "Time = " << runTime.timeName() << nl << endl;
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                advector->surf().reconstruct();

                mesh.update();

                if (mesh.changing())
                {
                    // gets recompute by surfaces forces
                    // gh = (g & mesh.C()) - ghRef;
                    // ghf = (g & mesh.Cf()) - ghRef;
                    advector->surf().mapAlphaField();
                    alpha2 = 1.0 - alpha1;
                    alpha2.correctBoundaryConditions();
                    rho == alpha1*rho1 + alpha2*rho2;
                    rho.correctBoundaryConditions();
                    rho.oldTime() = rho;
                    alpha2.oldTime() = alpha2;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            if(overwrite)
            {
                continue;
            }

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            surfForces.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "acousticCoupling.H"
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

//            #include "computeTimeAverage.H"

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    PetscFinalize();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
