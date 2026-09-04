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
    acousticWaveFoam

Group
    grpAcousticSolvers

Description
    Time-domain acoustic solver for the pressure wave equation with
    heterogeneous media, non-reflecting PML, and radiation-pressure outputs.

Author
    Jun Liu, MMA, TU Darmstadt, 30. Oct. 2025
    Email: liu@mma.tu-darmstadt.de
SourceFiles
    acousticWaveFoam.C

\*---------------------------------------------------------------------------*/

#include <petscksp.h>
#include "fvCFD.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "isoAdvection.H"
#include "cutFaceAdvect.H"
#include "surfaceIteratorPLIC.H"
#include "reconstructionSchemes.H"
#include "upwind.H"
#include "processorPolyPatch.H"
#include "processorLduInterface.H"
#include "processorBC.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"

#include "petscTimeBlockAssembly.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Acoustic solver solving the acoustic pressure wave equation."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
   // #include "readTransportProperties.H"
    #include "createFields.H"
    #include "createPMLFields.H"


    Info<< "\nStarting time loop\n" << endl;
    #include "computeAlphaf.H"
    #include "computePMLCoefs.H"

    globalIndex globalCells(mesh.nCells());
    Mat timeBlockMatrix = nullptr;
    Vec timeBlockSolution = nullptr;
    Vec timeBlockRhs = nullptr;
    KSP timeBlockKsp = nullptr;
    bool timeBlockMatrixAssembled = false;
    scalar timeBlockMatrixDeltaT = -GREAT;

    if (blockCoupledTimeIntegration)
    {
        PetscInitialize(&argc, &argv, nullptr, nullptr);
        initializeTimeBlockPetscSystem
        (
            mesh,
            mesh.nCells(),
            globalCells.size(),
            timeBlockMatrix,
            timeBlockSolution,
            timeBlockRhs,
            timeBlockKsp,
            timeBlockLinearSolver
        );
    }

    #include "createRadiationAveraging.H"

    rho = alpha1*rhol + (1 - alpha1)*rhog;
    compressibility = alpha1*kl + (1 - alpha1)*kg;
    invRhof = 1/(alphaf*rhol + (1 - alphaf)*rhog);

    while (runTime.run())
    {
        tmp<volScalarField> previousPressureLaplacian;
        tmp<volVectorField> previousPressureGradient;

        if
        (
            blockCoupledTimeIntegration
         || timeIntegrationMethod == "thetaNewmark"
        )
        {
            // Function-based patch coefficients are not retained by oldTime().
            // Capture the old spatial operators before advancing the clock.
            p.correctBoundaryConditions();
            previousPressureGradient = fvc::grad(p);

            if (blockCoupledTimeIntegration)
            {
                previousPressureLaplacian =
                    -rho*fvc::laplacian
                    (
                        invRhof,
                        p,
                        "laplacian(invRhof,p)"
                    );
            }
            else
            {
                previousPressureLaplacian =
                    -rho*fvc::laplacian
                    (
                        invRhof*physicalTimeWeightf,
                        p,
                        "laplacian(invRhof,p)"
                    );
            }
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (blockCoupledTimeIntegration)
        {
            #include "timeBlockEqn.H"
        }
        else
        {
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
        }

        if (blockCoupledTimeIntegration)
        {
            solve
            (
                fvm::ddt(rho, U)
             ==
                - timeIntegrationTheta*fvc::grad(p)
                - (1.0 - timeIntegrationTheta)*previousPressureGradient()
            );
        }
        else if (timeIntegrationMethod == "thetaNewmark")
        {
            pDot ==
                (p - p.oldTime())
               /(timeIntegrationTheta*runTime.deltaT())
              - ((1.0 - timeIntegrationTheta)/timeIntegrationTheta)
               *pDot.oldTime();

            solve
            (
                fvm::ddt(rho, U)
             ==
                - timeIntegrationTheta*fvc::grad(p)
                - (1.0 - timeIntegrationTheta)*previousPressureGradient()
            );
        }
        else
        {
            solve
            (
                fvm::ddt(rho, U) == -fvc::grad(p)
            );
        }
        prT ==
            0.5*compressibility*sqr(p)
           - 0.5*rho*magSqr(U);
        momFluxT == rho*(U*U);

        #include "updateRadiationFields.H"
        runTime.write();

        runTime.printExecutionTime(Info);
    }

    if (blockCoupledTimeIntegration)
    {
        KSPDestroy(&timeBlockKsp);
        VecDestroy(&timeBlockSolution);
        VecDestroy(&timeBlockRhs);
        MatDestroy(&timeBlockMatrix);
        PetscFinalize();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
