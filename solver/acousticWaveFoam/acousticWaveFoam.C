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

#include "fvCFD.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "isoAdvection.H"
#include "cutFaceAdvect.H"
#include "surfaceIteratorPLIC.H"
#include "reconstructionSchemes.H"
#include "upwind.H"
#include "processorPolyPatch.H"
#include "processorBC.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"

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
    scalar sampleCount = 0;
    
    rho = 1/(alpha1/rhol + (1 - alpha1)/rhog);
    invRhof = alphaf/rhol + (1 - alphaf)/rhog;

    while (runTime.run())
    {
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (pimple.correct())
        {
            #include "pEqn.H"
        }

        if (laplacianCorrection > 0.5)
        {
            pDot ==
                (p - p.oldTime())
               /(timeIntegrationTheta*runTime.deltaT())
              - ((1.0 - timeIntegrationTheta)/timeIntegrationTheta)
               *pDot.oldTime();

            solve
            (
                fvm::ddt(rho, U)
             == -fvc::grad
                 (
                     timeIntegrationTheta*p
                   + (1.0 - timeIntegrationTheta)*p.oldTime()
                 )
            );
        }
        else
        {
            solve
            (
                fvm::ddt(rho, U) == -fvc::grad(p)
            );
        }

        sampleCount += 1;
        prT == 0.5*(kl*alpha1 + kg*(1-alpha1))*p*p - 0.5*rho*(U&U);
        momFluxT == rho*(U*U);

        const scalar oldWeight = (sampleCount - 1)/sampleCount;
        const scalar newWeight = 1/sampleCount;
        pr == oldWeight*pr + newWeight*prT;
        momFlux == oldWeight*momFlux + newWeight*momFluxT;

        if
        (
            mesh.foundObject<volScalarField>("prTMean")
         && mesh.foundObject<volTensorField>("momFluxTMean")
        )
        {
            pr == mesh.lookupObject<volScalarField>("prTMean");
            momFlux == mesh.lookupObject<volTensorField>("momFluxTMean");
        }
        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
