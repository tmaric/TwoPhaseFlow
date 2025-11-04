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
   acousticTDSPMLFoam

Group
    grpAcousticSolvers

Description
    Time domain acoustic solver solving the acoustic pressure wave equation and calculating acoustic radiation pressure, 
    considering inhomogenous scattering and non-reflective PML.

Author
    Jun Liu, MMA, TU Darmstadt
    Email: liu@mma.tu-darmstadt.de
SourceFiles
    acousticTDSPMLFoam.C

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
//    #include "readTransportProperties.H"
    #include "createFields.H"
    #include "createPMLFields.H"


    Info<< "\nStarting time loop\n" << endl;
    #include "computeAlphaf.H"
    #include "computePMLCoefs.H"
    float iter=0.0;
    
    while (runTime.run())
    {
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (pimple.correct())
        {
            #include "pEqn.H"
        }

        iter=iter+1.0;
//        AcousticForcePotential == (iter - 1)/iter * AcousticForcePotential + 1/iter * (2*Foam::constant::mathematical::pi*pow(R,3)*( 1/(3.0*rho*sqr(cg))*p*p - rho/2*(U&U)));
        prT == 0.5*(kl*alpha1 + kg*(1-alpha1))*p*p - 0.5*rho*(U&U);
        momFluxT == rho*(U*U);
        if(runTime.timeIndex() > runTime.startTimeIndex() + 1)
        {
            if(mesh.foundObject<volScalarField>("prTMean") && mesh.foundObject<volTensorField>("momFluxTMean"))
            {
                pr == mesh.lookupObject<volScalarField>("prTMean");
                momFlux == mesh.lookupObject<volTensorField>("momFluxTMean");
            }
            else
            {
                        FatalErrorInFunction
                            << "Did not find field prTMean and momFluxTMean, please check if add below functionObject to controlDict: "
                            << nl
                            << "functions { " << nl
                            << "ARFAverage" << nl
                            << "    type            fieldAverage;" << nl
                            << "    libs            (fieldFunctionObjects);" << nl
                            << "    writeControl    writeTime;" << nl
                            << "" << nl
                            << "    fields" << nl
                            << "    (" << nl
                            << "        prT" << nl
                            << "        {" << nl
                            << "            mean            on;" << nl
                            << "            prime2Mean      off;" << nl
                            << "            base            time;" << nl
                            << "            windowType      approximate; // approximation of averaging a window, can save a lot space for old time fields." << nl
                            << "        window        #eval{$N/$f}; // set the average of field in last N periods. f is the frequency." << nl 
                            << "            windowName      movingAverageprTWindow;" << nl
                            << "            allowRestart    no;" << nl
                            << "        }" << nl
                            << "" << nl
                            << "        momFluxT" << nl
                            << "        {" << nl
                            << "            mean            on;" << nl
                            << "            prime2Mean      off;" << nl
                            << "            base            time;" << nl
                            << "            windowType      approximate;" << nl
                            << "        window        #eval{$N/$f};" << nl
                            << "            windowName      movingAverageWindow;" << nl
                            << "            allowRestart    no;" << nl
                            << "        } ); }" << nl
                            << nl           
                        
                            << exit(FatalError);
            }
        }
        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
