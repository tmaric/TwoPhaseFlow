/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 Tomislav Maric, TU Darmstadt
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
    testReconstructionSchemes

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "reconstructionSchemes.H"
#include "areaFractionReconstruction.H"

using namespace Foam::reconstruction;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    const fvSolution& fvSolutionDict(mesh); 
    const dictionary& solverDict = fvSolutionDict.subDict("solvers");
    const dictionary& alphaDict = solverDict.subDict("alpha.water*"); 

    // Reconstruct the interface using a scheme from fvSolution.solver.alpha.water* 
    autoPtr<reconstructionSchemes> recon = 
        reconstructionSchemes::New(alpha1, phi, U, alphaDict);
    recon->reconstruct();
    
    // Compute area fractions using available reconstruction data. 
    areaFractionReconstruction areaFractionRecon(alpha1, phi, U, alphaDict);
    
    areaFractionRecon.calcAreaFractions
    (
        alpha1, 
        recon->centre(), 
        recon->normal(), 
        recon->interfaceCell()
    );

    const surfaceScalarField& alphaf = areaFractionRecon.areaFractions();
    alphaf.write();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
