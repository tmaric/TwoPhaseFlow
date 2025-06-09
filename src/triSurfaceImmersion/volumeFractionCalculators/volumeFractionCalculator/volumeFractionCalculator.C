/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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

\*---------------------------------------------------------------------------*/

#include "surfaceFieldsFwd.H"
#include "surfaceInterpolate.H"
#include "volumeFractionCalculator.H"

#include "DynamicList.H"
#include "IOobject.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "fvc.H"
#include "pointIndexHit.H"
#include "pointMesh.H"

namespace Foam::TriSurfaceImmersion
{
defineTypeNameAndDebug(volumeFractionCalculator, 0);
defineRunTimeSelectionTable(volumeFractionCalculator, Dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
volumeFractionCalculator::volumeFractionCalculator(
    const dictionary& configDict, const fvMesh& mesh)
    : mesh_{mesh}, runTime_{mesh.time()}, pMesh_{mesh},
      writeGeometry_(configDict.get<Switch>("writeGeometry"))
{
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
autoPtr<volumeFractionCalculator> volumeFractionCalculator::New(
    const dictionary& configDict, const fvMesh& mesh)
{
    const word name = configDict.get<word>("algorithm");
    auto* ctorPtr = DictionaryConstructorTable(name);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            configDict,
            "volumeFractionCalculator",
            name,
            *DictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<volumeFractionCalculator>(ctorPtr(configDict, mesh));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void volumeFractionCalculator::bulkVolumeFraction(volScalarField& alpha)
{
    alpha = pos(this->sigDistCalc().cellSignedDist());
}


void volumeFractionCalculator::bulkAreaFraction(surfaceScalarField& alpha)
{
    alpha = pos(fvc::interpolate(this->sigDistCalc().cellSignedDist()));    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// ************************************************************************* //
