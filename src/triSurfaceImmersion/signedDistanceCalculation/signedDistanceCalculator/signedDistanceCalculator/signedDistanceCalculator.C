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

#include "signedDistanceCalculator.H"

#include "dictionary.H"

namespace Foam::TriSurfaceImmersion
{

defineTypeNameAndDebug(signedDistanceCalculator, 0);
defineRunTimeSelectionTable(signedDistanceCalculator, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * *

signedDistanceCalculator::signedDistanceCalculator(
    const dictionary& configDict, const fvMesh& mesh)
    : dict_{configDict}, mesh_{mesh}, pMesh_{mesh},
      narrowBandWidth_{configDict.get<scalar>("narrowBandWidth")},
      outOfNarrowBandValue_{configDict.get<scalar>("bulkValue")},
      cellSignedDist_{IOobject("cellSignedDist",
                          mesh.time().timeName(),
                          mesh,
                          IOobject::NO_READ,
                          IOobject::NO_WRITE),
          mesh,
          dimensionedScalar("cellSignedDist", dimLength, 0),
          "zeroGradient"},
      cellSignedDist0_{"cellSignedDist0", cellSignedDist_},
      cellNearestTriangle_{}, pointSignedDist_{IOobject("pointSignedDist",
                                                   mesh.time().timeName(),
                                                   mesh_,
                                                   IOobject::NO_READ,
                                                   IOobject::NO_WRITE),
                                  pMesh_,
                                  dimensionedScalar(
                                      "pointSignedDist", dimLength, 0),
                                  "zeroGradient"},
      pointNearestTriangle_{}
{
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<signedDistanceCalculator> signedDistanceCalculator::New(
    const dictionary& configDict, const fvMesh& mesh)
{
    const word name = configDict.get<word>("surfaceType");

    auto* ctorPtr = DictionaryConstructorTable(name);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            configDict,
            "signedDistanceCalculator",
            name,
            *DictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<signedDistanceCalculator>(ctorPtr(configDict,mesh));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void signedDistanceCalculator::outOfNarrowBandValue(const scalar value)
{
    outOfNarrowBandValue_ = value;
}


void signedDistanceCalculator::narrowBandWidth(const scalar width)
{
    narrowBandWidth_ = width;
}


void signedDistanceCalculator::writeFields() const
{
    cellSignedDist0_.write();
    cellSignedDist_.write();
    pointSignedDist_.write();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
