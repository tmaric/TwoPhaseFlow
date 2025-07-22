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

    // Compute bulk area fractions for boundary faces
    auto& aboundary = alpha.boundaryFieldRef();

    forAll(aboundary, patchI)
    {
        auto& patchField = aboundary[patchI];
        const auto& patch = patchField.patch();
        const auto& cf = patch.Cf();

        forAll(patchField, faceI)
        {
            patchField[faceI] = pos(this->sigDistCalc().signedDistance(cf[faceI]));
        }
    }
}


void volumeFractionCalculator::bulkAreaFraction(scalarField& alpha)
{
    const auto& cf = mesh_.faceCentres();

    auto bulkFractionsTmp = pos(fvc::interpolate(this->sigDistCalc().cellSignedDist()));
    const auto& bulkFractions = bulkFractionsTmp.cref();
    
    // Transfer bulk fractions of internal faces
    forAll(bulkFractions, faceI)
    {
        alpha[faceI] = bulkFractions[faceI];
    }

    // Compute bulk fractions of boundary faces
    for (label faceI = mesh_.nInternalFaces(); faceI != alpha.size(); ++faceI)
    {
        alpha[faceI] = pos(this->sigDistCalc().signedDistance(cf[faceI]));
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// ************************************************************************* //
