/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
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

#include "davofFaceClipping.H"
#include "davof.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
    defineTypeNameAndDebug(davofFaceClipping, 0);
    defineRunTimeSelectionTable(davofFaceClipping, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofFaceClipping::davofFaceClipping
(
    const word& type,
    const dictionary& dict,
    const reconstruction::davof& scheme
)
:
    scheme_(scheme),
    mesh_(scheme.mesh()),
    nCorrected_(0),
    nCapped_(0),
    nClamped_(0),
    nSkippedCos_(0),
    nSkippedChord_(0)
{}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::davofFaceClipping> Foam::davofFaceClipping::New
(
    const dictionary& dict,
    const reconstruction::davof& scheme
)
{
    const word modelType(dict.getOrDefault<word>("faceClipping", "plane"));

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "faceClipping",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<davofFaceClipping>(ctorPtr(dict, scheme));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::davofFaceClipping::resetCounters()
{
    nCorrected_ = 0;
    nCapped_ = 0;
    nClamped_ = 0;
    nSkippedCos_ = 0;
    nSkippedChord_ = 0;
}


void Foam::davofFaceClipping::report(Ostream& os) const
{
    os  << "sagitta corrected "
        << returnReduce(nCorrected_, sumOp<label>()) << " faces, "
        << returnReduce(nCapped_, sumOp<label>()) << " capped, "
        << returnReduce(nClamped_, sumOp<label>()) << " clamped, "
        << returnReduce(nSkippedCos_, sumOp<label>()) << " skipped (cos), "
        << returnReduce(nSkippedChord_, sumOp<label>()) << " skipped (chord)";
}


// ************************************************************************* //
