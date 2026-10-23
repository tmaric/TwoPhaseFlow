/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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

#include "davofVolumeRule.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(davofVolumeRule, 0);
    defineRunTimeSelectionTable(davofVolumeRule, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofVolumeRule::davofVolumeRule
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    tolRel_(dict.getOrDefault<scalar>("curvedTolRel", 1e-11)),
    maxDepth_(dict.getOrDefault<label>("curvedMaxDepth", 12)),
    nLeaves_(0),
    nFallback_(0)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::davofVolumeRule> Foam::davofVolumeRule::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    const word ruleType
    (
        dict.getOrDefault<word>("volumeRule", "splitQuadrature")
    );

    Info<< "isoAdvection: volumeRule " << ruleType << endl;

    auto* ctorPtr = dictionaryConstructorTable(ruleType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "volumeRule",
            ruleType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<davofVolumeRule>(ctorPtr(mesh, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::davofVolumeRule::splitCellLiquidVolume
(
    const label celli,
    const davof::cellSurface& S,
    const scalar tolRel
) const
{
    return davof::cellLiquidVolumeQuadric
    (
        mesh_, celli, [&](const point& x) { return S.value(x); },
        tolRel, maxDepth_, nLeaves_
    );
}


Foam::scalar Foam::davofVolumeRule::splitPrismLiquidVolume
(
    const label facei,
    const vector& d,
    const davof::cellSurface& S,
    const scalar tolRel
) const
{
    return davof::faceSweptLiquidVolume
    (
        mesh_, facei, d, [&](const point& x) { return S.value(x); },
        tolRel, maxDepth_, nLeaves_
    );
}


Foam::scalar Foam::davofVolumeRule::faceLiquidFraction
(
    const label facei,
    const vector& d,
    const davof::cellSurface& S,
    label& nFallback
) const
{
    return davof::faceCutDisplaced
    (
        mesh_, facei, d, [&](const point& x) { return S.value(x); }, nFallback
    );
}


// ************************************************************************* //
