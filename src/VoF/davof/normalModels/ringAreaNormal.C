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

#include "ringAreaNormal.H"
#include "davof.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofNormalModels
{
    defineTypeNameAndDebug(ringArea, 0);
    addToRunTimeSelectionTable(davofNormalModel, ringArea, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofNormalModels::ringArea::ringArea
(
    const dictionary& dict,
    const reconstruction::davof& scheme
)
:
    davofNormalModel(typeName, dict, scheme, "all"),
    mInterface_
    (
        IOobject
        (
            IOobject::groupName("davofRingArea", scheme.alpha1().group()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedVector(dimArea, Zero)
    ),
    filled_()
{}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::davofNormalModels::ringArea::fillInterfaceNormals()
{
    davofNormalModel::fillInterfaceNormals(scheme_, mInterface_, filled_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::davofNormalModels::ringArea::modelNormals
(
    const scalar,
    vectorField& n
)
{
    // The one implementation shared with the noNormalFallback ringArea
    // rule of the base class
    ringAreaNormals(scheme_, mInterface_, filled_, n);
}


// ************************************************************************* //
