/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "gravityConstantHeight.H"
#include "addToRunTimeSelectionTable.H"
#include "gravityMeshObject.H"
#include "reconstructionSchemes.H"

#include "plane.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gravityConstantHeight, 0);
    addToRunTimeSelectionTable(accelerationForceModel,gravityConstantHeight, components);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gravityConstantHeight::gravityConstantHeight
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    accelerationForceModel
    (
        typeName,
        dict,
        mesh
    ),
    gravityDict_(dict),
    g_
    (
        "gravity",
        dimAcceleration,
        vector(0,0,0)
    ),
    origin_
    (
        dict.get<vector>("origin")
    )
{
    // calculateAcc();
}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //

void Foam::gravityConstantHeight::calculateAcc()
{
        const fvMesh& mesh = acc_.mesh();
        const uniformDimensionedVectorField& g = meshObjects::gravity::New(mesh.time());
        g_.value() = g.value();
        auto constantIC = dimensionedVector("interfaceHeight",dimLength,origin_);

        acc_ = g_ & constantIC;
        accf_ = g_ & constantIC;
}

Foam::tmp<Foam::surfaceScalarField> Foam::gravityConstantHeight::accelerationForce()
{
    const fvMesh& mesh = acc_.mesh();
    const volScalarField& rho = mesh.lookupObject<volScalarField>("rho");
    return -accf()*fvc::snGrad(rho);
}




// ************************************************************************* //
