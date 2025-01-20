/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 AUTHOR,AFFILIATION
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

#include "vofForcesFunctionObject.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "fvcReconstruct.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(vofForcesFunctionObject, 0);
    addToRunTimeSelectionTable(functionObject, vofForcesFunctionObject, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::vofForcesFunctionObject::vofForcesFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    rho_(mesh_.lookupObject<volScalarField>("rho")),
    U_(mesh_.lookupObject<volVectorField>("U")),
    p_rgh_(mesh_.lookupObject<volScalarField>("p_rgh")),
    phi_(mesh_.lookupObject<surfaceScalarField>("phi")),
    rhoPhi_(mesh_.lookupObject<surfaceScalarField>("rhoPhi")),
    mixture_(U_, phi_),
    surfForces_(mixture_.alpha1(), phi_, U_, mixture_),
    ddtRhoU_
    (
        IOobject
        (
            "ddtRhoU", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("ddtRhoU", dimForce / dimVolume, pTraits<vector>::zero)
    ),
    inertialForce_ 
    (
        IOobject
        (
            "inertialForce", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("divRhoPhiU", dimForce / dimVolume, pTraits<vector>::zero)
    ),
    viscousForce_ 
    (
        IOobject
        (
            "viscousForce", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("divDevRhoReff", dimForce / dimVolume, pTraits<vector>::zero)
    ),
    surfaceTensionForce_ 
    (
        IOobject
        (
            "surfaceTensionForce", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("surfaceTensionForce", dimForce / dimVolume, pTraits<vector>::zero)
    ),
    accelerationForce_ 
    (
        IOobject
        (
            "accelerationForce", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("accelerationForce", dimForce / dimVolume, pTraits<vector>::zero)
    ),
    pressureForce_ 
    (
        IOobject
        (
            "pressureForce", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("pressureForce", dimForce / dimVolume, pTraits<vector>::zero)
    )
{
    read(dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::vofForcesFunctionObject::read(const dictionary& dict)
{
    return true;
}

bool Foam::functionObjects::vofForcesFunctionObject::execute()
{

    ddtRhoU_ = fvc::ddt(rho_, U_);

    inertialForce_ = fvc::div(rhoPhi_, U_);

    volTensorField gradU = fvc::grad(U_, "pointCellsLeastSquares");

    viscousForce_ = fvc::div(mixture_.mu() * (gradU + gradU.T()));

    surfForces_.correct();
    surfaceTensionForce_ = fvc::reconstruct(surfForces_.surfaceTensionForce() * 
                                            mesh_.magSf());

    accelerationForce_ = fvc::reconstruct(surfForces_.accelerationForce() * mesh_.magSf());

    pressureForce_ = fvc::reconstruct(-fvc::snGrad(p_rgh_) * mesh_.magSf());

    return true;
}


bool Foam::functionObjects::vofForcesFunctionObject::end()
{
    return true;
}


bool Foam::functionObjects::vofForcesFunctionObject::write()
{
    return true;
}


// ************************************************************************* //
