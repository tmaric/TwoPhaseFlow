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

#include "gravityVDirac.H"
#include "addToRunTimeSelectionTable.H"
#include "gravityMeshObject.H"
#include "reconstructionSchemes.H"
#include "surfaceInterpolationScheme.H"
#include "ListOps.H"
#include "processorPolyPatch.H"

#include "plane.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gravityVDirac, 0);
    addToRunTimeSelectionTable(accelerationForceModel,gravityVDirac, components);
}
/*
Foam::vector Foam::gravityVDirac::closestDist(const point p, const vector n ,const vector centre)
{
    vector normal = n/mag(n);
    return p - normal*((p - centre) & normal);
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gravityVDirac::gravityVDirac
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
    mixture_(mesh.lookupObject<volVectorField>("U"), mesh.lookupObject<surfaceScalarField>("phi")),
   diracV_
    (
        IOobject
        (
            "diracV",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("diracV_", dimAcceleration, vector(0.0, 0.0, 0.0)),
        "zeroGradient"
    ),
    diracVf_
    (
        IOobject
        (
            "diracVf",
            mesh.time().timeName(),
            mesh
        ),
    mesh,
    dimensionedVector("diracVf_", dimAcceleration, vector(0.0, 0.0, 0.0))
    )
    //cutFace_(mesh)
{
    // calculateAcc();
}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //

void Foam::gravityVDirac::calculateAcc()
{
    const fvMesh& mesh = diracV_.mesh();
    const uniformDimensionedVectorField& g = meshObjects::gravity::New(mesh.time());
    g_.value() = g.value();
    dimensionedScalar hRef("hRef",dimLength, gravityDict_.lookupOrDefault("hRef",0));
    
    dimensionedScalar ghRef
    (
        mag(g_.value()) > SMALL
      ? g_ & (cmptMag(g_.value())/mag(g_.value()))*hRef
      : dimensionedScalar("ghRef", g_.dimensions()*dimLength, 0)
    );

    if(mesh.foundObject<reconstructionSchemes>("reconstructionScheme"))
    {
        reconstructionSchemes& surf = mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

        surf.reconstruct(false);

        const volVectorField& faceCentre = surf.centre();
        const volVectorField& faceNormal = surf.normal();
        
        const DynamicField<label>& interfaceLabels = surf.interfaceLabels();
        forAll(interfaceLabels, i)
        {
            label celli = interfaceLabels[i];

            if (mag(faceCentre[celli]) != 0.0)
            {
                diracV_[celli] =((g_ & faceCentre[celli]).value() - ghRef.value())*(faceNormal[celli]/mesh.V()[celli]); 
            }
        }

        diracV_.correctBoundaryConditions();
    }
    else
    {
        Info << "notFound" << endl;
    }
}

Foam::tmp<Foam::surfaceScalarField> Foam::gravityVDirac::accelerationForce()
{
    const dimensionedScalar& rhoJump=mag(mixture_.rho1() - mixture_.rho2());
    const fvMesh& mesh = diracV_.mesh();

    return -(fvc::flux(diracV_))/mesh.magSf()*rhoJump; //fvc::snGrad(rho);
}




// ************************************************************************* //
