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

#include "gravityRecon.H"
#include "addToRunTimeSelectionTable.H"
#include "gravityMeshObject.H"
#include "reconstructionSchemes.H"

#include "plane.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gravityRecon, 0);
    addToRunTimeSelectionTable(accelerationForceModel,gravityRecon, components);
}

Foam::vector Foam::gravityRecon::closestDist(const point p, const vector n ,const vector centre)
{
    vector normal = n/mag(n);
    // normal = vector(0,0,-1);
    // vector closeP = p - normal*((p - centre) & normal);

    // return centre;
    return p - normal*((p - centre) & normal);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gravityRecon::gravityRecon
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
    )


{
    // calculateAcc();
}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //

void Foam::gravityRecon::calculateAcc()
{
    const fvMesh& mesh = acc_.mesh();
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh);
    // g_.value() = gravityDict_.lookup("gravity");
    const uniformDimensionedVectorField& g = meshObjects::gravity::New(mesh.time());
    g_.value() = g.value();
    dimensionedScalar hRef("hRef",dimLength, gravityDict_.lookupOrDefault("hRef",0));


    dimensionedScalar ghRef
    (
        mag(g_.value()) > SMALL
      ? g_ & (cmptMag(g_.value())/mag(g_.value()))*hRef
      : dimensionedScalar("ghRef", g_.dimensions()*dimLength, 0)
    );


    acc_ = (g_ & mesh.C()) - ghRef;
    accf_ = ((g_ & mesh.Cf()) - ghRef);

    if(mesh.foundObject<reconstructionSchemes>("reconstructionScheme"))
    {
        reconstructionSchemes& surf = mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

        surf.reconstruct(false);

        const volVectorField& faceCentre = surf.centre();
        const volVectorField& faceNormal = surf.normal();

        boolList nextToInterface(mesh.nCells(),false);
        markInterfaceRegion markIF(mesh);

        markIF.markCellsNearSurf(surf.interfaceCell(),1,nextToInterface);

        exchangeFields.setUpCommforZone(nextToInterface,true);

        Map<vector > mapCentres =
            exchangeFields.getDatafromOtherProc(nextToInterface,faceCentre);
        Map<vector > mapNormal =
            exchangeFields.getDatafromOtherProc(nextToInterface,faceNormal);

        const labelListList& stencil = exchangeFields.getStencil();

        forAll(surf.interfaceCell(),celli)
        {
            if(mag(faceNormal[celli]) != 0)
            {
                vector closeP = closestDist(mesh.C()[celli],-faceNormal[celli],faceCentre[celli]);
                acc_[celli] = closeP & g_.value();
            }
            else if(nextToInterface[celli])
            {
                // the to the interface
                vector averageSurfP = vector::zero;
                scalar avgWeight = 0;
                const point p = mesh.C()[celli];

                forAll(stencil[celli],i)
                {
                    const label& gblIdx = stencil[celli][i];
                    vector n = -exchangeFields.getValue(faceNormal,mapNormal,gblIdx);
                    if (mag(n) != 0)
                    {
                        n /= mag(n);
                        vector c =
                            exchangeFields.getValue(faceCentre,mapCentres,gblIdx);
                        vector distanceToIntSeg = (c - p);
                        vector closeP = closestDist(p,n,c);
                        scalar weight = 0;

                        if (mag(distanceToIntSeg) != 0)
                        {
                            distanceToIntSeg /= mag(distanceToIntSeg);
                            weight = sqr(mag(distanceToIntSeg & n));
                        }
                        else // exactly on the center
                        {
                            weight = 1;
                        }
                        averageSurfP += closeP * weight;
                        avgWeight += weight;
                    }
                }

                if (avgWeight != 0)
                {
                    averageSurfP /= avgWeight;
                    acc_[celli] = averageSurfP & g_.value();
                }

            }
            else
            {
                // do nothing
            }

        }

        acc_.correctBoundaryConditions();
        accf_ = fvc::interpolate(acc_);
        Info << "### find min/max of accf_: " << gMin(accf_) << "/" << gMax(accf_) << endl;

    }
    else
    {
        Info << "notFound" << endl;
    }

}

Foam::tmp<Foam::surfaceScalarField> Foam::gravityRecon::accelerationForce()
{
    const fvMesh& mesh = acc_.mesh();
    const volScalarField& rho = mesh.lookupObject<volScalarField>("rho");
    return -accf()*fvc::snGrad(rho);
}




// ************************************************************************* //
