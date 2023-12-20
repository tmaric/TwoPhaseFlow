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

#include "gravityReconImprove.H"
#include "addToRunTimeSelectionTable.H"
#include "gravityMeshObject.H"
#include "reconstructionSchemes.H"

#include "plane.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gravityReconImprove, 0);
    addToRunTimeSelectionTable(accelerationForceModel,gravityReconImprove, components);
}

Foam::vector Foam::gravityReconImprove::closestDist(const point p, const vector n ,const vector centre)
{
    if (mag(n) != 0)
    {
        point testIC = vector(0, 0, 0.5145);
        // normal = vector(0,0,-1);
        // return centre; //p - normal*((p - centre) & normal);
        return testIC;
    } else {
        point testIC = vector(0, 0, 0.5145); 
        return testIC;
        //return centre; 
        }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gravityReconImprove::gravityReconImprove
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

void Foam::gravityReconImprove::calculateAcc()
{
        const fvMesh& mesh = acc_.mesh();
        const uniformDimensionedVectorField& g = meshObjects::gravity::New(mesh.time());
        g_.value() = g.value();
        auto constantIC = dimensionedVector("interfaceHeight",dimLength,vector(0, 0, 0.5145));

        acc_ = g_ & constantIC;
        accf_ = g_ & constantIC;       
/*    const fvMesh& mesh = acc_.mesh();
    // zoneDistribute& exchangeFields = zoneDistribute::New(mesh);
    // g_.value() = gravityDict_.lookup("gravity");

    // assume except for interface cells' faces, g&x on all other faces are equalt to zero
    acc_ = g_ & mesh.C();
    accf_ = g_ & mesh.Cf();
    const uniformDimensionedVectorField& g = meshObjects::gravity::New(mesh.time());
    g_.value() = g.value();
    dimensionedScalar hRef("hRef",dimLength, gravityDict_.lookupOrDefault("hRef",0));


    dimensionedScalar ghRef
    (
        mag(g_.value()) > SMALL
      ? g_ & (cmptMag(g_.value())/mag(g_.value()))*hRef
      : dimensionedScalar("ghRef", g_.dimensions()*dimLength, 0)
    );


//    acc_ = (g_ & mesh.C()) - ghRef;
//    accf_ = ((g_ & mesh.Cf()) - ghRef);

    if(mesh.foundObject<reconstructionSchemes>("reconstructionScheme"))
    {
        reconstructionSchemes& surf = mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

        surf.reconstruct(false);

        const volVectorField& faceCentre = surf.centre();
        const volVectorField& faceNormal = surf.normal();

        // needed mesh information
        const auto& own = mesh.faceOwner();
        const auto& neigh = mesh.faceNeighbour();
        const surfaceVectorField& Cf = mesh.Cf();

        scalarField& accfIn = accf_.primitiveFieldRef();

        // update internal field 
        forAll(accfIn, fi)
        {

            label ownf = own[fi];
            label neighf = neigh[fi];

             if (surf.interfaceCell()[ownf] && !(surf.interfaceCell()[neighf]))
            {
                // vector closeP = faceCentre[ownf];
                vector closeP = closestDist( Cf[fi], faceNormal[ownf], faceCentre[ownf]);
                // Info << "## faceCentre of [" << ownf << "]: " << faceCentre[ownf] << "; closeP[" << fi << "]: " << closeP << endl; 
                accf_[fi] = closeP & g_.value();
            }
            else if (surf.interfaceCell()[neighf] && surf.interfaceCell()[ownf])
            {
                //vector closeP = 0.5*(faceCentre[ownf] + faceCentre[neighf]);
                vector closeP = 0.5* (closestDist( Cf[fi], faceNormal[ownf], faceCentre[ownf]) 
                                   + closestDist( Cf[fi], faceNormal[neighf], faceCentre[neighf]));
                // Info << "## faceCentre of [" << ownf << "]: " << faceCentre[ownf] << " and faceCentre of ["<< neighf << "]: " << faceCentre[neighf] << "; closeP[" << fi << "]: " << closeP << endl;
                accf_[fi] = closeP & g_.value();
            }
            else if (surf.interfaceCell()[neighf] && !(surf.interfaceCell()[ownf])) 
            {
                // vector closeP = faceCentre[neighf];
                vector closeP = closestDist( Cf[fi], faceNormal[neighf], faceCentre[neighf]);
                // Info << "## faceCentre of [" << neighf << "]: " << faceCentre[neighf] << "; closeP[" << fi << "]: " << closeP << endl;
                // Info << "## unite faceNormal of ["<< neighf << "]: " << faceNormal[neighf]/mag(faceNormal[neighf]) << endl;
                accf_[fi] = closeP & g_.value();
            }
            else
            { continue; }
        }

        surfaceScalarField::Boundary& accfBf = accf_.boundaryFieldRef();
        const surfaceVectorField::Boundary& CfBf = Cf.boundaryField();

        forAll(accfBf, pi)
        {
            fvsPatchScalarField& paccf = accfBf[pi];
            const fvsPatchVectorField& pCf = CfBf[pi];

            if(paccf.coupled())
            {
                forAll(paccf, i)
                {
                    //special handle for cyclic and processor BC
                    // for future
                }
            }
            else
            {
                forAll(paccf, i)
                {
                    const label celli = paccf.patch().faceCells()[i];

                    if(surf.interfaceCell()[celli])
                    {
                        // vector closeP = faceCentre[celli]; 
                        vector closeP = closestDist( pCf[i], faceNormal[celli], faceCentre[celli]);
                        // Info << "## faceCentre of [" << celli << "]: " << faceCentre[celli] << "; closeP[" << i + paccf.patch().start() << "]: " << closeP << endl;
                        paccf[i] = closeP & g_.value();
                    }
                }
            }
        }

    }
    else
    {
        Info << "notFound" << endl;
    }*/

}

Foam::tmp<Foam::surfaceScalarField> Foam::gravityReconImprove::accelerationForce()
{
    const fvMesh& mesh = acc_.mesh();
    const volScalarField& rho = mesh.lookupObject<volScalarField>("rho");
    return -accf()*fvc::snGrad(rho);
}




// ************************************************************************* //
