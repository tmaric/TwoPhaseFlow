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

#include "gravityDirac.H"
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
    defineTypeNameAndDebug(gravityDirac, 0);
    addToRunTimeSelectionTable(accelerationForceModel,gravityDirac, components);
}
/*
Foam::vector Foam::gravityDirac::closestDist(const point p, const vector n ,const vector centre)
{
    vector normal = n/mag(n);
    return p - normal*((p - centre) & normal);
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gravityDirac::gravityDirac
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
//    snGradAlpha_(fvc::snGrad(mixture_.alpha1())), //"snGradAlpha",dimless/dimLength, fvc::snGrad(mixture_.alpha1())),
    mixture_(mesh.lookupObject<volVectorField>("U"), mesh.lookupObject<surfaceScalarField>("phi")),
    snGradAlpha_(fvc::snGrad(mixture_.alpha1())),
   /* snGradAlpha_(
        IOobject
        (
            "snGradAlpha",
            mesh.time().timeName(),
            mesh
        ),
    mesh,
    dimensionedScalar("snGradAlpha_", dimless/dimLength, 0.0)
    ),*/
    cutFace_(mesh)
{
    // calculateAcc();
}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::gravityDirac::calculateDiracInteg (
            const label faceID,
            scalar accf,
            const scalar snGradAlpha,
            const scalar rmagSf,
            const vector Sf,
            const vector faceNormal,
            const vector faceCentre,
            const vector cellFaceInterfaceNormal
        )
{

    vector fN = faceNormal;
    vector fIN = cellFaceInterfaceNormal;
    vector nf = Sf*rmagSf;

    if((mag(faceNormal) != 0) && (mag(cellFaceInterfaceNormal) != 0))
    {
        fN /= mag(fN);
        fIN /= mag(fIN);
// only for test, delete latter
        fN = vector(0,0,-1);
        fIN = vector(0,0,-1);
//    Info << "### print interfaceNormal: " << fN 
//         << " print cellFaceInterfaceNormal: " << fIN << endl;
        label cutStatus = cutFace_.calcSubFace(faceID, fN, faceCentre);
        
        if (cutStatus == 0)
        {
            if (cutFace_.surfacePoints().size() == 2)
                {
                    // ### METHOD 1:  calculate from cell face geometric integration of Delta sampling
                    
                    point Lc = 0.5*(cutFace_.surfacePoints()[0] + cutFace_.surfacePoints()[1]);
                    dimensionedScalar Llen = mag(cutFace_.surfacePoints()[0] - cutFace_.surfacePoints()[1]);
                    //Info << "print cutted line center z: " << Lc.z() 
                    // << " print cutted line length: " << Llen.value() << endl;
                    accf = (rmagSf*Llen*(g_ & Lc)*(fIN & nf)).value();
                    

                    // ### METHOD 2: calculate from snGradAlpha
                    // accf = ((g_ & faceCentre)*snGradAlpha).value();
                    
                    // ### METHOD 3: direction correction --> doesnot work
                    // accf = ((g_ & faceCentre)*snGradAlpha*(fIN & nf)).value();
                    // snGradAlpha = 1;
                }
        }
        else {
                // ### METHOD 2 face Centre replaced by closest point on interface
                 accf = ((g_ & faceCentre)*snGradAlpha).value();
                // ### METHOD 3
                //    accf = ((g_ & faceCentre)*snGradAlpha*(fIN & nf)).value(); 
                // snGradAlpha = 1; 
            }
    }
    return accf;
}


void Foam::gravityDirac::calculateAcc()
{
    const fvMesh& mesh = acc_.mesh();
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
    snGradAlpha_ == fvc::snGrad(mesh.lookupObject<volScalarField>("alpha.phase1"));


    if(mesh.foundObject<reconstructionSchemes>("reconstructionScheme"))
    {
        reconstructionSchemes& surf = mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

        surf.reconstruct(false);

        const volVectorField& faceCentre = surf.centre();
        const volVectorField& faceNormal = surf.normal();

        // JUN: interpolate cell interface normal to cell face, correct the boundaryField/patch faces value.
        surfaceVectorField cellFaceNormal(fvc::interpolate(faceNormal));
        tmp<surfaceInterpolationScheme<vector>> interp(surfaceInterpolationScheme<vector>::New
        (
            mesh,
            mesh.interpolationScheme("interpolate(" + cellFaceNormal.name() + ")")
        ));

        surfaceVectorField::Boundary& cellFaceNormalbf = cellFaceNormal.boundaryFieldRef(); 
        forAll(cellFaceNormalbf, pi)
        {
            surfaceScalarField weights(interp.ref().weights(faceNormal));
            fvsPatchScalarField pWeights = weights.boundaryField()[pi];
            if (cellFaceNormal.boundaryField()[pi].coupled())
            {
                cellFaceNormalbf[pi] = pWeights*faceNormal.boundaryField()[pi].patchInternalField()
                  + (1.0 - pWeights)*faceNormal.boundaryField()[pi].patchNeighbourField();

            }
            else
            {
                cellFaceNormalbf[pi] = faceNormal.boundaryField()[pi].patchInternalField();

            }
        }

        // needed mesh information
        const auto& own = mesh.faceOwner();
        const auto& neigh = mesh.faceNeighbour();
        const surfaceScalarField& rmagSf(1/mesh.magSf());
        surfaceScalarField::Boundary rmagSfb(mesh.magSf().boundaryField());
        rmagSfb = 1/mesh.magSf().boundaryField();
        
        const surfaceVectorField& Sf(mesh.Sf());
        scalarField& accfIn = accf_.primitiveFieldRef();
        snGradAlpha_.setOriented();

        // update internal field 
        forAll(accfIn, fi)
        {

            label ownf = own[fi];
            label neighf = neigh[fi];

             if (surf.interfaceCell()[ownf] && !(surf.interfaceCell()[neighf]))
            {
                //Info << "### print faceid " << fi << "; snGradAlpha_: " << snGradAlpha_[fi] << "; accf_ " << accf_[fi] << endl;
                accf_[fi] = calculateDiracInteg(fi, accf_[fi], snGradAlpha_[fi], rmagSf[fi], Sf[fi], faceNormal[ownf], faceCentre[ownf], cellFaceNormal[fi]);
                snGradAlpha_[fi] = 1;
                //Info << "### print faceid " << fi << "; snGradAlpha_: " << snGradAlpha_[fi] << "; accf_ " << accf_[fi] << endl;
            }
            else if (surf.interfaceCell()[neighf] && surf.interfaceCell()[ownf])
            {
                //Info << "### print faceid " << fi << "; snGradAlpha_: " << snGradAlpha_[fi] << "; accf_ " << accf_[fi] << endl;
                accf_[fi] = 0.5*(calculateDiracInteg(fi, accf_[fi], snGradAlpha_[fi], rmagSf[fi], Sf[fi], faceNormal[ownf], faceCentre[ownf], cellFaceNormal[fi]) + 
                        calculateDiracInteg(fi, accf_[fi], snGradAlpha_[fi], rmagSf[fi], Sf[fi], faceNormal[neighf], faceCentre[neighf], cellFaceNormal[fi]));
                snGradAlpha_[fi] = 1;
                //Info << "### print faceid " << fi << "; snGradAlpha_: " << snGradAlpha_[fi] << "; accf_ " << accf_[fi] << endl;
            }
            else if (surf.interfaceCell()[neighf] && !(surf.interfaceCell()[ownf])) 
            {
                //Info << "### print faceid " << fi << "; snGradAlpha_: " << snGradAlpha_[fi] << "; accf_ " << accf_[fi] << endl;
                accf_[fi] = calculateDiracInteg(fi, accf_[fi], snGradAlpha_[fi], rmagSf[fi], Sf[fi], faceNormal[neighf], faceCentre[neighf], cellFaceNormal[fi]);
                snGradAlpha_[fi] = 1;
                //Info << "### print faceid " << fi << "; snGradAlpha_: " << snGradAlpha_[fi] << "; accf_ " << accf_[fi] << endl;
            }
            else
            { continue; }
        }

        
        surfaceScalarField::Boundary& accfBf = accf_.boundaryFieldRef();
        surfaceScalarField::Boundary& snGradAlphaBf = snGradAlpha_.boundaryFieldRef();

        forAll(accfBf, pi)
        {
            fvsPatchScalarField& paccf = accfBf[pi];
            fvsPatchScalarField& psnGradAlpha = snGradAlphaBf[pi];

            if(paccf.coupled())
            {
                forAll(paccf, i)
                {
                    //const label celli = paccf.patch().faceCells()[i];
                    // jun: get/send information, still need expolare 
                }
            }
            else
            {
                forAll(paccf, i)
                {
                    const label celli = paccf.patch().faceCells()[i];

                    if(surf.interfaceCell()[celli])
                    {
                        label fi = i + paccf.patch().start();
                        // Info << "### print faceid " << fi << "; snGradAlpha_: " << psnGradAlpha[i] << "; accf_ " << paccf[i] << endl;
                        paccf[i] = calculateDiracInteg(
                                                        fi, 
                                                        paccf[i], 
                                                        psnGradAlpha[i], 
                                                        rmagSfb[pi][i], 
                                                        Sf.boundaryField()[pi][i], 
                                                        faceNormal[celli], 
                                                        faceCentre[celli], 
                                                        cellFaceNormal.boundaryField()[pi][i]
                                                       );
                        psnGradAlpha[i] = 1;
                        // Info << "### print faceid " << fi << "; snGradAlpha_: " << psnGradAlpha[i] << "; accf_ " << paccf[i] << endl;
                    }
                    //Info << "### print patchid: " << pi << ", accf " << paccf[i] << endl;
                }
            }
        }
        //forAll(snGradAlphaBf[0], fi)
        //{ Info << "### print snGradAlphaBf on patch 0 : " << snGradAlphaBf[0][fi] << endl;}

    }
    else
    {
        Info << "notFound" << endl;
    }
}

Foam::tmp<Foam::surfaceScalarField> Foam::gravityDirac::accelerationForce()
{
    const dimensionedScalar& rhoJump = mixture_.rho1() - mixture_.rho2();

    //snGradAlpha_.setOriented();
    return -accf()*rhoJump*snGradAlpha_; //fvc::snGrad(rho);
}




// ************************************************************************* //
