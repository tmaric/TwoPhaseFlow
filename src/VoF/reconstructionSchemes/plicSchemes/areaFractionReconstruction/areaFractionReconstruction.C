/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Tomislav Maric, TU Darmstadt 
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

#include "areaFractionReconstruction.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "linear.H"
#include "cutFacePLIC.H"
#include "processorPolyPatch.H"
#include "volFieldsFwd.H"
#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(areaFractionReconstruction, 0);
    addToRunTimeSelectionTable(reconstructionSchemes, areaFractionReconstruction, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::areaFractionReconstruction::areaFractionReconstruction
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    gradAlpha 
    (
        alpha1,
        phi,
        U,
        dict
    ), 
    mesh_(alpha1.mesh()),
    runTime_(mesh_.time()),
    surfCellTol_(modelDict().lookupOrDefault<scalar>("surfCellTol", 1e-8)),
    alphaf_ 
    (
        IOobject
        (
            IOobject::groupName("alphaf", alpha1.group()),
            runTime_.timeName(), 
            mesh_, 
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_, 
        dimensionedScalar(dimless, Zero)
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reconstruction::areaFractionReconstruction::calcAreaFractions
(
    const volScalarField& alpha1, 
    const volVectorField& centres,
    const volVectorField& normals, 
    const boolList& interfaceCell 
)
{
    // Compute initial area fractions by linear interpolation to correctly  
    // set the values in the bulk of each phase.
    alphaf_ = linear<scalar>(mesh_).interpolate(alpha1); 
    
    const auto& own = mesh_.owner();
    const auto& nei = mesh_.neighbour();
    const auto& magSf = mesh_.magSf();

    cutFacePLIC cutter(mesh_);

    // Calculate face-area fractions of internal cell-faces. 
    forAll(own, faceI)
    {
        // Initialize face-owner and face-neighbor area fractions to 0. 
        scalar alphafOwn = 0; 
        scalar alphafNei = 0;
        
        // If the owner-cell contains the interface.  
        const auto ownCellI = own[faceI];
        if (interfaceCell[ownCellI])
        {
            // Calculate area fraction of the owner interface.
            cutter.calcSubFace(faceI, normals[ownCellI], centres[ownCellI]);
            alphafOwn = mag(cutter.subFaceArea()) / magSf[faceI];
        }
            
        // If the neighbor-cell contains the interface. 
        const auto neiCellI = nei[faceI];
        if (interfaceCell[neiCellI])
        {
            // Calculate the area fraction of the neighbor interface.
            cutter.calcSubFace(faceI, normals[neiCellI], centres[neiCellI]);
            alphafNei = mag(cutter.subFaceArea()) / magSf[faceI];
        }

        // If both areas are wetted
        if ((alphafOwn > surfCellTol_) && (alphafNei > surfCellTol_))
        {
            // Final area fraction is the average owner/neighbor fraction. 
            alphaf_[faceI] = 0.5*(alphafOwn + alphafNei);
        }       
        // Else if interface from the face-owner cuts the face
        else if ((alphafOwn > surfCellTol_) && ((alphafNei < surfCellTol_))) 
        {
            alphaf_[faceI] = alphafOwn; 
        }
        // Else if interface from the face-neighbor cuts the face
        else if ((alphafNei > surfCellTol_) && (alphafOwn < surfCellTol_)) 
        {
            alphaf_[faceI] = alphafNei; 
        }
        else if ((alphafNei < surfCellTol_) && (alphafOwn < surfCellTol_)) // No geometrical intersection 
        {
            // If the face belongs to a full cell. 
            if ((alpha1[ownCellI] > (1 - surfCellTol_)) || 
                (alpha1[neiCellI] > (1 - surfCellTol_))) 
            {  
                alphaf_[faceI] = 1;
            }
            else // the face is not wetted at all. 
            {
                alphaf_[faceI] = 0;
            }
        }
    }

    // For all alpha1_ boundary patches 
    const auto& alpha1BoundaryField = alpha1.boundaryField(); 
    const auto& normalBoundaryField = normals.boundaryField(); 
    const auto& centreBoundaryField = centres.boundaryField(); 
    const auto& magSfBoundaryField = magSf.boundaryField();
    const auto& faceOwner = mesh_.faceOwner();
    forAll(alpha1BoundaryField, patchI)
    {
        // For all faces in the boundary patch
        const auto& alpha1PatchField = alpha1BoundaryField[patchI];
        const auto& normalPatchField = normalBoundaryField[patchI];
        const auto& centrePatchField = centreBoundaryField[patchI];
        const auto& magSfPatchField  = magSfBoundaryField[patchI];
        
        // Get the patch-neighbor alpha, PLIC centres, 
        // and PLIC normals for coupled patches. 
        tmp<scalarField> alpha1PatchNeiFieldTmp;
        tmp<vectorField> centrePatchNeiFieldTmp;
        tmp<vectorField> normalPatchNeiFieldTmp;
        
        const auto& meshPatch = mesh_.boundary()[patchI]; 
        // TODO(TM): implement cyclic boundaries. 
        const bool isProcessorPatch = isA<processorPolyPatch>(meshPatch);
        if (isProcessorPatch)
        {
            alpha1PatchNeiFieldTmp = alpha1PatchField.patchNeighbourField();
            centrePatchNeiFieldTmp  = centrePatchField.patchNeighbourField();
            normalPatchNeiFieldTmp = normalPatchField.patchNeighbourField();
        }

        // Compute area fractions for a boundary patch
        forAll(alpha1PatchField, faceI)
        {
            // Initialize face-owner and face-neighbor area fractions to 0. 
            scalar alphafOwn = 0; 
            scalar alphafNei = 0;
            
            // If the face-owner cell contains the interface.  
            label faceG = meshPatch.start() + faceI;
            const auto ownCellI = faceOwner[faceG];
            if (interfaceCell[ownCellI])
            {
                // Calculate area fraction of the owner interface.
                cutter.calcSubFace
                (
                    faceG, 
                    normalPatchField[faceI], 
                    centrePatchField[faceI]
                );
                alphafOwn = mag(cutter.subFaceArea()) / magSfPatchField[faceI];
            }
                
            // If the patch is coupled, each face-owner cell has a patch-neighbor. 
            if (isProcessorPatch)
            {
                const scalarField& alpha1PatchNeiField = alpha1PatchNeiFieldTmp.cref(); 
                const vectorField& centrePatchNeiField = centrePatchNeiFieldTmp.cref();
                const vectorField& normalPatchNeiField = normalPatchNeiFieldTmp.cref();

                bool faceNeiCellHasInterface = 
                     (surfCellTol_ < alpha1PatchNeiField[faceI]) &&
                     (alpha1PatchNeiField[faceI] < (1 - surfCellTol_));
                
                if (faceNeiCellHasInterface)
                {
                    // Calculate the area fraction of the neighbor interface.
                    cutter.calcSubFace
                    (
                        faceG, 
                        normalPatchNeiField[faceI],
                        centrePatchNeiField[faceI] 
                    );
                    alphafNei = mag(cutter.subFaceArea()) / magSfPatchField[faceI];
                }
            }

            // If both areas are wetted
            if ((alphafOwn > surfCellTol_) && (alphafNei > surfCellTol_))
            {
                // Final area fraction is the average owner/neighbor fraction. 
                alphaf_[faceI] = 0.5*(alphafOwn + alphafNei);
            }  
            else if ((alphafOwn > surfCellTol_) && (alphafNei < surfCellTol_)) 
            {
                alphaf_[faceI] = alphafOwn; 
            }
            // Else if interface from the face-neighbor cuts the face
            else if ((alphafNei > surfCellTol_) && (alphafOwn < surfCellTol_)) 
            {
                alphaf_[faceI] = alphafNei; 
            }
            else if ((alphafNei < surfCellTol_) && (alphafOwn < surfCellTol_) && 
                     isProcessorPatch) // No geometrical intersection across processor patches.
            {
                const scalarField& alpha1PatchNeiField = alpha1PatchNeiFieldTmp.cref(); 
                // If the face belongs to a full cell. 
                if ((alpha1PatchField[faceI] > (1 - surfCellTol_)) || 
                    (alpha1PatchNeiField[faceI] > (1 - surfCellTol_))) 
                {  
                    alphaf_[faceI] = 1;
                }
                else // the face is not wetted at all. 
                {
                    alphaf_[faceI] = 0;
                }
            }
        }
    }
}

void Foam::reconstruction::areaFractionReconstruction::calcAreaFractions()
{
    calcAreaFractions(alpha1_, centre_, normal_, interfaceCell_);
}

void Foam::reconstruction::areaFractionReconstruction::calcAreaFractions
(
    const volScalarField& alpha1,
    const reconstructionSchemes& recon
)
{
    calcAreaFractions
    (
        alpha1, 
        recon.centre(), 
        recon.normal(), 
        recon.interfaceCell()
    );
}

void Foam::reconstruction::areaFractionReconstruction::reconstruct(bool forceUpdate)
{
    // Compute normals using Youngs' algorithm.
    gradAlpha::reconstruct(forceUpdate); 
    // Calculate area fractions using internal Youngs' data. 
    calcAreaFractions();
    // TODO(TM): compute the interface area-normal vectors from Youngs' data.
    -fvc::surfaceSum(alphaf_ * mesh_.Sf());
    // TODO(TM): position the interface.
}

// ************************************************************************* //
