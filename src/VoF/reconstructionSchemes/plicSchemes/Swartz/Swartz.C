/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 TU Darmstadt 
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

Class
    Foam::reconstruction::Swartz

Description
    Reconstructs an interface (centre and normal vector) consisting of planes
    to match the internal fluid distribution in cells. 

    Reference:
    \verbatim
        Marić, T., Marschall, H., & Bothe, D. (2018). 
        An enhanced un-split face-vertex flux-based VoF method. 
        Journal of Computational Physics, 371, 967–993.
        https://doi.org/10.1016/j.jcp.2018.03.048
    \endverbatim

Authors
    Tomislav Maric, TU Darmstadt (2022)

SourceFiles
    Swartz.H
    Swartz.C

\*---------------------------------------------------------------------------*/

#include "Swartz.H"
#include "StencilGenerators.H" 
#include "fvcGrad.H"
                               
#include "processorFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(Swartz, 0);
    addToRunTimeSelectionTable(reconstructionSchemes,Swartz, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::Swartz::Swartz
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    reconstructionSchemes
    (
        typeName,
        alpha1,
        phi,
        U,
        dict
    ),
    iteration_(modelDict().lookupOrDefault("iterations" , 3)),
    isoFaceTol_(modelDict().lookupOrDefault("isoFaceTol" , 1e-08)),
    critAngle_(modelDict().lookupOrDefault("criticalAngle" , 30)),
    youngsNormals_
    (
        IOobject
        (
            IOobject::groupName("youngsNormals", alpha1.group()),
            alpha1_.mesh().time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(dimless, Zero)
    ),
    swartzNormals_
    (
        IOobject
        (
            IOobject::groupName("swartzNormals", alpha1.group()),
            alpha1_.mesh().time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(dimless, Zero)
    ),
    procPatchIds_(),
    sIterPLIC_(youngsNormals_.mesh(), isoFaceTol_)
{
    // Get process boundary patch IDs
    const fvMesh& mesh = U.mesh();
    const auto& meshBoundary = mesh.boundary(); 
    forAll(meshBoundary, patchI)
    {
        const fvPatch& patch = meshBoundary[patchI]; 
        // FIXME: should be coupled, not only processor. TM.
        if (isA<processorFvPatch>(patch))
        {
            procPatchIds_.append(patchI); 
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reconstruction::Swartz::improveNormals()
{
    CellCornerCellStencil<std::set<label>> cellCornerStencil; 

    // Compute the Swartz normals in serial. 
    const fvMesh& mesh = swartzNormals_.mesh();
    forAll(interfaceLabels_, ifaceC)
    {
        const label cellC = interfaceLabels_[ifaceC];

        const vector& xc = centre_[cellC];
        const vector& ncYoungs = youngsNormals_[cellC];

        auto stencil = cellCornerStencil.calcStencil(cellC, mesh); 
            
        // For all interface-cell corner neighbors. 
        for(const auto cellN : stencil)
        { 
            // If the neighbor cell is an interface cell.
            if (interfaceCell_[cellN] && (cellN != cellC))
            {
                // If the angle between the normals is smaller than critical. 
                // TODO: Add the critical angle as dictionary parameter. TM.
                if ((ncYoungs & youngsNormals_[cellN])
                     < cos(critAngle_ * M_PI / 180.))
                {
                    // Get the interface plane position.  
                    const vector& xn = centre_[cellN];
                    
                    // Centroid->neihbor centroid vector. 
                    vector xcxn = xn - xc;
                    
                    // Rotate the xcxn vector 90 degrees w.r.t the rotation 
                    // plane given by the existing interface normal.
                    quaternion qrot
                    (
                        (xcxn) ^ ncYoungs, 
                        constant::mathematical::pi * 0.5
                    ); 

                    // Add the rotated xcxn to the local Swartz normal.
                    swartzNormals_[cellC] += qrot.transform(xcxn);
                }
            }
        } // End for all interface-cell corner neighbors 
    } // End for all interface cells

    // Add contributions to the Swartz normals across MPI boundaries. 
    const auto& meshBoundary = mesh.boundary(); 
    const auto& faceOwner = mesh.faceOwner(); 
    const auto& centreBdryField = centre_.boundaryField();
    const auto& alphaBdryField = alpha1_.boundaryField();
    forAll(procPatchIds_, patchI)
    {
        // MPI process patch label.  
        const auto procPatchI = procPatchIds_[patchI]; 

        // Get the patch-neighbor PLIC centres.
        const auto& centrePatchField = centreBdryField[procPatchI]; 
        auto centrePatchNeiFieldTmp = 
            centrePatchField.patchNeighbourField();
        const auto& centrePatchNeiField = centrePatchNeiFieldTmp();

        // Get the patch-neighbor alpha values.
        const auto& alphaPatchField = alphaBdryField[procPatchI]; 
        auto alphaPatchNeiFieldTmp = 
            alphaPatchField.patchNeighbourField();
        const auto& alphaPatchNeiField = alphaPatchNeiFieldTmp();

        // For all faces in of the MPI process boundary patch.
        const fvPatch& patch = meshBoundary[procPatchI]; 
        const auto& pPatch = patch.patch(); 
        forAll(pPatch, faceI)
        {
            const label cellC = faceOwner[patch.start() + faceI];
            // If the owner (local) cell contains an interface element. 
            if (interfaceCell_[cellC])
            {
                // Get the owner-cell interface polygon centroid.
                const auto& xc = centre_[cellC];

                // Get the owner-cell interface Youngs normal.
                const auto& ncYoungs = youngsNormals_[cellC];

                // Get the faceCornerFace stencil.
                FaceCornerFaceStencil<std::set<label>> faceStencil; 
                const auto stencil = faceStencil.calcStencil(faceI, pPatch);
                
                // For all faces in the faceCornerFace stencil. 
                for(const auto& faceJ : stencil)
                {
                    // If the face-neighbor is an interface cell 
                    if ((alphaPatchNeiField[faceJ] > isoFaceTol_) && 
                        (alphaPatchNeiField[faceJ] < 1 - isoFaceTol_))
                    {
                        // Get the face-neighbor centroid.
                        const auto& xn = centrePatchNeiField[faceJ];  
                        
                        // Centroid->neihbor centroid vector. 
                        vector xcxn = xn - xc;
                        
                        // Rotate the xcxn vector 90 degrees w.r.t the rotation 
                        // plane given by the existing interface normal.
                        quaternion qrot
                        (
                            (xcxn) ^ ncYoungs, 
                            constant::mathematical::pi * 0.5
                        ); 

                        // Add the rotated xcxn to the Swartz normal.
                        swartzNormals_[cellC] += qrot.transform(xcxn);
                    }
                } // End loop over faceCornerFace stencil elements
            } // End if face-owner is an interface cell
        } // End for all faces in MPI process patch
    } // End for all MPI process patches 

    // Normalize the Swartz normal field. 
    forAll(interfaceLabels_, ifaceC)
    {
        label cellC = interfaceLabels_[ifaceC];
        swartzNormals_[cellC] /= mag(swartzNormals_[cellC]);
    }

}

void Foam::reconstruction::Swartz::positionInterface(const volVectorField& normals)
{
    forAll(interfaceCell_, ifaceC)
    {
        label cellC = interfaceCell_[ifaceC];

        sIterPLIC_.vofCutCell
        (
            cellC,
            alpha1_[cellC],
            isoFaceTol_,
            100,
            normals[cellC]
        );

        if (sIterPLIC_.cellStatus() == 0)
        {
            normal_[cellC] = sIterPLIC_.surfaceArea();
            centre_[cellC] = sIterPLIC_.surfaceCentre();
            if (mag(normals[cellC]) == 0)
            {
                normal_[cellC] = Zero;
                centre_[cellC] = Zero;
            }
        }
        else
        {
            normal_[cellC] = Zero;
            centre_[cellC] = Zero;
        }
    }
}

void Foam::reconstruction::Swartz::reconstruct(bool forceUpdate)
{
    // Estimate Youngs' normal using LSQ grad(alpha)
    auto alphaGradTmp = fvc::grad(alpha1_, "pointCellsLeastSquares"); 
    auto alphaMagGradTmp = mag(alphaGradTmp()); 
    alphaMagGradTmp.ref() += 
        dimensionedScalar("SMALL", Foam::pow(dimLength, -1), SMALL);
    youngsNormals_ = alphaGradTmp() / alphaMagGradTmp();

    positionInterface(youngsNormals_);

    // Initialize Swartz normals to Youngs normals
    swartzNormals_ = youngsNormals_;
    
    for(label it = 0; it < iteration_; ++it)
    {
        // Use Swartz algorithm to improve the normal.
        improveNormals();
        // Use Swartz normals to position the interface. 
        // Set normal_ to area-normal vectors given by the PLIC polygons.
        // Set center_ to the centroids of PLIC polygons.
        positionInterface(swartzNormals_);
    }
}

// TODO: enable, shared with plicRDF. TM. 
void Foam::reconstruction::Swartz::mapAlphaField() const
{
    //addProfilingInFunction(geometricVoF);
    // without it we seem to get a race condition
    //mesh_.C();

    //cutCellPLIC cutCell(mesh_);

    //forAll(normal_, celli)
    //{
        //if (mag(normal_[celli]) != 0)
        //{
            //vector n = normal_[celli]/mag(normal_[celli]);
            //scalar cutValue = (centre_[celli] - mesh_.C()[celli]) & (n);
            //cutCell.calcSubCell
            //(
                //celli,
                //cutValue,
                //n
            //);
            //alpha1_[celli] = cutCell.VolumeOfFluid();

        //}
    //}
    //alpha1_.correctBoundaryConditions();
    //alpha1_.oldTime () = alpha1_;
    //alpha1_.oldTime().correctBoundaryConditions();
}


// ************************************************************************* //
