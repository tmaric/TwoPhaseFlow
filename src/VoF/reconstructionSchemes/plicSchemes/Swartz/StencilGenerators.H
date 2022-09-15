/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022 TU Darmstadt 
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

Class
    Foam::WideCellStencilGenerator

Description
    Cell stencil generator using cell-corner-point neighbor connectivity. 

Authors
    Tomislav Maric maric@mma.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#ifndef StencilGenerators_H
#define StencilGenerators_H

#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

/*---------------------------------------------------------------------------*\
                         WideCellStencilGenerator 
\*---------------------------------------------------------------------------*/

template<typename Stencil = std::set<label>>
struct CellCornerCellStencil 
{
    public:

        Stencil calcStencil(const label& cellLabel, const fvMesh& mesh) const
        {
            Stencil result; 

            const cellList& cells = mesh.cells(); 
            const faceList& faces = mesh.faces();

            // Get the labels of points that belong to the cell with cellLabel.
            const labelList cellPointLabels = cells[cellLabel].labels(faces); 
            const labelListList& pointCells = mesh.pointCells(); 
            
            forAll(cellPointLabels, pointI)
            {
                // Get the labels of all cells that interect a point of the 
                // cell with cellLabel.
                const labelList& pointCellLabels = 
                    pointCells[cellPointLabels[pointI]];

                // Add all the point-cell labels to the set of unique labels. 
                forAll(pointCellLabels, cellI)
                    result.insert(pointCellLabels[cellI]);
            }

            return result;
        }
}; 

template<typename Stencil = std::set<label>>
struct FaceCornerFaceStencil 
{
    public:

        Stencil calcStencil(const label& faceLabel, const primitivePatch& primPatch) const
        {
            Stencil result; 
            const auto& facePointFaces = primPatch.pointFaces()[faceLabel];

            for(const auto& faceL : facePointFaces)
                result.insert(faceL);

            return result;
        }
}; 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
