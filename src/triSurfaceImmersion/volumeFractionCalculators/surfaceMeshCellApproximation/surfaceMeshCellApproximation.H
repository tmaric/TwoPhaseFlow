/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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
    Foam::TriSurfaceImmersion::surfaceMeshCellApproximation

Description
    Implements the surface-mesh-cell-approximation algorithm (SMCA) using
    tetrahedral refinement and tet volume fraction approximation
    to compute volume fractions of interface cells.

SourceFiles
    surfaceMeshCellApproximation.cpp
    surfaceMeshCellApproximationI.H

\*---------------------------------------------------------------------------*/

#ifndef surfaceMeshCellApproximation_H
#define surfaceMeshCellApproximation_H

#include <vector>
#include <volFieldsFwd.H>

#include "AdaptiveTetCellRefinement.H"
#include "AdaptiveTriFaceRefinement.H"
#include "signedDistanceCalculator.H"
#include "surfaceFieldsFwd.H"
#include "volumeFractionCalculator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion
{

/*---------------------------------------------------------------------------*\
               Class surfaceMeshCellApproximation Declaration
\*---------------------------------------------------------------------------*/

class surfaceMeshCellApproximation : public volumeFractionCalculator
{
private:
    
    // Private Types
    //- A triple containing set of indexed tetrahedra, points and sig.distances
    using cellDecompositionTuple = std::
        tuple<std::vector<indexedTet>, std::vector<point>, std::vector<scalar>>;
    //- A triple containing set of indexed triangles, points and signed distances
    using faceDecompositionTuple = std::
        tuple<std::vector<indexedTri>, std::vector<point>, std::vector<scalar>>;

    // Private Data
    //- Pointer to the signed distance calculator holding interface information
    autoPtr<signedDistanceCalculator> sigDistCalcPtr_;

    //- Set of interface cell IDs
    std::vector<label> interfaceCellIDs_;

    //- Relative global volume error threshold to determine refinement level
    scalar relVolumeErrorThreshold_;

    //- Maximum refinement level allowed for tetrahedral refinement
    label maxAllowedRefinementLevel_;

    //- Maximum refinement level actually used for tetrahedral refinement
    label maxUsedRefinementLevel_ = 0;

    // Private Member Functions
    //- Decompose the cell into tetrahedra using its centre and face centres
    cellDecompositionTuple decomposeCell(label cellID) const;

    //- Decompose the face into triangles at its centroid
    faceDecompositionTuple decomposeFace(label faceID) const;

    //- Number of tetrahedra a cell is decomposed into
    label nTets(label cellID) const;

    //- Flag faces of interface cells by setting a negative area fraction value
    void flagInterfaceCellFaces(surfaceScalarField& alpha) const;

    //- Compute volume fractions for interface cells.
    //  Returns the maximum refinement level used.
    label interfaceCellVolumeFraction(volScalarField& alpha,
            bool writeDecomposition);

    //- Compute area fractions for faces of interface cells.
    //  Returns the maximum refinement level used.
    label intersectedFacesAreaFraction(surfaceScalarField& alpha,
            bool writeDecomposition);

    //- Compute the target refinement level
    void determineRefinementLevel(volScalarField& alpha);


public:

    // Static Data Members
    TypeName("SMCA");


    // Constructors
    surfaceMeshCellApproximation(
        const dictionary& configDict, const fvMesh& mesh);


    // Member Functions
    //- Return average number of triangles per cell
    inline scalar nTrianglesPerCell() const override;

    //- Return number of cells which may be intersected by the interface
    inline label nIntersectedCells() const override;

    //- Return the maximum refinement level used for tetrahedral refinement
    inline label maxRefinementLevel() const override;

    //- Return reference to the signed distance calculator
    const signedDistanceCalculator& sigDistCalc() const override;

    //- Compute the volume fraction field alpha
    void calcVolumeFraction(volScalarField& alpha) override;

    //- Compute the area fraction field alpha
    void calcAreaFraction(surfaceScalarField& alpha) override;

    //- Find all cells which may be intersected by the interface
    void findIntersectedCells() override;

    //- Write additional fields used during the computation of volume fractions
    void writeFields() const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceMeshCellApproximationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
