/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 Tobias Tolle, TU Darmstadt
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
    Foam::TriSurfaceImmersion::triSurfaceDistCalc

Author
    Tobias Tolle
    tolle@mma.tu-darmstadt.de
    Mathematical Modelling and Analysis Group
    Thermo-Fluids and Interfaces
    TU Darmstadt
    Germany

Description
    Given an oriented triangular surface, compute the signed distance of a given
    point or point cloud to the surface.
    The absolute distance is computed as the Euclidean distance between a given
    point and its closest point on the surface.
    The sign (inside/outside) is determined using normal information
    interpolated from the triangle vertices by a linear shape function.

SourceFiles
    triSurfaceDistCalc.cpp
    triSurfaceDistCalcI.H

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceDistCalc_H
#define triSurfaceDistCalc_H

#include "signedDistanceCalculator.H"

#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "vectorField.H"

#include "searchDistanceCalculator.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion
{

/*---------------------------------------------------------------------------*\
                         Class signedDistanceCalculator Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceDistCalc : public signedDistanceCalculator
{
    // Private Data
    //- Calculator to compute the search radii used by the octree search
    searchDistanceCalculator searchDistCalc_;

    //- Triangulated surface
    triSurface surface_;

    //- Octree based search to find closest triangle to a point
    triSurfaceSearch surfaceSearch_;

    //- Vertex normals computed as angle weighted average of adjacent triangles
    vectorField vertexNormals_;


    // Private Member Functions
    //- Compute normals at vertices from triangle normals
    void computeVertexNormals();

    //- Compute signed distances at cell centres and cell corners
    void computeSignedDistances();

public:

    // Static Data Members
    TypeName("triSurface");


    // Constructors
    triSurfaceDistCalc(const dictionary& configDict, const fvMesh& mesh);


    // Member Functions
    //- Compute signed distances for pf and closest point information
    //  Compute signed distances for the points pf. For each point, the
    //  closests point and triangle are also computed. Uses the squared
    //  search distance radii searchDistSqr to determine the narrow band.
    //  Points outisde the narrow band are set to outOfSearchDomain.
    scalarField signedDistance(
        DynamicList<pointIndexHit>& pointToNearestTriangle,
        const pointField& pf,
        const scalarField& searchDistSqr,
        scalar outOfSearchDomain = 0.0) const;

    //- Compute signed distances for points pf.
    //  Compute signed distances to pf in a narrow band defined by the
    //  squared search radii searchDistSqr.
    scalarField signedDistance(const pointField& pf,
        const scalarField& searchDistSqr,
        scalar outOfSearchDomain = 0.0) const;

    //- Return signed distance to p and closest point information.
    std::tuple<pointIndexHit, scalar> signedDistance(
        const point& p, scalar searchDistSqr) const;

    //- Return the normal vector interpolated to a surface point.
    vector normalAtSurface(const pointIndexHit& hitInfo) const;

    //- Reference to the vertex normal field.
    inline const vectorField& vertexNormals() const;

    //- Refrence to the octree search.
    inline const triSurfaceSearch& surfaceSearch() const;

    //- Refrence to the triangulated surface
    inline const triSurface& surface() const;

    //- Return signed distance to p
    scalar signedDistance(const point& p) const override;

    //- Return the minimum edge length of the surface
    scalar referenceLength() const override;

    //- Return the number of surface triangles
    label nSurfaceElements() const override;

    //- Return enclosed volume for *closed* surfaces.
    //  Computes the enclosed volume using surface divergence. The result
    //  only make sense if the surface is *closed*.
    scalar surfaceEnclosedVolume() const override;

    //- Write signed distance and squared search radii fields
    void writeFields() const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triSurfaceDistCalcI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
