/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019 AUTHOR,AFFILIATION
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
    Foam::TriSurfaceImmersion::triAreaFractionCalculator

Description
    Compute area fractions and absolute areas of triangles intersected
    by an interface based on the signed distances of its vertices. The part
    of a triangle that is located on the positive side of the interface
    (defined by the interface normal orientation) contributes. This means
    a triangle with negative signed distances at all its vertices has a
    volume fraction of 0.
    This class implements the model from

    \verbatim
        Detrixhe, M., & Aslam, T. D. (2016).
        From level set to volume of fluid and back again at second‐order
        accuracy.
        International Journal for Numerical Methods in Fluids, 80(4), 231-255.
    \endverbatim

SourceFiles
    triAreaFractionCalculator.cpp

\*---------------------------------------------------------------------------*/

#ifndef triAreaFractionCalculator_H
#define triAreaFractionCalculator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include <array>

#include "AdaptiveTriFaceRefinement.H"


namespace Foam::TriSurfaceImmersion
{

/*---------------------------------------------------------------------------*\
                Class triAreaFractionCalculator Declaration
\*---------------------------------------------------------------------------*/

class triAreaFractionCalculator
{
    // Private data

    //- Temporary storage of vertex signed distances for a triangle
    mutable std::array<scalar, 3> signedDistanceBuffer_;

    // Private Member Functions

    //- Count the number of vertices of a triangle which have negative
    //distances
    label countNegativeDistancesTri() const;


public:
    
    // Member Functions

    //- Compute the triangle's area
    static scalar area(const indexedTri& t, const std::vector<point>& p);

    //- Area fraction of the triangle from signed distances at vertices
    scalar areaFraction
    (
        const indexedTri& tri,
        const std::vector<scalar>& signedDistance
    ) const;

    //- Area of the triangle located on the positive side of the interface
    scalar omegaPlusArea
    (
        const indexedTri& tri,
        const std::vector<scalar>& signedDistance,
        const std::vector<point>& points
    ) const;

    //- Compute area fractions for given triangles
    std::vector<scalar> areaFractions
    (
        const std::vector<indexedTri>& tris,
        const std::vector<scalar>& signedDistance
    ) const;

    //- Compute area on positive side of interface for given triangles
    std::vector<scalar> omegaPlusAreas
    (
        const std::vector<indexedTri>& tris,
        const std::vector<scalar>& signedDistance,
        const std::vector<point>& points
    ) const;

    //- Accumulated area of given triangles which is located on the positive
    //  side of the interface
    scalar accumulatedOmegaPlusArea
    (
        const std::vector<indexedTri>& tris,
        const std::vector<scalar>& signedDistance,
        const std::vector<point>& points
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
