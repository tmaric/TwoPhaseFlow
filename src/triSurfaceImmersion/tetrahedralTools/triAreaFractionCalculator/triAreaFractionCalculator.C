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

\*---------------------------------------------------------------------------*/

#include "triAreaFractionCalculator.H"

#include <algorithm>
#include <cassert>
#include <numeric>

namespace Foam::TriSurfaceImmersion
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
label triAreaFractionCalculator::countNegativeDistancesTri() const
{
    return std::count_if(signedDistanceBuffer_.begin(),
        signedDistanceBuffer_.end(),
        [](const scalar d) { return d < 0.0; });
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar triAreaFractionCalculator::area
(
    const indexedTri& t, const std::vector<point>& p
)
{
    return mag((p[t[1]] - p[t[0]]) ^ (p[t[2]] - p[t[0]])) / 2.0;
}


scalar triAreaFractionCalculator::areaFraction
(
    const indexedTri& tri,
    const std::vector<scalar>& signedDistance
) const
{
    // This function implements the actual model of Detrixhe and Aslam (TT)
    scalar areaFraction = 0.0;

    for (int idx = 0; idx != 3; ++idx)
    {
        signedDistanceBuffer_[idx] = signedDistance[tri[idx]];
    }

    std::sort(signedDistanceBuffer_.begin(), signedDistanceBuffer_.end());
    const auto& d = signedDistanceBuffer_;

    auto negative_distances = countNegativeDistancesTri();

    // TODO: continue here with model from Detrixhe and Aslam for triangles
    // Take care to have the correct signs and fractions with respect to
    // +/- distances.
    // Also, check the equations in the paper for plausibility.
    // The paper is located in src/triSurfaceImmersion
    if (negative_distances == 3)
    {
        areaFraction = 0.0;
    }
    else if (negative_distances == 2)
    {
        areaFraction = d[2]*d[2]/((d[2] - d[0])*(d[2] - d[1]));
    }
    else if (negative_distances == 1)
    {
        areaFraction = 1.0 -  d[0]*d[0]/((d[1] - d[0])*(d[2] - d[0]));
    }
    else
    {
        areaFraction = 1.0;
    }

    assert(areaFraction >= 0.0 && areaFraction <= 1.0);
    return areaFraction;
}


scalar triAreaFractionCalculator::omegaPlusArea
(
    const indexedTri& tri,
    const std::vector<scalar>& signedDistance,
    const std::vector<point>& points
) const
{
    return areaFraction(tri, signedDistance) * area(tri, points);
}


std::vector<scalar> triAreaFractionCalculator::areaFractions
(
    const std::vector<indexedTri>& tris,
    const std::vector<scalar>& signedDistance
) const
{
    std::vector<scalar> areaFractions(tris.size());

    forAll(tris, triI)
    {
        areaFractions[triI] = areaFraction(tris[triI], signedDistance);
    }

    return areaFractions;
}


std::vector<scalar> triAreaFractionCalculator::omegaPlusAreas
(
    const std::vector<indexedTri>& tris,
    const std::vector<scalar>& signedDistance,
    const std::vector<point>& points
) const
{
    std::vector<scalar> plus_areas(tris.size());

    forAll(tris, triI)
    {
        plus_areas[triI] =
            omegaPlusArea(tris[triI], signedDistance, points);
    }

    return plus_areas;
}


scalar triAreaFractionCalculator::accumulatedOmegaPlusArea
(
    const std::vector<indexedTri>& tris,
    const std::vector<scalar>& signedDistance,
    const std::vector<point>& points
) const
{
    auto v = omegaPlusAreas(tris, signedDistance, points);
    return std::accumulate(v.begin(), v.end(), 0.0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
