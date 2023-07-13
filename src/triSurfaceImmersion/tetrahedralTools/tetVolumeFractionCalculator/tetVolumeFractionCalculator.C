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

#include "tetVolumeFractionCalculator.H"

#include <algorithm>
#include <cassert>
#include <numeric>

namespace Foam::TriSurfaceImmersion
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
label tetVolumeFractionCalculator::countNegativeDistances() const
{
    return std::count_if(signedDistanceBuffer_.begin(),
        signedDistanceBuffer_.end(),
        [](const scalar d) { return d < 0.0; });
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar tetVolumeFractionCalculator::volume(
    const indexedTet& t, const std::vector<point>& p)
{
    return mag((p[t[1]] - p[t[0]]) &
               ((p[t[2]] - p[t[0]]) ^ (p[t[3]] - p[t[0]]))) /
        6.0;
}


scalar tetVolumeFractionCalculator::volumeFraction(
    const indexedTet& tet, const std::vector<scalar>& signedDistance) const
{
    // This function implements the actual model of Detrixhe and Aslam (TT)
    scalar volFraction = 0.0;

    for (int idx = 0; idx != 4; ++idx)
    {
        signedDistanceBuffer_[idx] = signedDistance[tet[idx]];
    }

    std::sort(signedDistanceBuffer_.begin(), signedDistanceBuffer_.end());
    const auto& d = signedDistanceBuffer_;

    auto negative_distances = countNegativeDistances();

    if (negative_distances == 4)
    {
        volFraction = 0.0;
    }
    else if (negative_distances == 3)
    {
        volFraction =
            std::pow(d[3], 3) / ((d[3] - d[0]) * (d[3] - d[1]) * (d[3] - d[2]));
    }
    else if (negative_distances == 2)
    {
        volFraction =
            (d[0] * d[1] * (d[2] * d[2] + d[2] * d[3] + d[3] * d[3]) +
                d[2] * d[3] * (d[2] * d[3] - (d[0] + d[1]) * (d[2] + d[3]))) /
            ((d[0] - d[2]) * (d[1] - d[2]) * (d[0] - d[3]) * (d[1] - d[3]));
    }
    else if (negative_distances == 1)
    {
        volFraction = 1.0 +
            std::pow(d[0], 3) / ((d[1] - d[0]) * (d[2] - d[0]) * (d[3] - d[0]));
    }
    else
    {
        volFraction = 1.0;
    }

    assert(volFraction >= 0.0 && volFraction <= 1.0);
    return volFraction;
}


scalar tetVolumeFractionCalculator::omegaPlusVolume(const indexedTet& tet,
    const std::vector<scalar>& signedDistance,
    const std::vector<point>& points) const
{
    return volumeFraction(tet, signedDistance) * volume(tet, points);
}


std::vector<scalar> tetVolumeFractionCalculator::volumeFractions(
    const std::vector<indexedTet>& tets,
    const std::vector<scalar>& signedDistance) const
{
    std::vector<scalar> volumeFractions(tets.size());

    forAll(tets, tetI)
    {
        volumeFractions[tetI] = volumeFraction(tets[tetI], signedDistance);
    }

    return volumeFractions;
}


std::vector<scalar> tetVolumeFractionCalculator::omegaPlusVolumes(
    const std::vector<indexedTet>& tets,
    const std::vector<scalar>& signedDistance,
    const std::vector<point>& points) const
{
    std::vector<scalar> plus_volumes(tets.size());

    forAll(tets, tetI)
    {
        plus_volumes[tetI] =
            omegaPlusVolume(tets[tetI], signedDistance, points);
    }

    return plus_volumes;
}


scalar tetVolumeFractionCalculator::accumulatedOmegaPlusVolume(
    const std::vector<indexedTet>& tets,
    const std::vector<scalar>& signedDistance,
    const std::vector<point>& points) const
{
    auto v = omegaPlusVolumes(tets, signedDistance, points);
    return std::accumulate(v.begin(), v.end(), 0.0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
