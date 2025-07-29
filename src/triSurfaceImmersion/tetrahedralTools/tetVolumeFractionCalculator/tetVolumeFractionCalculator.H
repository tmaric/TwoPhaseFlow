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
    Foam::TriSurfaceImmersion::tetVolumeFractionCalculator

Description
    Compute volume fractions and absolute volumes of tetrahedra intersected
    by an interface based on the signed distances of its vertices. The part
    of a tetrahedron that is located on the positive side of the interface
    (defined by the interface normal orientation) contributes. This means
    a tetrahedron with negative sidnged distances at all its vertices has a
    volume fraction of 0.
    This class implements the model from

    \verbatim
        Detrixhe, M., & Aslam, T. D. (2016).
        From level set to volume of fluid and back again at second‐order
        accuracy.
        International Journal for Numerical Methods in Fluids, 80(4), 231-255.
    \endverbatim

SourceFiles
    tetVolumeFractionCalculator.cpp

\*---------------------------------------------------------------------------*/

#ifndef tetVolumeFractionCalculator_H
#define tetVolumeFractionCalculator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include <array>

#include "AdaptiveTetCellRefinement.H"


namespace Foam::TriSurfaceImmersion
{

/*---------------------------------------------------------------------------*\
                Class tetVolumeFractionCalculator Declaration
\*---------------------------------------------------------------------------*/

class tetVolumeFractionCalculator
{
    // Private data

    //- Temporary storage of vertex signed distances for a tetrahedron
    mutable std::array<scalar, 4> signedDistanceBuffer_;

    // Private Member Functions

    //- Count the number of vertices of a tetrahedron which have negative
    //distances
    label countNegativeDistances() const;


public:
    
    // Member Functions

    //- Compute the tetrahedron's volume
    static scalar volume(const indexedTet& t, const std::vector<point>& p);

    //- Volume fraction of the tetrahedron from signed distances at vertices
    scalar volumeFraction(
        const indexedTet& tet, const std::vector<scalar>& signedDistance) const;

    //- Volume of the tetrahedron located on the positive side of the interface
    scalar omegaPlusVolume(const indexedTet& tet,
        const std::vector<scalar>& signedDistance,
        const std::vector<point>& points) const;

    //- Compute volume fractions for given tetrahedra
    std::vector<scalar> volumeFractions(const std::vector<indexedTet>& tets,
        const std::vector<scalar>& signedDistance) const;

    //- Compute volume on positive side of interface for given tetrahedra
    std::vector<scalar> omegaPlusVolumes(const std::vector<indexedTet>& tets,
        const std::vector<scalar>& signedDistance,
        const std::vector<point>& points) const;

    //- Accumulated volume of given tetrahedra which is located on the positive
    //  side of the interface
    scalar accumulatedOmegaPlusVolume(const std::vector<indexedTet>& tets,
        const std::vector<scalar>& signedDistance,
        const std::vector<point>& points) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
