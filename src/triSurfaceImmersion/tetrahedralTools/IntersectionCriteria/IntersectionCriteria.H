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
    Foam::TriSurfaceImmersion::boundingBallCriterion
    Foam::TriSurfaceImmersion::signCriterion

Description
    Criteria to decide whether a cell is intersected by an interface based
    on the level set values or the signed distances defined at the cell
    corner points. 
    The classes themselves are only used for tag dispatch. The actual
    criteria are implemented in the overloaded functions
    'considerIntersected'.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef intersectionCriteria_H
#define intersectionCriteria_H

#include "fvCFD.H"

#include <algorithm>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion
{

// Structs used for tag dispatch
struct boundingBallCriterion{};
struct signCriterion{};

// boundingBallCriterion implementation
template<class IndexedPolyhedron,
    class PointContainer,
    class LevelSetValueContainer>
bool considerIntersected(const point& refPoint,
    scalar maxDistSqr,
    const IndexedPolyhedron& poly,
    const PointContainer& points,
    const LevelSetValueContainer& values,
    const boundingBallCriterion& criterion)
{
    return std::any_of(poly.begin(), poly.end(),
        [&](label pID){return Foam::magSqr(points[pID] - refPoint) >= maxDistSqr;}
    );
}


template<class IndexSet, class ValueContainer>
std::tuple<scalar,label> maxDistSqrAndPointID(const IndexSet& ids,
    const ValueContainer& values)
{
    label maxID = *(std::max_element(ids.begin(), ids.end(),
        [&](label a, label b){
            return (Foam::magSqr(values[a]) < Foam::magSqr(values[b]));}
    ));

    return std::make_tuple(Foam::magSqr(values[maxID]), maxID);
}


template<class IndexedPolyhedron,
    class PointContainer,
    class LevelSetValueContainer>
bool considerIntersected(const IndexedPolyhedron& poly,
    const PointContainer& points,
    const LevelSetValueContainer& values,
    const boundingBallCriterion& criterion)
{
    const auto [maxDistSqr, pID] = maxDistSqrAndPointID(poly, values);

    return considerIntersected(points[pID],
        maxDistSqr, poly, points, values, criterion);
}


template<class LevelSet>
scalar levelSetValue(const LevelSet& ls,
    const point& p,
    const boundingBallCriterion& criterion)
{
    return ls.signedDistance(p);
}


// signCriterion implementation 
template<class IndexedPolyhedron,
    class PointContainer,
    class LevelSetValueContainer>
bool considerIntersected(const point& refPoint,
    scalar refValue,
    const IndexedPolyhedron& poly,
    const PointContainer& points,
    const LevelSetValueContainer& values,
    const signCriterion& criterion)
{
    return std::any_of(poly.begin(), poly.end(),
        [&](label pID){return values[pID]*refValue < 0.0;}
    );
}


template<class IndexedPolyhedron,
    class PointContainer,
    class LevelSetValueContainer>
bool considerIntersected(const IndexedPolyhedron& poly,
    const PointContainer& points,
    const LevelSetValueContainer& values,
    const signCriterion& criterion)
{
    return considerIntersected(points[poly[0]],
        values[poly[0]], poly, points, values, criterion);
}


template<class LevelSet>
scalar levelSetValue(const LevelSet& ls,
    const point& p,
    const signCriterion& criterion)
{
    return ls.value(p);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //