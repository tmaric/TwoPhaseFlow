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

\*---------------------------------------------------------------------------*/

#include "triSurfaceDistCalc.H"

#include "addToRunTimeSelectionTable.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "triSurfaceTools.H"

#include "insideOutsidePropagation.H"

namespace Foam::TriSurfaceImmersion
{

defineTypeNameAndDebug(triSurfaceDistCalc, 0);
addToRunTimeSelectionTable(
    signedDistanceCalculator, triSurfaceDistCalc, Dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void triSurfaceDistCalc::computeVertexNormals()
{
    /* The vertex normals are computed as a weighted sum of normals
     * of the adjacent triangles. The weights are the triangle angles at the
     * vertex at hand.
     * See
     *      "Computing Vertex Normals from Polygonal Facets",
     *      G. Thürmer & C. Wüthrich (2012)
     *      https://doi.org/10.1080/10867651.1998.10487487
     * for details. There is a proof that this gives correct inside/outside
     * information by J. Baerentzen and H. Aanaes, Technical University of
     * Denmark
     */
    const auto& tri_normals = surface_.faceNormals();
    const auto& vertex_to_faces = surface_.pointFaces();
    const auto& vertices = surface_.localPoints();

    forAll(vertices, v_id)
    {
        label vid_a{v_id};
        label vid_b{0};
        label vid_c{0};

        for (const auto fid : vertex_to_faces[v_id])
        {
            triSurfaceTools::otherVertices(surface_, fid, vid_a, vid_b, vid_c);
            vector v1{vertices[vid_b] - vertices[vid_a]};
            vector v2{vertices[vid_c] - vertices[vid_a]};
            scalar alpha{Foam::acos(
                std::clamp((v1 & v2) / (mag(v1) * mag(v2)), -1.0, 1.0))};

            vertexNormals_[v_id] += alpha * tri_normals[fid];
        }

        vertexNormals_[v_id] /= mag(vertexNormals_[v_id]) + SMALL;
    }
}


void triSurfaceDistCalc::computeSignedDistances()
{
    cellSignedDist0_.primitiveFieldRef() = signedDistance(cellNearestTriangle_,
        this->mesh().C(),
        searchDistCalc_.cellSqrSearchDist(),
        this->outOfNarrowBandValue());

    cellSignedDist_ =
        insideOutsidePropagation::propagateInsideOutside(cellSignedDist0_);

    pointSignedDist_.primitiveFieldRef() = signedDistance(pointNearestTriangle_,
        this->mesh().points(),
        searchDistCalc_.pointSqrSearchDist(),
        this->outOfNarrowBandValue());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

triSurfaceDistCalc::triSurfaceDistCalc(
    const dictionary& configDict, const fvMesh& mesh)
    : signedDistanceCalculator{configDict, mesh}, searchDistCalc_{mesh,
                                                      this->narrowBandWidth()},
      surface_{mesh.time().path() + "/" + configDict.get<fileName>("surfaceFile")},
      surfaceSearch_{surface_}, vertexNormals_{
                                    surface_.nPoints(), vector{0, 0, 0}}
{
    computeVertexNormals();
    computeSignedDistances();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalarField triSurfaceDistCalc::signedDistance(
    DynamicList<pointIndexHit>& pointToNearestTriangle,
    const pointField& pf,
    const scalarField& searchDistSqr,
    const scalar outOfSearchDomain) const
{
    scalarField distances(pf.size(), outOfSearchDomain);

    // Use the octree and the square search distance to build the
    // surface mesh / point field proximity information list.
    pointToNearestTriangle.reserve(pf.size());
    surfaceSearch_.findNearest(pf, searchDistSqr, pointToNearestTriangle);

    forAll(pointToNearestTriangle, p_id)
    {
        const pointIndexHit& hitInfo = pointToNearestTriangle[p_id];

        if (hitInfo.hit())
        {
            vector delta_v{pf[p_id] - hitInfo.hitPoint()};
            auto snormal = normalAtSurface(hitInfo);
            distances[p_id] = mag(delta_v) * sign(snormal & delta_v);
        }
    }

    return distances;
}


scalarField triSurfaceDistCalc::signedDistance(const pointField& pf,
    const scalarField& searchDistSqr,
    const scalar outOfSearchDomain) const
{
    DynamicList<pointIndexHit> pointToNearestTriangle{};
    pointToNearestTriangle.reserve(pf.size());

    return signedDistance(
        pointToNearestTriangle, pf, searchDistSqr, outOfSearchDomain);
}


std::tuple<pointIndexHit, scalar> triSurfaceDistCalc::signedDistance(
    const point& p, const scalar searchDistSqr) const
{
    scalar distance{1e15};

    const auto hitInfo =
        surfaceSearch_.nearest(p, vector{Foam::sqrt(searchDistSqr), 0, 0});

    if (hitInfo.hit())
    {
        vector delta_v{p - hitInfo.hitPoint()};
        auto snormal = normalAtSurface(hitInfo);
        distance = mag(delta_v) * sign(snormal & delta_v);
    }

    return std::make_tuple(hitInfo, distance);
}


scalar triSurfaceDistCalc::signedDistance(const point& p) const
{
    return std::get<1>(signedDistance(p, 1e15));
}


vector triSurfaceDistCalc::normalAtSurface(const pointIndexHit& hitInfo) const
{
    vector normal{0, 0, 0};

    const auto& fnormals = surface_.faceNormals();
    const auto& v = surface_.localPoints();

    // Transformation to local coordinate system:
    // - origin: point a (first point of triangle)
    // - x-axis: point a to point b (second point of triangle)
    // - y-axis: cross product of triangle normal and ex
    const auto tri_hit = surface_.localFaces()[hitInfo.index()];
    const auto a_to_b = v[tri_hit[1]] - v[tri_hit[0]];
    const auto a_to_c = v[tri_hit[2]] - v[tri_hit[0]];
    const auto a_to_hit = hitInfo.hitPoint() - v[tri_hit[0]];
    const auto ex = a_to_b / mag(a_to_b);
    const auto ey = fnormals[hitInfo.index()] ^ ex;
    const vector b{a_to_b & ex, 0, 0};
    const vector c{a_to_c & ex, a_to_c & ey, 0};
    const vector h_l{a_to_hit & ex, a_to_hit & ey, 0};

    // This is a linear shape function for a triangle where the vertices are:
    // - 1) the origin (0,0)
    // - 2) a point on the x-axis (t,0)
    // - 3) an arbitrary point (u,v)
    normal = vertexNormals_[tri_hit[0]] *
            (1.0 - h_l.x() / b.x() + h_l.y() / c.y() * (c.x() / b.x() - 1.0)) +
        vertexNormals_[tri_hit[1]] *
            (h_l.x() / b.x() - h_l.y() * c.x() / (b.x() * c.y())) +
        vertexNormals_[tri_hit[2]] * (h_l.y() / c.y());

    return normal;
}


scalar triSurfaceDistCalc::referenceLength() const
{
    // Use minimum edge length as reference length
    const auto& edges = surface_.edges();
    const auto& vertices = surface_.localPoints();

    scalar min_length = 1.0e15;

    for (const auto& anEdge : edges)
    {
        min_length = min(
            min_length, mag(vertices[anEdge.start()] - vertices[anEdge.end()]));
    }

    return min_length;
}


label triSurfaceDistCalc::nSurfaceElements() const
{
    return surface_.size();
}


scalar triSurfaceDistCalc::surfaceEnclosedVolume() const
{
    scalar Vsurf = 0;

    const auto& surfacePoints = surface_.points();
    forAll(surface_, triangleI)
    {
        const auto& Sf = surface_.Sf()[triangleI];
        const auto& triangle = surface_[triangleI];
        Vsurf += dot(-Sf, // Surface normals are oriented into the phase.
            (surfacePoints[triangle[0]] + surfacePoints[triangle[1]] +
                surfacePoints[triangle[2]]));
    }

    Vsurf = 1. / 9. * mag(Vsurf);

    return Vsurf;
}


void triSurfaceDistCalc::writeFields() const
{
    this->signedDistanceCalculator::writeFields();
    searchDistCalc_.writeFields();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
