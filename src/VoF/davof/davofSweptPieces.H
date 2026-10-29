/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
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

Description
    The pieces of a swept volume in the donor cell and its neighbours, each
    clipped by its own cell's PLIC plane, without constructing the cut
    polygons (the working note leia/docs/davof/davof-swept-volume,
    multicell.py, 2026-10-29; the fork's fluxIntegration sweptFaceCells).

    A piece is the swept volume Omega_f (a closed, outward-oriented list of
    planar triangles) intersected with a set of half-spaces: the cell's
    face planes (the hexahedron's six, (x - x0).n <= 0 inside with n the
    outward unit normal) and, for an interface cell, the cell's PLIC
    half-space (liquid where (x - x0).n <= 0). The piece's boundary is the
    clipped boundary triangles of Omega_f plus one cap per plane that cuts
    it. The closure of that boundary,

        sum_T a_T^- + sum_i A_i n_i = 0,

    is three equations for the k cap areas A_i, solvable when k <= 3 and the
    normals are independent (a face-neighbour piece: one cell plane and the
    PLIC plane; an edge neighbour: two and the PLIC plane); x.n_i is
    constant on cap i, so the divergence theorem for x gives

        V = (1/3) [ sum_T p_T . a_T^-  +  sum_i A_i (x0_i . n_i) ].

    The least-squares solution of the overdetermined (k < 3) or square
    (k = 3) system is taken; its residual is the closure residual of the
    piece; a rank deficiency (k > 3, or a PLIC plane parallel to a cell
    face) is reported and the caller falls back. Pieces in FULL cells are
    wet in total (no PLIC plane), in EMPTY cells dry. Verified in Python
    against signed tetrahedra clipped as convex polyhedra with explicit caps
    (data/multicell.csv: 1e-16 relative on the donor, face-neighbour and
    edge-neighbour pieces of planar and warped faces); the fork's
    apps/test/testSweptPieces reproduces the two translation configurations
    of that table.

SourceFiles
    (header only)

\*---------------------------------------------------------------------------*/

#ifndef davofSweptPieces_H
#define davofSweptPieces_H

#include "point.H"
#include "vector.H"
#include "FixedList.H"
#include "List.H"

namespace Foam
{
namespace davof
{

//- A half-space (x - x0).n <= 0, n a unit vector
struct halfSpace
{
    vector n;
    point x0;
};

//- A polygon of at most 16 vertices (a triangle clipped by up to 13
//  planes keeps at most 16)
struct clipPolygon
{
    FixedList<point, 16> p;
    label n = 0;

    inline void setTriangle(const point& a, const point& b, const point& c)
    {
        p[0] = a; p[1] = b; p[2] = c; n = 3;
    }

    //- The area vector (Newell), zero below three vertices
    inline vector areaVector() const
    {
        vector a(Zero);
        for (label i = 0; i < n; ++i)
        {
            a += p[i] ^ p[(i + 1) % n];
        }
        return 0.5*a;
    }
};


//- Clip the polygon by the half-space in place (Sutherland-Hodgman against
//  one plane; a vertex with (x - x0).n <= 0 is kept). Returns true where the
//  plane crossed the polygon (a cap edge exists), false where the polygon
//  was entirely on one side (kept, or removed with n = 0).
inline bool clipByHalfSpace(clipPolygon& poly, const halfSpace& hs)
{
    if (poly.n < 3) return false;
    FixedList<scalar, 16> s;
    bool anyIn = false, anyOut = false;
    for (label i = 0; i < poly.n; ++i)
    {
        s[i] = (poly.p[i] - hs.x0) & hs.n;
        if (s[i] <= 0) anyIn = true; else anyOut = true;
    }
    if (!anyOut) return false;
    if (!anyIn) { poly.n = 0; return false; }
    FixedList<point, 16> out;
    label m = 0;
    for (label i = 0; i < poly.n; ++i)
    {
        const label j = (i + 1) % poly.n;
        if (s[i] <= 0) out[m++] = poly.p[i];
        if ((s[i] <= 0) != (s[j] <= 0))
        {
            out[m++] = poly.p[i] + (s[i]/(s[i] - s[j]))*(poly.p[j] - poly.p[i]);
        }
    }
    poly.p = out;
    poly.n = (m >= 3) ? m : 0;
    return true;
}


//- The result of one piece
struct pieceResult
{
    //- The liquid volume of the piece (signed with the list's orientation)
    scalar volume = 0;

    //- The closure residual |N A + sum a| over sum |a_T| of the wet
    //  boundary (round-off for a consistent cap system)
    scalar closureResidual = 0;

    //- The number of cap planes that cut the piece, and whether their
    //  system was solvable (k <= 3, independent normals)
    label nCaps = 0;
    bool ok = true;

    //- The coverage: the total area of the boundary triangles' parts inside
    //  the cell planes (before the PLIC clip), for the check that the
    //  pieces of the ring cells cover the whole boundary
    scalar coverage = 0;

    //- The clipped triangles (partially cut by any plane)
    label nCut = 0;

    //- A degenerate piece: boundary triangles that lie IN a cell plane
    //  (the face's own fan in the plane of the downwind neighbour's face,
    //  a lateral triangle in a side plane under a translation normal to
    //  the face) pass the inside test with a zero signed distance and form
    //  a slab of zero thickness whose boundary does not close (no cap was
    //  cut); its volume is zero and it is counted
    bool degenerate = false;
};


//- The liquid volume of Omega_f (the triangle list T, outward) inside the
//  cell given by cellPlanes, clipped by the PLIC half-space plic when
//  wetPlane is true (an interface cell), wet in total otherwise (a full
//  cell). xRef is the origin of the moments (any point; the face centre in
//  the flux, so that a slab in the face plane carries no moment even
//  before the degenerate rule). See the file header.
inline pieceResult sweptPieceVolume
(
    const UList<FixedList<point, 3>>& T,
    const UList<halfSpace>& cellPlanes,
    const bool wetPlane,
    const halfSpace& plic,
    const point& xRef
)
{
    pieceResult r;
    const label nc = cellPlanes.size();
    if (nc > 7)
    {
        // more cell planes than the cap bookkeeping holds (a polyhedron):
        // reported as not solvable, the caller falls back
        r.ok = false;
        return r;
    }
    FixedList<bool, 8> cut(false);
    bool plicCut = false;
    vector sumA(Zero);
    scalar moment = 0, sumMag = 0;
    forAll(T, k)
    {
        clipPolygon poly;
        poly.setTriangle(T[k][0], T[k][1], T[k][2]);
        bool wasCut = false;
        for (label i = 0; i < nc; ++i)
        {
            if (clipByHalfSpace(poly, cellPlanes[i])) { cut[i] = true; wasCut = true; }
            if (poly.n < 3) break;
        }
        if (poly.n < 3) continue;
        r.coverage += mag(poly.areaVector());
        if (wetPlane)
        {
            if (clipByHalfSpace(poly, plic)) { plicCut = true; wasCut = true; }
            if (poly.n < 3) continue;
        }
        if (wasCut) ++r.nCut;
        const vector a = poly.areaVector();
        sumA += a;
        sumMag += mag(a);
        moment += (poly.p[0] - xRef) & a;
    }

    // the cap system N A = -sum a, N the 3 x k matrix of the cap normals
    FixedList<vector, 8> N;
    FixedList<scalar, 8> d;   // x0_i . n_i
    label k = 0;
    for (label i = 0; i < nc; ++i)
    {
        if (cut[i]) { N[k] = cellPlanes[i].n; d[k] = (cellPlanes[i].x0 - xRef) & cellPlanes[i].n; ++k; }
    }
    if (plicCut) { N[k] = plic.n; d[k] = (plic.x0 - xRef) & plic.n; ++k; }
    r.nCaps = k;
    if (sumMag <= VSMALL)
    {
        r.volume = 0;
        return r;
    }
    if (k == 0 && mag(sumA) > 1e-9*sumMag)
    {
        // no cap and no closure: a slab of zero thickness in a cell plane
        r.degenerate = true;
        r.closureResidual = mag(sumA)/sumMag;
        r.volume = 0;
        return r;
    }
    if (k > 3)
    {
        r.ok = false;
        r.closureResidual = mag(sumA)/sumMag;
        return r;
    }
    FixedList<scalar, 3> A(scalar(0));
    if (k > 0)
    {
        // the normal equations G A = -N^T sum a, G = N^T N (k x k)
        scalar G[3][3], rhs[3];
        for (label i = 0; i < k; ++i)
        {
            rhs[i] = -(N[i] & sumA);
            for (label j = 0; j < k; ++j) G[i][j] = N[i] & N[j];
        }
        // the Gram determinant against its scale (unit normals): a
        // rank deficiency below 1e-10
        scalar det = 1;
        if (k == 1) det = G[0][0];
        else if (k == 2) det = G[0][0]*G[1][1] - G[0][1]*G[1][0];
        else
        {
            det = G[0][0]*(G[1][1]*G[2][2] - G[1][2]*G[2][1])
                - G[0][1]*(G[1][0]*G[2][2] - G[1][2]*G[2][0])
                + G[0][2]*(G[1][0]*G[2][1] - G[1][1]*G[2][0]);
        }
        if (mag(det) < 1e-10)
        {
            r.ok = false;
            r.closureResidual = mag(sumA)/sumMag;
            return r;
        }
        if (k == 1)
        {
            A[0] = rhs[0]/G[0][0];
        }
        else if (k == 2)
        {
            A[0] = (rhs[0]*G[1][1] - rhs[1]*G[0][1])/det;
            A[1] = (G[0][0]*rhs[1] - G[1][0]*rhs[0])/det;
        }
        else
        {
            // Cramer's rule
            for (label c = 0; c < 3; ++c)
            {
                scalar M[3][3];
                for (label i = 0; i < 3; ++i)
                    for (label j = 0; j < 3; ++j)
                        M[i][j] = (j == c) ? rhs[i] : G[i][j];
                A[c] =
                (
                    M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
                  - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
                  + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0])
                )/det;
            }
        }
    }
    vector res = sumA;
    scalar capTerm = 0;
    for (label i = 0; i < k; ++i)
    {
        res += A[i]*N[i];
        capTerm += A[i]*d[i];
    }
    r.closureResidual = mag(res)/sumMag;
    if (r.closureResidual > 1e-6)
    {
        // the cap system does not close the piece: a slab with a spurious
        // cut, or an inconsistent list; nothing is taken from it
        r.degenerate = true;
        r.volume = 0;
        return r;
    }
    r.volume = (moment + capTerm)/3;
    return r;
}


//- The wetted area vector of the pre-image fan triangles T (in the face's
//  orientation) inside the cell planes, clipped by the PLIC half-space when
//  wetPlane; coverage as above
inline vector sweptPieceArea
(
    const UList<FixedList<point, 3>>& T,
    const UList<halfSpace>& cellPlanes,
    const bool wetPlane,
    const halfSpace& plic,
    scalar& coverage
)
{
    vector a(Zero);
    coverage = 0;
    forAll(T, k)
    {
        clipPolygon poly;
        poly.setTriangle(T[k][0], T[k][1], T[k][2]);
        forAll(cellPlanes, i)
        {
            clipByHalfSpace(poly, cellPlanes[i]);
            if (poly.n < 3) break;
        }
        if (poly.n < 3) continue;
        coverage += mag(poly.areaVector());
        if (wetPlane)
        {
            clipByHalfSpace(poly, plic);
            if (poly.n < 3) continue;
        }
        a += poly.areaVector();
    }
    return a;
}

} // End namespace davof
} // End namespace Foam

#endif

// ************************************************************************* //
