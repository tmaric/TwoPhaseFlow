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

    THE CAP SET, ORDER-INDEPENDENT (2026-11-02). Which planes carry a cap
    decides k and so whether the caller falls back to the one-donor clip of
    the whole face. Until 2026-11-01 a plane counted as cut when it crossed
    the polygon at the moment it was clipped, i.e. after the planes before
    it in the cell's face list; a plane whose new edge a later plane removed
    again still counted, and the clipped remnant of zero area (a triangle
    touching a cell plane in one edge) set a flag as well. The cap set then
    depended on the order of the cell's faces, which differs between the
    serial and a decomposed mesh: the N = 32 sixteen-rank run (fork
    davof-swept-par, job 55812211) took the pieces for one face at step 20
    where serial fell back, through the zero-volume piece of an edge
    neighbour whose boundary triangles lay in that neighbour's face plane.
    Two rules make the cap set a property of the piece itself, i.e. of the
    set Omega_f intersected with the half-spaces, not of the clipping order:
    (a) a boundary triangle lying IN a cell plane (every vertex within 1e-10
        of its longest edge) is a contact, not a boundary of the piece, when
        Omega_f lies outside the cell there: it is dropped; with Omega_f
        inside it bounds the piece and is not clipped by that plane. The
        side is the DONOR's ((x_donor - x0_i).n_i > 0: outside), as in the
        coplanar rule of the face update (sweptPieceArea, donorSide): such
        triangles arise where the motion is parallel to a cell plane that
        holds the face or an edge of it (the face's own fan, the lateral
        quads of the edges in a side plane), and the swept volume then lies
        on the face's, i.e. the donor's, side of that plane. The triangle's
        own orientation (a_T.n_i < 0: outside) is not used for this: on a
        twisted swept volume a lateral quad in a side plane is a bow tie
        whose two triangles point to opposite sides of the same plane. This
        replaces the zero signed distance, which accepted such a triangle in
        the cells on both sides of the plane;
    (b) a plane carries a cap when the FINAL wet polygon of some triangle
        (after every cell plane and the PLIC plane) has an edge that this
        plane produced, longer than 1e-12 of the triangle's longest edge, on
        a polygon whose area exceeds 1e-14 of that edge squared (round-off;
        a remnant below it is dropped with its area). The edge labels are
        carried through the clipping (clipByHalfSpace): an exit point starts
        the edge along the plane, an entry point continues the input edge.
    In exact arithmetic the final polygons and their plane-produced edges
    are those of the intersection, whatever the order; the remaining
    switches are the round-off ties of a genuinely degenerate geometry. A
    piece whose triangles are all contacts has no coverage and no rank
    status, so a zero-volume piece no longer decides the fallback.

    Polygon clipping by a plane: I. E. Sutherland, G. W. Hodgman, Reentrant
    polygon clipping, Communications of the ACM 17 (1974) 32-42,
    doi:10.1145/360767.360802 (one plane at a time, with the edge labels
    added here).

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

//- A half-space (x - x0).n <= 0, n a unit vector. axis >= 0 marks a plane
//  whose normal has exactly two zero components (an axis-aligned cell face
//  on a hexahedral mesh; key sweptHexClip, 2026-11-05): its signed distance
//  is the one coordinate difference times that component, which is bit for
//  bit the dot product (the two other products are exact zeros and
//  0 + 0 + x = x in IEEE arithmetic), without the two multiplications and
//  two additions. -1 (the default): the dot product.
struct halfSpace
{
    vector n;
    point x0;
    label axis = -1;
};

//- The axis of an axis-aligned normal (exactly two zero components), else -1
inline label axisOf(const vector& n)
{
    if (n.y() == 0 && n.z() == 0 && n.x() != 0) return 0;
    if (n.x() == 0 && n.z() == 0 && n.y() != 0) return 1;
    if (n.x() == 0 && n.y() == 0 && n.z() != 0) return 2;
    return -1;
}

//- The signed distance (x - x0).n of the half-space (the axis path where set)
inline scalar sdist(const halfSpace& hs, const point& p)
{
    if (hs.axis >= 0)
    {
        return (p[hs.axis] - hs.x0[hs.axis])*hs.n[hs.axis];
    }
    return (p - hs.x0) & hs.n;
}

//- The half-space clips performed (every clipWith call: one polygon against
//  one plane), a counter the caller resets per step and reports (the cost
//  measure of the clipping rules, key sweptClip, 2026-11-05)
inline label& clipCount()
{
    static label n = 0;
    return n;
}

//- A polygon of at most 16 vertices (a triangle clipped by up to 13
//  planes keeps at most 16)
struct clipPolygon
{
    FixedList<point, 16> p;

    //- The label of the edge p[i] -> p[i+1]: -1 for a part of an input
    //  triangle's edge, otherwise the id of the plane that produced it
    //  (clipByHalfSpace; 2026-11-02, the order-independent cap set)
    FixedList<label, 16> e;
    label n = 0;

    inline void setTriangle(const point& a, const point& b, const point& c)
    {
        p[0] = a; p[1] = b; p[2] = c; n = 3;
        e[0] = e[1] = e[2] = -1;
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
//  plane crossed the polygon, false where the polygon was entirely on one
//  side (kept, or removed with n = 0). The edge labels follow the clip: a
//  kept vertex keeps its edge's label, an exit point (p[i] in, p[i+1] out)
//  starts the new edge along the plane (label id), an entry point (p[i]
//  out, p[i+1] in) continues the rest of edge i (its label). Whether the
//  plane bounds the final polygon is read from those labels after all
//  clips (sweptPieceVolume), not from the return value.
//  Since 2026-11-04 the clip is in two steps with the same arithmetic: the
//  signed distances of the vertices (planeDistances, shared with the
//  contact and coplanar tests that need them too) and the clip into a
//  second polygon (clipWith; the caller alternates two buffers instead of
//  copying the result back), so that the result is bit for bit the one of
//  the single-step clip of 2026-11-02.
inline void planeDistances
(
    const clipPolygon& poly,
    const halfSpace& hs,
    FixedList<scalar, 16>& s,
    bool& anyIn,
    bool& anyOut,
    scalar& sMaxAbs
)
{
    anyIn = false; anyOut = false; sMaxAbs = 0;
    for (label i = 0; i < poly.n; ++i)
    {
        s[i] = sdist(hs, poly.p[i]);
        if (s[i] <= 0) anyIn = true; else anyOut = true;
        sMaxAbs = max(sMaxAbs, mag(s[i]));
    }
}

//- The clip of poly by the plane with the distances s (anyIn and anyOut
//  from planeDistances) into out. Returns true where the plane crossed the
//  polygon (out holds the result, with n = 0 below three vertices), false
//  where poly was entirely on one side: kept (unchanged) or removed (n set
//  to 0 in poly).
inline bool clipWith
(
    clipPolygon& poly,
    const FixedList<scalar, 16>& s,
    const bool anyIn,
    const bool anyOut,
    clipPolygon& out,
    const label id = -1
)
{
    ++clipCount();
    if (poly.n < 3) return false;
    if (!anyOut) return false;
    if (!anyIn) { poly.n = 0; return false; }
    label m = 0;
    for (label i = 0; i < poly.n; ++i)
    {
        const label j = (i + 1) % poly.n;
        const bool inI = (s[i] <= 0);
        const bool inJ = (s[j] <= 0);
        if (inI)
        {
            out.p[m] = poly.p[i];
            out.e[m] = poly.e[i];
            ++m;
        }
        if (inI != inJ)
        {
            out.p[m] = poly.p[i] + (s[i]/(s[i] - s[j]))*(poly.p[j] - poly.p[i]);
            out.e[m] = inI ? id : poly.e[i];
            ++m;
        }
    }
    out.n = (m >= 3) ? m : 0;
    return true;
}

//- The single-step clip (the distances and the clip, the result back in
//  poly); kept for reference, the routines below use the two steps
inline bool clipByHalfSpace
(
    clipPolygon& poly,
    const halfSpace& hs,
    const label id = -1
)
{
    if (poly.n < 3) return false;
    FixedList<scalar, 16> s;
    bool anyIn, anyOut; scalar sMaxAbs;
    planeDistances(poly, hs, s, anyIn, anyOut, sMaxAbs);
    clipPolygon out;
    if (!clipWith(poly, s, anyIn, anyOut, out, id)) return false;
    poly = out;
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

    //- The boundary triangles dropped as contacts: lying in a cell plane
    //  with Omega_f (the donor) outside the cell (the file header, rule a)
    label nContact = 0;

    //- A degenerate piece: boundary triangles that lie IN a cell plane
    //  (the face's own fan in the plane of the downwind neighbour's face,
    //  a lateral triangle in a side plane under a translation normal to
    //  the face) pass the inside test with a zero signed distance and form
    //  a slab of zero thickness whose boundary does not close (no cap was
    //  cut); its volume is zero and it is counted. Since 2026-11-02 such
    //  triangles are contacts (rule a) where Omega_f lies outside the cell
    //  and never enter the piece; the flag remains for a closure that fails
    //  otherwise
    bool degenerate = false;
};


//- The liquid volume of Omega_f (the triangle list T, outward) inside the
//  cell given by cellPlanes, clipped by the PLIC half-space plic when
//  wetPlane is true (an interface cell), wet in total otherwise (a full
//  cell). xRef is the origin of the moments (any point; the face centre in
//  the flux, so that a slab in the face plane carries no moment even
//  before the degenerate rule). xSide (the donor's centre in the flux) is
//  the side of a cell plane on which a triangle lying in that plane bounds
//  Omega_f (rule a); without it the triangle's own orientation decides.
//  crossing true restores the cap set before 2026-11-02 (a plane counts
//  when it crosses the polygon at its turn in the clipping order, no
//  contact rule, no remnant filter; fork key sweptCapSet crossing, kept for
//  the comparison). See the file header.
inline pieceResult sweptPieceVolume
(
    const UList<FixedList<point, 3>>& T,
    const UList<halfSpace>& cellPlanes,
    const bool wetPlane,
    const halfSpace& plic,
    const point& xRef,
    const point* xSide = nullptr,
    const bool crossing = false
)
{
    pieceResult r;
    const label nc = cellPlanes.size();
    if (nc > 13)
    {
        // more planes than the cap bookkeeping holds (a polyhedron; the
        // polygon buffer holds a triangle clipped by 13 planes): reported
        // as not solvable, the caller falls back. Until 2026-11-05 the bound
        // was 7 (a hexahedron's six planes); the joint clip (sweptClip
        // hybrid) adds the touched interface cells' PLIC planes to a piece's
        // cell planes.
        r.ok = false;
        return r;
    }
    // the cap planes: cut[i] for the cell plane i, plicCut for the PLIC
    // plane (id nc), from the edges of the final wet polygons (the file
    // header, rule b)
    FixedList<bool, 16> cut(false);
    bool plicCut = false;
    vector sumA(Zero);
    scalar moment = 0, sumMag = 0;
    // the two polygon buffers of the clips (alternated, no copy back)
    clipPolygon bufA, bufB;
    FixedList<scalar, 16> s;
    forAll(T, k)
    {
        const FixedList<point, 3>& t = T[k];
        const scalar L =
            max(mag(t[1] - t[0]), max(mag(t[2] - t[1]), mag(t[0] - t[2])));
        if (!(L > 0)) continue;
        clipPolygon* cur = &bufA;
        clipPolygon* nxt = &bufB;
        cur->setTriangle(t[0], t[1], t[2]);
        bool anyIn, anyOut; scalar sMaxAbs;
        // one clip of the current polygon by hs into the other buffer
        auto clip = [&](const halfSpace& hs, const label id) -> bool
        {
            planeDistances(*cur, hs, s, anyIn, anyOut, sMaxAbs);
            if (clipWith(*cur, s, anyIn, anyOut, *nxt, id)) { Swap(cur, nxt); return true; }
            return false;
        };

        if (crossing)
        {
            // the cap set before 2026-11-02 (order-dependent)
            bool wasCut = false;
            for (label i = 0; i < nc; ++i)
            {
                if (clip(cellPlanes[i], -1)) { cut[i] = true; wasCut = true; }
                if (cur->n < 3) break;
            }
            if (cur->n < 3) continue;
            r.coverage += mag(cur->areaVector());
            if (wetPlane)
            {
                if (clip(plic, -1)) { plicCut = true; wasCut = true; }
                if (cur->n < 3) continue;
            }
            if (wasCut) ++r.nCut;
            const vector a = cur->areaVector();
            sumA += a;
            sumMag += mag(a);
            moment += (cur->p[0] - xRef) & a;
            continue;
        }

        // rule a: a triangle lying in a cell plane is a contact (dropped)
        // where Omega_f lies outside the cell, and is not clipped by that
        // plane where Omega_f lies inside; the side is the donor's (xSide),
        // else the triangle's orientation. The distances of the three
        // vertices to every plane are kept (s3): they are the distances the
        // first clips would compute again (the same expression) while the
        // polygon is still the triangle, and they settle a plane the whole
        // triangle lies beyond (removed) or within (kept) with a margin of
        // 1e-12 L that no crossing point's round-off reaches (2026-11-04).
        const vector aT = cur->areaVector();
        FixedList<bool, 16> inPlane(false);
        FixedList<FixedList<scalar, 3>, 16> s3;
        FixedList<scalar, 16> s3min, s3max;
        bool contact = false;
        bool beyond = false;
        for (label i = 0; i < nc; ++i)
        {
            scalar sMax = 0;
            s3min[i] = GREAT; s3max[i] = -GREAT;
            for (label j = 0; j < 3; ++j)
            {
                s3[i][j] = sdist(cellPlanes[i], t[j]);
                sMax = max(sMax, mag(s3[i][j]));
                s3min[i] = min(s3min[i], s3[i][j]); s3max[i] = max(s3max[i], s3[i][j]);
            }
            if (sMax <= 1e-10*L)
            {
                const scalar sSide =
                    xSide ? sdist(cellPlanes[i], *xSide) : 0;
                const bool outside =
                    (mag(sSide) > 1e-10*L) ? (sSide > 0) : ((aT & cellPlanes[i].n) < 0);
                if (outside) { contact = true; break; }
                inPlane[i] = true;
            }
            else if (s3min[i] > 1e-12*L)
            {
                // the triangle lies beyond this plane: its clip, or an
                // earlier plane's, leaves nothing of it
                beyond = true;
            }
        }
        if (contact)
        {
            ++r.nContact;
            continue;
        }
        if (beyond) continue;

        bool wasCut = false;
        bool pristine = true;   // the polygon is still the triangle
        for (label i = 0; i < nc; ++i)
        {
            if (inPlane[i]) continue;
            if (pristine)
            {
                // the triangle's own distances (the clip's values)
                anyIn = anyOut = false;
                for (label j = 0; j < 3; ++j)
                {
                    s[j] = s3[i][j];
                    if (s[j] <= 0) anyIn = true; else anyOut = true;
                }
                if (clipWith(*cur, s, anyIn, anyOut, *nxt, i)) { Swap(cur, nxt); wasCut = true; pristine = false; }
            }
            else if (s3max[i] < -1e-12*L)
            {
                continue;   // the triangle lies within this plane: kept
            }
            else if (clip(cellPlanes[i], i))
            {
                wasCut = true;
            }
            if (cur->n < 3) break;
        }
        if (cur->n < 3) continue;
        // a remnant of zero area (a triangle touching the cell in an edge
        // or a point) is dropped, edges and all
        const scalar areaTol = 1e-14*L*L;
        const vector ac = cur->areaVector();
        if (mag(ac) <= areaTol) continue;
        r.coverage += mag(ac);
        if (wetPlane)
        {
            if (clip(plic, nc)) wasCut = true;
            if (cur->n < 3) continue;
        }
        const clipPolygon& poly = *cur;
        const vector a = poly.areaVector();
        if (mag(a) <= areaTol) continue;

        // rule b: the planes that produced an edge of the final polygon
        for (label j = 0; j < poly.n; ++j)
        {
            const label id = poly.e[j];
            if (id < 0) continue;
            if (mag(poly.p[(j + 1) % poly.n] - poly.p[j]) <= 1e-12*L) continue;
            if (id == nc) plicCut = true; else cut[id] = true;
        }
        if (wasCut) ++r.nCut;
        sumA += a;
        sumMag += mag(a);
        moment += (poly.p[0] - xRef) & a;
    }

    // the cap system N A = -sum a, N the 3 x k matrix of the cap normals
    FixedList<vector, 16> N;
    FixedList<scalar, 16> d;   // x0_i . n_i
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
//  wetPlane; coverage as above.
//
//  A COPLANAR part (a polygon whose vertices all lie within 1e-10 of its
//  size of a cell plane: the pre-image of a face parallel to the velocity
//  lies in the face's own plane, which is a face plane of the two cells on
//  either side of it) passes the inside test of BOTH cells with a zero
//  signed distance, so without a rule it is counted twice, once under each
//  cell's PLIC plane (the defect found 2026-10-29: the face fraction of
//  every face parallel to the velocity was the sum of two fractions). With
//  nRef non-zero a coplanar part belongs to the cell whose outward normal
//  on that plane points along nRef (the donor's side: nRef is the donor's
//  outward normal through the face) and is dropped for the cell across;
//  with nRef zero it is kept in both (the defective rule, kept for the
//  comparison). aCell returns the part inside the cell planes before the
//  PLIC clip; coplanarHit whether a coplanar part was met.
inline vector sweptPieceArea
(
    const UList<FixedList<point, 3>>& T,
    const UList<halfSpace>& cellPlanes,
    const bool wetPlane,
    const halfSpace& plic,
    scalar& coverage,
    const vector& nRef = vector::zero,
    vector* aCell = nullptr,
    bool* coplanarHit = nullptr
)
{
    vector a(Zero);
    coverage = 0;
    if (aCell) *aCell = Zero;
    if (coplanarHit) *coplanarHit = false;
    const bool tie = (magSqr(nRef) > 0);
    clipPolygon bufA, bufB;
    FixedList<scalar, 16> s;
    forAll(T, k)
    {
        clipPolygon* cur = &bufA;
        clipPolygon* nxt = &bufB;
        cur->setTriangle(T[k][0], T[k][1], T[k][2]);
        // the longest edge of the triangle: the bound of the size L of every
        // polygon the clips leave (2026-11-04)
        const scalar Le =
            max(mag(T[k][1] - T[k][0]), max(mag(T[k][2] - T[k][1]), mag(T[k][0] - T[k][2])));
        bool anyIn, anyOut; scalar sMaxAbs;
        forAll(cellPlanes, i)
        {
            // the coplanar test: every vertex within 1e-10 of the polygon's
            // size of the plane; the distances serve the clip below. A
            // polygon whose largest distance exceeds 1e-10 Le is not
            // coplanar (its size is at most Le); only the others measure
            // their size
            planeDistances(*cur, cellPlanes[i], s, anyIn, anyOut, sMaxAbs);
            if (cur->n >= 3 && sMaxAbs <= 1e-10*Le)
            {
                scalar L = 0;
                for (label j = 0; j < cur->n; ++j)
                {
                    L = max(L, mag(cur->p[j] - cur->p[0]));
                }
                if (L > 0 && sMaxAbs <= 1e-10*L)
                {
                    if (coplanarHit) *coplanarHit = true;
                    if (tie && (cellPlanes[i].n & nRef) <= 0)
                    {
                        cur->n = 0;     // the part belongs to the cell across
                        break;
                    }
                    continue;           // inside (the donor's side, or no rule)
                }
            }
            if (clipWith(*cur, s, anyIn, anyOut, *nxt, -1)) Swap(cur, nxt);
            if (cur->n < 3) break;
        }
        if (cur->n < 3) continue;
        const vector ac = cur->areaVector();
        coverage += mag(ac);
        if (aCell) *aCell += ac;
        if (wetPlane)
        {
            planeDistances(*cur, plic, s, anyIn, anyOut, sMaxAbs);
            if (clipWith(*cur, s, anyIn, anyOut, *nxt, -1)) Swap(cur, nxt);
            if (cur->n < 3) continue;
        }
        a += cur->areaVector();
    }
    return a;
}

} // End namespace davof
} // End namespace Foam

#endif

// ************************************************************************* //
