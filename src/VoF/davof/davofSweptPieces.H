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

//- A half-space (x - x0).n <= 0, n a unit vector. An axis-aligned plane
//  (two components of n exactly zero, the third +-1; the faces of a
//  hexahedral mesh) carries its axis and sign (setAxis, key
//  sweptAxisPlanes, 2026-11-05): its signed distance is then the one
//  coordinate difference, sgn (p[a] - x0[a]), which is bit for bit the dot
//  product (p - x0).n (the two zero products and their sums vanish exactly,
//  the product by +-1 is exact) at a fifth of the cost
struct halfSpace
{
    vector n;
    point x0;
    label axis = -1;
    scalar sgn = 1;

    inline scalar dist(const point& p) const
    {
        return (axis < 0) ? ((p - x0) & n) : sgn*(p[axis] - x0[axis]);
    }
};

//- The axis of an axis-aligned half-space (none otherwise)
inline void setAxis(halfSpace& hs)
{
    hs.axis = -1;
    if (hs.n.y() == 0 && hs.n.z() == 0 && (hs.n.x() == 1 || hs.n.x() == -1)) { hs.axis = 0; hs.sgn = hs.n.x(); }
    else if (hs.n.x() == 0 && hs.n.z() == 0 && (hs.n.y() == 1 || hs.n.y() == -1)) { hs.axis = 1; hs.sgn = hs.n.y(); }
    else if (hs.n.x() == 0 && hs.n.y() == 0 && (hs.n.z() == 1 || hs.n.z() == -1)) { hs.axis = 2; hs.sgn = hs.n.z(); }
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
        s[i] = hs.dist(poly.p[i]);
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

    //- sweptExactCaps (2026-11-06): the piece took the exact path (no cap
    //  system), and, for a piece the closure solved, |V_exact - V_closure|
    //  (the validation of the exact path against the closure; -1 unset)
    bool exactUsed = false;
    scalar exactDiff = -1;

    //- The coverage: the total area of the boundary triangles' parts inside
    //  the cell planes (before the PLIC clip), for the check that the
    //  pieces of the ring cells cover the whole boundary
    scalar coverage = 0;

    //- The clipped triangles (partially cut by any plane)
    label nCut = 0;

    //- The half-space clips run (planeDistances + clipWith, or the clip of
    //  the pristine triangle) and the triangles entered (sweptClassStats,
    //  2026-11-05; counters only)
    label nClip = 0;
    label nTri = 0;

    //- The wetted area vector (in the list's orientation) and the coverage
    //  of the pre-image fan triangles (every fourth triangle, k % 4 == 1,
    //  of sweptTriangles' list) after the same clips: the face update from
    //  the volume clip (sweptFaceUpdateFromVolume, 2026-11-05)
    vector aPre = Zero;
    scalar covPre = 0;

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
//- The exact volume of the piece for ANY number of caps (the deformation
//  diagnosis of 2026-11-06, control (a), key sweptExactCaps of isoAdvection):
//  the closed triangle list of the swept volume is decomposed into the
//  signed tetrahedra (x_ref, T) about one point (the volume of a closed
//  oriented triangulation as the sum of signed tetrahedra, the divergence
//  theorem face by face: C. Zhang, T. Chen, Efficient feature extraction
//  for 2D/3D objects in mesh representation, Proceedings of the IEEE
//  International Conference on Image Processing (ICIP 2001), vol. 3,
//  935-938, doi:10.1109/ICIP.2001.958278); every tetrahedron is clipped as
//  a convex polyhedron by the cell's half-spaces and, where the cell holds
//  a plane, by the PLIC half-space, face by face with the polygon clipping
//  of I. E. Sutherland, G. W. Hodgman, Reentrant polygon clipping, Commun.
//  ACM 17 (1974) 32-42, doi:10.1145/360767.360802, the cap of each clip
//  closed from the crossing points ordered about their centroid in the
//  plane (a convex section); the volume of each clipped polyhedron is
//  (1/3) sum_faces (p_0 - x_ref) . a_face. The identity
//  sum_T vol(tet_T ∩ K) = vol(P ∩ K) holds for every convex K because the
//  signed indicators of the tetrahedra about one point sum to the
//  indicator of the closed volume P (the winding number), so no cap system
//  is solved: k > 3 caps and a PLIC plane parallel to a cell face need no
//  fallback. The same sign convention as sweptPieceVolume (the list's
//  orientation). Cost: one convex clip per triangle and plane; used by
//  sweptPieceVolume only where its cap system is not solvable.
inline scalar sweptPieceExactVolume
(
    const UList<FixedList<point, 3>>& T,
    const UList<halfSpace>& cellPlanes,
    const bool wetPlane,
    const halfSpace& plic,
    const point& xRef
)
{
    struct poly
    {
        FixedList<point, 24> p;
        label n = 0;
    };
    const label nc = cellPlanes.size();
    scalar vol = 0;
    forAll(T, k)
    {
        const FixedList<point, 3>& t = T[k];
        const scalar L =
            max(mag(t[1] - t[0]), max(mag(t[2] - t[1]), mag(t[0] - t[2])));
        if (!(L > 0)) continue;
        const vector aT = 0.5*((t[1] - t[0]) ^ (t[2] - t[0]));
        if (mag(aT) <= 1e-16*L*L) continue;
        const scalar vs = ((t[0] - xRef) & aT)/3;
        if (mag(vs) <= 1e-30) continue;
        const scalar sgn = (vs > 0) ? 1 : -1;
        // the base (A, B, C) with the apex on its negative side, the four
        // faces outward
        point A = t[0], B = t[1], C = t[2];
        if (sgn < 0) Swap(B, C);
        FixedList<poly, 16> P;
        label nP = 4;
        auto set3 = [](poly& q, const point& u, const point& v, const point& w)
        {
            q.p[0] = u; q.p[1] = v; q.p[2] = w; q.n = 3;
        };
        set3(P[0], A, B, C);
        set3(P[1], xRef, B, A);
        set3(P[2], xRef, C, B);
        set3(P[3], xRef, A, C);
        bool anyLeft = true;
        auto clipAll = [&](const halfSpace& hs)
        {
            FixedList<point, 48> cross;
            label nCross = 0;
            anyLeft = false;
            for (label i = 0; i < nP; ++i)
            {
                poly& q = P[i];
                if (q.n < 3) continue;
                poly out;
                for (label j = 0; j < q.n; ++j)
                {
                    const point& p0 = q.p[j];
                    const point& p1 = q.p[(j + 1) % q.n];
                    const scalar s0 = hs.dist(p0), s1 = hs.dist(p1);
                    if (s0 <= 0 && out.n < 24) out.p[out.n++] = p0;
                    if ((s0 <= 0) != (s1 <= 0))
                    {
                        const point x = p0 + (p1 - p0)*(s0/(s0 - s1));
                        if (out.n < 24) out.p[out.n++] = x;
                        if (nCross < 48) cross[nCross++] = x;
                    }
                }
                if (out.n < 3) out.n = 0;
                q = out;
                if (q.n >= 3) anyLeft = true;
            }
            if (nCross >= 3 && nP < 16)
            {
                // the plane's outward normal (the kept side is dist <= 0)
                vector nrm = hs.n;
                if (hs.axis >= 0)
                {
                    nrm = Zero;
                    nrm[hs.axis] = hs.sgn;
                }
                nrm /= max(mag(nrm), VSMALL);
                point cen(Zero);
                for (label i = 0; i < nCross; ++i) cen += cross[i];
                cen /= nCross;
                // a basis in the plane
                vector u = (mag(nrm.x()) < 0.9) ? vector(1, 0, 0) : vector(0, 1, 0);
                u -= (u & nrm)*nrm;
                u /= max(mag(u), VSMALL);
                const vector v = nrm ^ u;
                FixedList<scalar, 48> ang;
                FixedList<label, 48> idx;
                for (label i = 0; i < nCross; ++i)
                {
                    const vector d = cross[i] - cen;
                    ang[i] = Foam::atan2(d & v, d & u);
                    idx[i] = i;
                }
                // insertion sort by the angle
                for (label i = 1; i < nCross; ++i)
                {
                    const label key = idx[i];
                    label j = i - 1;
                    while (j >= 0 && ang[idx[j]] > ang[key])
                    {
                        idx[j + 1] = idx[j];
                        --j;
                    }
                    idx[j + 1] = key;
                }
                poly& cap = P[nP];
                cap.n = 0;
                for (label i = 0; i < nCross && cap.n < 24; ++i)
                {
                    cap.p[cap.n++] = cross[idx[i]];
                }
                // orient the cap outward
                vector ac(Zero);
                for (label i = 0; i < cap.n; ++i)
                {
                    ac += cap.p[i] ^ cap.p[(i + 1) % cap.n];
                }
                if ((ac & nrm) < 0)
                {
                    for (label i = 0; i < cap.n/2; ++i)
                    {
                        Swap(cap.p[i], cap.p[cap.n - 1 - i]);
                    }
                }
                ++nP;
                anyLeft = true;
            }
        };
        for (label i = 0; i < nc && anyLeft; ++i)
        {
            clipAll(cellPlanes[i]);
        }
        if (wetPlane && anyLeft)
        {
            clipAll(plic);
        }
        if (!anyLeft) continue;
        scalar v3 = 0;
        for (label i = 0; i < nP; ++i)
        {
            const poly& q = P[i];
            if (q.n < 3) continue;
            vector a(Zero);
            for (label j = 0; j < q.n; ++j)
            {
                a += (q.p[j] - xRef) ^ (q.p[(j + 1) % q.n] - xRef);
            }
            v3 += (q.p[0] - xRef) & (0.5*a);
        }
        vol += sgn*v3/3;
    }
    return vol;
}


inline pieceResult sweptPieceVolume
(
    const UList<FixedList<point, 3>>& T,
    const UList<halfSpace>& cellPlanes,
    const bool wetPlane,
    const halfSpace& plic,
    const point& xRef,
    const point* xSide = nullptr,
    const bool crossing = false,
    const bool exactCaps = false
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
    // the cap planes: cut[i] for the cell plane i, plicCut for the PLIC
    // plane (id nc), from the edges of the final wet polygons (the file
    // header, rule b)
    FixedList<bool, 8> cut(false);
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
            ++r.nClip;
            planeDistances(*cur, hs, s, anyIn, anyOut, sMaxAbs);
            if (clipWith(*cur, s, anyIn, anyOut, *nxt, id)) { Swap(cur, nxt); return true; }
            return false;
        };
        ++r.nTri;

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
            if (k % 4 == 1) r.covPre += mag(cur->areaVector());
            if (wetPlane)
            {
                if (clip(plic, -1)) { plicCut = true; wasCut = true; }
                if (cur->n < 3) continue;
            }
            if (wasCut) ++r.nCut;
            const vector a = cur->areaVector();
            if (k % 4 == 1) r.aPre += a;
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
        FixedList<bool, 8> inPlane(false);
        FixedList<FixedList<scalar, 3>, 8> s3;
        FixedList<scalar, 8> s3min, s3max;
        bool contact = false;
        bool beyond = false;
        for (label i = 0; i < nc; ++i)
        {
            scalar sMax = 0;
            s3min[i] = GREAT; s3max[i] = -GREAT;
            for (label j = 0; j < 3; ++j)
            {
                s3[i][j] = cellPlanes[i].dist(t[j]);
                sMax = max(sMax, mag(s3[i][j]));
                s3min[i] = min(s3min[i], s3[i][j]); s3max[i] = max(s3max[i], s3[i][j]);
            }
            if (sMax <= 1e-10*L)
            {
                const scalar sSide = xSide ? cellPlanes[i].dist(*xSide) : 0;
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
                ++r.nClip;
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
        if (k % 4 == 1) r.covPre += mag(ac);
        if (wetPlane)
        {
            if (clip(plic, nc)) wasCut = true;
            if (cur->n < 3) continue;
        }
        const clipPolygon& poly = *cur;
        const vector a = poly.areaVector();
        if (mag(a) <= areaTol) continue;
        if (k % 4 == 1) r.aPre += a;

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
        if (exactCaps)
        {
            // control (a): the exact piece volume, no cap system
            r.volume = sweptPieceExactVolume(T, cellPlanes, wetPlane, plic, xRef);
            r.closureResidual = 0;
            r.exactUsed = true;
            return r;
        }
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
            if (exactCaps)
            {
                r.volume = sweptPieceExactVolume(T, cellPlanes, wetPlane, plic, xRef);
                r.closureResidual = 0;
                r.exactUsed = true;
                return r;
            }
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
    if (exactCaps)
    {
        // the validation of the exact path on the pieces the closure solves
        r.exactDiff =
            mag(sweptPieceExactVolume(T, cellPlanes, wetPlane, plic, xRef) - r.volume);
    }
    return r;
}


//- The piece of a cell that ONE plane alone cuts, in closed form (round 4,
//  gain 5, key sweptSliverClosedForm, 2026-11-05): the cell holds no PLIC
//  plane and every other cell plane has no vertex of the swept volume
//  outside, so the piece is Omega_f cut by hs. Each triangle's part inside
//  hs is a triangle, a quadrilateral or nothing (davof::triNegativePolygon,
//  the one-donor clip's closed form, no polygon clipping), the volume is
//  (1/3) sum_T (x_T - x0) . a_T about the plane's point x0, where the cap's
//  moment (x_cap - x0) . a_cap vanishes (the cap lies in the plane): no cap
//  system. The contact rule and the remnant filter of sweptPieceVolume are
//  kept (a triangle lying in the plane with xSide outside is dropped, with
//  xSide inside kept whole; a part below 1e-14 L^2 is dropped). Round-off
//  against the clip: the parts' area vectors come from the closed form
//  instead of the clipped polygon, and the moments are taken about x0
//  instead of the face centre with the cap solved. aPre and covPre as in
//  sweptPieceVolume (k % 4 == 1). No citation: the identity is the
//  divergence theorem on the piece's closed boundary.
inline pieceResult sweptPieceOnePlane
(
    const UList<FixedList<point, 3>>& T,
    const halfSpace& hs,
    const point* xSide = nullptr
)
{
    pieceResult r;
    scalar threeV = 0;
    forAll(T, k)
    {
        const FixedList<point, 3>& t = T[k];
        const scalar L =
            max(mag(t[1] - t[0]), max(mag(t[2] - t[1]), mag(t[0] - t[2])));
        if (!(L > 0)) continue;
        ++r.nTri;
        const scalar s0 = hs.dist(t[0]), s1 = hs.dist(t[1]), s2 = hs.dist(t[2]);
        const scalar sMax = max(mag(s0), max(mag(s1), mag(s2)));
        const vector aFull = 0.5*((t[1] - t[0]) ^ (t[2] - t[0]));
        vector a(Zero);
        if (sMax <= 1e-10*L)
        {
            // rule a: a triangle lying in the plane
            const scalar sSide = xSide ? hs.dist(*xSide) : 0;
            const bool outside = (mag(sSide) > 1e-10*L) ? (sSide > 0) : ((aFull & hs.n) < 0);
            if (outside) { ++r.nContact; continue; }
            a = aFull;
        }
        else if (s0 > 0 && s1 > 0 && s2 > 0)
        {
            continue;
        }
        else if (s0 <= 0 && s1 <= 0 && s2 <= 0)
        {
            a = aFull;
        }
        else
        {
            point xw;
            triNegativePolygon(t[0], t[1], t[2], s0, s1, s2, a, xw);
            ++r.nCut;
        }
        ++r.nClip;
        const scalar areaTol = 1e-14*L*L;
        if (mag(a) <= areaTol) continue;
        r.coverage += mag(a);
        if (k % 4 == 1) { r.covPre += mag(a); r.aPre += a; }
        // the moment about x0: t[0] lies in the part's plane
        threeV += (t[0] - hs.x0) & a;
    }
    r.volume = threeV/3;
    r.nCaps = 1;
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
    bool* coplanarHit = nullptr,
    label* nClipOut = nullptr
)
{
    vector a(Zero);
    coverage = 0;
    if (aCell) *aCell = Zero;
    if (coplanarHit) *coplanarHit = false;
    label nClip = 0;
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
            ++nClip;
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
            ++nClip;
            planeDistances(*cur, plic, s, anyIn, anyOut, sMaxAbs);
            if (clipWith(*cur, s, anyIn, anyOut, *nxt, -1)) Swap(cur, nxt);
            if (cur->n < 3) continue;
        }
        a += cur->areaVector();
    }
    if (nClipOut) *nClipOut += nClip;
    return a;
}

} // End namespace davof
} // End namespace Foam

#endif

// ************************************************************************* //
