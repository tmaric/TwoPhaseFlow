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
    The semi-Lagrangian swept-face flux of the DAVOF dual advection
    (fluxIntegration sweptFace and sweptFaceCells, 2026-10-28/29; the working
    note leia/docs/davof/davof-swept-volume, whose Python prototypes
    swept_volume.py and multicell.py these routines follow).

    The engineering problem: the prism flux (fluxIntegration prism) assumes
    that the face translates by the displacement of its centre over the step,
    so under a non-uniform velocity its flux volume is not the volume the face
    sweeps. Here the vertices p_i and the centre x_f of the face f at t^{n+1}
    are traced BACKWARD along the velocity to the pre-image f' at t^n (p_i',
    x_f'), as in the face-vertex flux computation of T. Maric, H. Marschall,
    D. Bothe, An enhanced un-split face-vertex flux-based VoF method,
    J. Comput. Phys. 371 (2018) 967-993, doi:10.1016/j.jcp.2018.03.048. The
    swept volume is the closed list of 4F planar triangles
        (x_f, p_i, p_{i+1})         the face fan,
        (x_f', p_{i+1}', p_i')      the pre-image fan, reversed,
        (p_i, p_i', p_{i+1}'), (p_i, p_{i+1}', p_{i+1})   the lateral quad,
    outward for phi_f = U_f.S_f > 0 with S_f owner -> neighbour; for
    phi_f < 0 every triangle is reversed. The signed volume
    (1/6) sum_T p.(q x r) checks the orientation and the list is flipped where
    it is negative (counted).

    sweptFace (ONE DONOR CELL): every triangle is clipped against the donor
    cell's plane (liquid where (x - x0).n <= 0; the wet polygon's area vector
    a_T^- by davof::triNegativePolygon), and

        V_f^-          = (1/3) sum_T (p_T - x0).a_T^-     (p_T a vertex of T)
        alpha_f^{n+1}  = (a_{f'}^- . A_{f'}) / |A_{f'}|^2

    with a_{f'}^- the sum of the wet area vectors of the pre-image fan and
    A_{f'} the pre-image's area vector: the flux and the new face fraction
    come from the same clips of the same donor plane. The closed form is the
    divergence theorem for x on the truncated volume, where the cap integral
    of x.n is (x0.n) times the cap area and the cap area follows from the
    closure of the wet boundary (the note's eqs. 3-5); with x0 on the plane
    the cap term vanishes identically (the hydrostatic form). The faces whose
    pre-image vertices (or centre) leave the donor cell are counted per step
    ("out of donor"; the sweeping faces only since 4c7f3b5).

    sweptFaceCells (THE PIECEWISE CLIP, 2026-10-29): the swept volume is cut
    into the pieces inside the donor and inside every cell sharing a vertex
    with it (the vertex ring; one layer of neighbours, enough at Co 0.2), by
    clipping every boundary triangle against each cell's six face planes
    (davofSweptPieces.H). A piece in a FULL cell (alpha^n >= 1 - 1e-8) is
    wet in total, in an EMPTY cell (alpha^n <= 1e-8) dry, in an INTERFACE
    cell (a cell with a surface) clipped by THAT cell's PLIC plane of t^n,
    the donor's pieces by the donor's plane; a partial cell without a
    surface (a wisp) is taken wet where alpha^n >= 1/2 and dry otherwise
    (counted). The volume of a piece follows from the clipped triangles and
    the cap areas solved from closure (k <= 3 caps: the cut cell planes and
    the PLIC plane), V = (1/3)[sum_T p_T.a_T^- + sum_i A_i (x0_i.n_i)]. The
    flux volume is the sum over the pieces; the face fraction is the same
    clip of the pre-image fan alone (2D, no caps). Logged per step: the
    pieces, the faces whose boundary is not covered by the ring (a part
    beyond one layer: that part is skipped, nothing is invented, counted
    with the largest uncovered share), the rank-deficient pieces (more than
    three caps, or a PLIC plane parallel to a cell face: the face falls back
    to the one-donor clip, counted), the largest closure residual of the
    pieces, the partial cells without a surface.

    dVf carries the sign of phi_f, as the prism flux does, and both cells of
    a face see the same volume (conservation by antisymmetry). The existing
    bounding and clipping of alpha are untouched.

    SCOPE AND LIMITS (stated here and in etc/README-davof.md):

    - PLANES ONLY. The formula needs x.n constant on the cap. Every surface
      used must be a plane: a cellSurface without the paraboloid flag, or
      one whose curvature tensor is at round-off, max|b_ij| h <= sweptPlaneTol
      (default 1e-6; on the translated plane the ring tensor is 1e-10 /m at
      h = 2e-4 m, b h = 2e-14). A true paraboloid or a sphere surface is a
      FatalError, never a silent fallback: use clippingSurface plane, or
      fluxIntegration prism with the paraboloid.

    - TRAJECTORIES. Backward midpoint rule (explicit RK2, second order in
      dt): x' = x - dt U(x - dt/2 U(x)). The vertex velocity U(p_i) is the
      volPointInterpolation value at the vertex, the face centre's is the
      cell-point interpolate of U at x_f in the donor (the prism's
      displacement d = U(x_f) dt, so for a uniform U the pre-image is the
      prism's f - d and the swept volume IS the prism). The midpoint velocity
      is the cell-point interpolate in the cell that holds the midpoint (the
      donor, a cell sharing a vertex with it, local or on another rank, else
      the local mesh's octree search); a midpoint outside the mesh, or in a
      cell of another rank beyond the donor's vertex ring, takes the Euler
      step x - dt U(x) (exact for a uniform U as well; counted, "Euler
      traces").

    - CELLS AND PARALLEL RUNS (2026-10-30). The ring is the donor's
      cell-point-cell stencil. In parallel its cells on other ranks are
      fetched once per step (sweptExchange: OpenFOAM's zoneDistribute
      stencil with global cell indices, H. Scheufler, J. Roenby, J. Comput.
      Phys. 383 (2019) 1-23, doi:10.1016/j.jcp.2019.01.009, and one
      mapDistributeBase transfer of the packed cells, davofSweptRemote.H):
      alpha^n, the surface, the face planes, the bounding box and the
      tetrahedra of the cell-point interpolation with the point velocities,
      so that a piece or a midpoint in a remote cell is computed from the
      same numbers as in serial. The remote cells follow the local ones in
      the ring (the sums differ from serial in their order only, round-off).
      A processor face is computed on its upwind side (the donor's rank).
      In serial nothing is exchanged and the path is the local one
      unchanged. Not covered: cells beyond one layer (as in serial),
      cyclic patches across ranks, the sweptClassify and sweptFaceDump
      diagnostics (local cells only).

    - DIAGNOSTICS per step in the Info line: the faces computed, the sign
      flips, the faces out of the donor, the largest closure residual of the
      lists |sum_T a_T| / sum_T |a_T| (round-off for a closed list by
      construction), the number of clipped triangles; for sweptFaceCells the
      piece counts above.

SourceFiles
    isoAdvectionSwept.C

\*---------------------------------------------------------------------------*/

#include "isoAdvection.H"
#include "davofSimplexGeometry.H"
#include "davofSweptPieces.H"
#include "interpolationCellPoint.H"
#include "pointFields.H"
#include "FixedList.H"
#include "boundBox.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "zoneDistribute.H"
#include "mapDistributeBase.H"
#include "polyMeshTetDecomposition.H"
#include "tetIndices.H"
#include <chrono>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

//- The point-in-cell test by the face planes (convex cells), with a
//  tolerance in units of the cell size h
inline bool inCellTol
(
    const Foam::fvMesh& mesh, const Foam::point& x, const Foam::label c,
    const Foam::scalar tolRel, const Foam::scalar h
)
{
    const Foam::vectorField& Cf = mesh.faceCentres();
    const Foam::vectorField& Sf = mesh.faceAreas();
    const Foam::labelList& own = mesh.faceOwner();
    for (const Foam::label cf : mesh.cells()[c])
    {
        const Foam::vector nOut = (own[cf] == c) ? Sf[cf] : -Sf[cf];
        if (((x - Cf[cf]) & nOut) > tolRel*h*Foam::mag(nOut))
        {
            return false;
        }
    }
    return true;
}


//- The plane of a cell surface (a plane, or a paraboloid with a tensor at
//  round-off): the raw normal (out of the liquid) and a point on it; a
//  sphere or a true paraboloid is a FatalError
inline void planeOf
(
    const Foam::davof::cellSurface& S, const Foam::label celli,
    const Foam::scalar h, const Foam::scalar tol, const char* flux,
    Foam::vector& nRaw, Foam::point& x0
)
{
    using namespace Foam;
    if (S.sphere)
    {
        FatalErrorInFunction
            << "fluxIntegration " << flux << ": the surface of cell " << celli
            << " is a sphere (clippingSurface exactSphere). The swept-face"
            << " flux is implemented for a plane only (x.n is constant on the"
            << " cap); use clippingSurface plane or fluxIntegration prism."
            << exit(FatalError);
    }
    const scalar bMax = max(mag(S.b11), max(mag(S.b12), mag(S.b22)));
    if (S.paraboloid && bMax*h > tol)
    {
        FatalErrorInFunction
            << "fluxIntegration " << flux << ": the surface of cell " << celli
            << " is a paraboloid (curvature tensor b11 " << S.b11 << " b12 "
            << S.b12 << " b22 " << S.b22 << " /m, max|b| h = " << bMax*h
            << " above sweptPlaneTol " << tol << "). The"
            << " swept-face flux is implemented for a plane only (x.n is"
            << " constant on the cap); use clippingSurface plane, or"
            << " fluxIntegration prism with the paraboloid."
            << exit(FatalError);
    }
    nRaw = S.paraboloid ? (S.e3 - S.a1*S.e1 - S.a2*S.e2) : S.e3;
    x0 = S.o + S.c*S.e3;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::advection::isoAdvection::sweptTriangles
(
    const label facei,
    const label donor,
    const scalar dt,
    const scalar phiF,
    const bool sweeps,
    const interpolationCellPoint<vector>& UInterp,
    const pointVectorField& Up,
    List<FixedList<point, 3>>& T,
    pointField& q,
    point& xq,
    label& preSign,
    scalar& vTotal
)
{
    const face& f = mesh_.faces()[facei];
    const pointField& P = mesh_.points();
    const vectorField& Cf = mesh_.faceCentres();
    const labelListList& cellPoints = mesh_.cellPoints();
    const labelListList& pointCells = mesh_.pointCells();
    const label F = f.size();
    const scalar h = Foam::cbrt(mesh_.V()[donor]);

    // The cell that holds a point of the trajectory: the donor, else the
    // cells sharing a vertex with the donor (a midpoint lies within Co h/2
    // of the donor's boundary, inside that ring for Co < 1), else the mesh's
    // octree search (CELL_TETS; the tree and the tet base points are built
    // once per step in curvedFlux, on every rank; the face-plane search of
    // polyMesh::findCell scans every cell and is not used); in parallel the
    // cells of other ranks in the donor's vertex ring (sweptExchange) are
    // searched after the local ring and before the octree, and returned as
    // -2 - k (k their index in sweptRemote_); -1 outside the mesh, or
    // beyond the ring on another rank
    const auto remoteIter = sweptRemoteRing_.cfind(donor);
    const labelList* remoteRing = remoteIter.good() ? &remoteIter.val() : nullptr;
    auto cellOf = [&](const point& x) -> label
    {
        if (inCellTol(mesh_, x, donor, 1e-12, h)) return donor;
        for (const label pi : cellPoints[donor])
        {
            for (const label c : pointCells[pi])
            {
                if (c != donor && inCellTol(mesh_, x, c, 1e-12, h)) return c;
            }
        }
        if (remoteRing)
        {
            for (const label k : *remoteRing)
            {
                if (davof::inRemoteCellTol(sweptRemote_[k], x, 1e-12, h)) return -2 - k;
            }
        }
        return mesh_.findCell(x, polyMesh::CELL_TETS);
    };

    // The velocity at a point of the trajectory: the cell-point interpolate
    // in the cell that holds the point (a remote cell's from its packed
    // tetrahedra, davof::interpolateRemote); nothing where no cell holds it
    auto velocityNear = [&](const point& x, bool& found) -> vector
    {
        const label c = cellOf(x);
        found = (c >= 0 || c <= -2);
        if (c >= 0) return UInterp.interpolate(x, c);
        if (c <= -2)
        {
            sweptFaceRemote_ = true;
            return davof::interpolateRemote(sweptRemote_[-2 - c], x);
        }
        return vector(Zero);
    };

    // The backward midpoint trace from x with the velocity u1 at x
    auto traceBack = [&](const point& x, const vector& u1) -> point
    {
        const point xm = x - 0.5*dt*u1;
        bool found = false;
        const vector um = velocityNear(xm, found);
        if (!found) ++nSweptEulerTrace_;
        return found ? (x - dt*um) : (x - dt*u1);
    };

    // The pre-image: the vertices by the point velocities, the centre by
    // the cell-point interpolate in the donor (the prism's displacement)
    const point xf = Cf[facei];
    q.setSize(F);
    forAll(f, i)
    {
        q[i] = traceBack(P[f[i]], Up[f[i]]);
    }
    xq = traceBack(xf, UInterp.interpolate(xf, donor));

    // The out-of-donor count: any traced point outside the donor, of the
    // faces that sweep a volume (a face parallel to the velocity carries no
    // flux, and its traced vertices leave the donor along its side faces by
    // construction; the studies of 2026-10-28 counted those faces too)
    if (sweeps)
    {
        bool outOfDonor = !inCellTol(mesh_, xq, donor, 1e-9, h);
        forAll(q, i)
        {
            if (!inCellTol(mesh_, q[i], donor, 1e-9, h)) outOfDonor = true;
        }
        if (outOfDonor) ++nSweptOutOfDonor_;
    }

    // The triangle list, outward for phi_f > 0 (the face fan along S_f, the
    // pre-image fan reversed, the lateral quads split along one diagonal)
    T.setSize(4*F);
    forAll(f, i)
    {
        const point& a = P[f[i]];
        const point& b = P[f[f.fcIndex(i)]];
        const point& aq = q[i];
        const point& bq = q[f.fcIndex(i)];
        T[4*i][0] = xf; T[4*i][1] = a;  T[4*i][2] = b;
        T[4*i + 1][0] = xq; T[4*i + 1][1] = bq; T[4*i + 1][2] = aq;
        T[4*i + 2][0] = a;  T[4*i + 2][1] = aq; T[4*i + 2][2] = bq;
        T[4*i + 3][0] = a;  T[4*i + 3][1] = bq; T[4*i + 3][2] = b;
    }
    // the orientation of the stored pre-image triangles against the face's
    // own order (x_f', p_i', p_{i+1}'): -1 as built
    preSign = -1;
    auto flipAll = [&]()
    {
        forAll(T, k) Swap(T[k][1], T[k][2]);
        preSign = -preSign;
    };
    if (phiF < 0) flipAll();

    vTotal = 0;
    if (sweeps)
    {
        // the signed volume about x_f (translation invariant; the local
        // origin avoids the cancellation of the global coordinates)
        scalar sixV = 0;
        forAll(T, k)
        {
            sixV += (T[k][0] - xf) & ((T[k][1] - xf) ^ (T[k][2] - xf));
        }
        if (sixV < 0)
        {
            flipAll();
            ++nSweptFlips_;
            sixV = -sixV;
        }
        vTotal = sixV/6;
        // the closure residual of the list
        vector sumA(Zero);
        scalar sumMag = 0;
        forAll(T, k)
        {
            const vector a = 0.5*((T[k][1] - T[k][0]) ^ (T[k][2] - T[k][0]));
            sumA += a;
            sumMag += mag(a);
        }
        if (sumMag > VSMALL)
        {
            sweptClosureMax_ = max(sweptClosureMax_, mag(sumA)/sumMag);
        }
    }
    ++nSweptFaces_;
}


void Foam::advection::isoAdvection::sweptFaceFlux
(
    const label facei,
    const label donor,
    const davof::cellSurface& S,
    const scalar dt,
    const scalar phiF,
    const bool sweeps,
    const interpolationCellPoint<vector>& UInterp,
    const pointVectorField& Up,
    scalar& af,
    vector& aWet,
    scalar& vLiq,
    pointField* qOut,
    point* xqOut
)
{
    const scalar h = Foam::cbrt(mesh_.V()[donor]);
    af = 0;
    aWet = Zero;
    vLiq = 0;

    // The donor's plane: (x - x0).nRaw = 0, liquid where the value is
    // negative (the surface's own sign convention: e3 points out of the
    // liquid). A paraboloid with the slopes a1, a2 and a zero tensor is the
    // plane with the normal e3 - a1 e1 - a2 e2 through o + c e3.
    vector nRaw;
    point x0;
    planeOf(S, donor, h, sweptPlaneTol_, "sweptFace", nRaw, x0);
    auto sdist = [&](const point& x) -> scalar
    {
        return (x - x0) & nRaw;
    };

    List<FixedList<point, 3>> T;
    pointField q;
    point xq;
    label preSign;
    scalar vTotal;
    sweptTriangles(facei, donor, dt, phiF, sweeps, UInterp, Up, T, q, xq, preSign, vTotal);
    if (qOut) *qOut = q;
    if (xqOut) *xqOut = xq;

    // The clips: the wet area vector of every triangle, the volume by the
    // closed form with the apex on the plane, the pre-image's wet and full
    // area vectors from its fan (k % 4 == 1)
    scalar threeV = 0;
    vector aPreWet(Zero), aPreFull(Zero);
    forAll(T, k)
    {
        const FixedList<point, 3>& t = T[k];
        const scalar s0 = sdist(t[0]), s1 = sdist(t[1]), s2 = sdist(t[2]);
        const vector aFull = 0.5*((t[1] - t[0]) ^ (t[2] - t[0]));
        vector aw(Zero);
        if (s0 <= 0 && s1 <= 0 && s2 <= 0)
        {
            aw = aFull;
        }
        else if (s0 > 0 && s1 > 0 && s2 > 0)
        {
            aw = Zero;
        }
        else
        {
            point xw;
            davof::triNegativePolygon(t[0], t[1], t[2], s0, s1, s2, aw, xw);
            ++nSweptClipTris_;
        }
        threeV += (t[0] - x0) & aw;
        if (k % 4 == 1)
        {
            aPreWet += aw;
            aPreFull += aFull;
        }
    }
    vLiq = sweeps ? threeV/3 : 0;
    // the kinematic consistency: the liquid share of the swept volume
    // times the flux phi_f dt (sweptScaleToPhi)
    if (sweptScaleToPhi_ && sweeps && vTotal > VSMALL)
    {
        vLiq *= mag(phiF)*dt/vTotal;
    }

    // The face fraction of the pre-image in the face's orientation
    aWet = preSign*aPreWet;
    const vector Apre = preSign*aPreFull;
    const scalar A2 = magSqr(Apre);
    af = (A2 > VSMALL) ? (aWet & Apre)/A2 : 0;
}


void Foam::advection::isoAdvection::sweptFaceCellsFlux
(
    const label facei,
    const label donor,
    const davof::cellSurface& S,
    const scalar dt,
    const scalar phiF,
    const bool sweeps,
    const interpolationCellPoint<vector>& UInterp,
    const pointVectorField& Up,
    scalar& af,
    vector& aWet,
    scalar& vLiq,
    pointField* qOut,
    point* xqOut
)
{
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();
    const labelList& own = mesh_.faceOwner();
    const cellList& cells = mesh_.cells();
    const labelListList& cellPoints = mesh_.cellPoints();
    const labelListList& pointCells = mesh_.pointCells();
    const scalarField& a0 = alpha1_.oldTime().primitiveField();
    const scalar h = Foam::cbrt(mesh_.V()[donor]);
    af = 0;
    aWet = Zero;
    vLiq = 0;

    List<FixedList<point, 3>> T;
    pointField q;
    point xq;
    label preSign;
    scalar vTotal;
    sweptTriangles(facei, donor, dt, phiF, sweeps, UInterp, Up, T, q, xq, preSign, vTotal);
    if (qOut) *qOut = q;
    if (xqOut) *xqOut = xq;

    // the pre-image fan in the face's orientation (x_f', p_i', p_{i+1}')
    const label F = q.size();
    List<FixedList<point, 3>> Tpre(F);
    forAll(q, i)
    {
        Tpre[i][0] = xq; Tpre[i][1] = q[i]; Tpre[i][2] = q[(i + 1) % F];
    }
    vector aPreFull(Zero);
    forAll(Tpre, i)
    {
        aPreFull += 0.5*((Tpre[i][1] - Tpre[i][0]) ^ (Tpre[i][2] - Tpre[i][0]));
    }
    scalar sumMagT = 0;
    forAll(T, k)
    {
        sumMagT += mag(0.5*((T[k][1] - T[k][0]) ^ (T[k][2] - T[k][0])));
    }
    scalar sumMagPre = 0;
    forAll(Tpre, k)
    {
        sumMagPre += mag(0.5*((Tpre[k][1] - Tpre[k][0]) ^ (Tpre[k][2] - Tpre[k][0])));
    }

    // The ring: the donor and the cells sharing a vertex with it whose
    // bounding box meets the swept volume's (a cell that does not touch
    // the swept volume holds no piece; at Co 0.2 the swept volume of a face
    // touches the donor and the cells about the face's edges and vertices)
    boundBox sweptBox(T[0][0], T[0][0]);
    forAll(T, k)
    {
        sweptBox.add(T[k][0]); sweptBox.add(T[k][1]); sweptBox.add(T[k][2]);
    }
    sweptBox.inflate(1e-9);
    const pointField& P = mesh_.points();
    auto touches = [&](const label c) -> bool
    {
        boundBox cellBox(P, cellPoints[c], false);
        return cellBox.overlaps(sweptBox);
    };
    DynamicList<label> ring(32);
    ring.append(donor);
    for (const label pi : cellPoints[donor])
    {
        for (const label c : pointCells[pi])
        {
            if (c != donor && !ring.found(c) && touches(c)) ring.append(c);
        }
    }
    // in parallel the ring cells of other ranks (sweptExchange), after the
    // local ones, as -2 - k (k their index in sweptRemote_)
    {
        const auto remoteIter = sweptRemoteRing_.cfind(donor);
        if (remoteIter.good())
        {
            for (const label k : remoteIter.val())
            {
                const davof::remoteCell& rc = sweptRemote_[k];
                if (boundBox(rc.bbMin, rc.bbMax).overlaps(sweptBox)) ring.append(-2 - k);
            }
        }
    }
    const point xf = Cf[facei];

    // The coplanar rule of the face update (sweptCoplanar, davofSweptPieces.H):
    // donorSide gives a part of the pre-image lying in a cell plane to the
    // cell on the donor's side of that plane, nRef the donor's outward
    // normal through the face; both keeps it in both cells (the defect)
    const vector SfOut = (own[facei] == donor) ? Sf[facei] : -Sf[facei];
    const vector nRef =
        (sweptCoplanar_ == "donorSide") ? SfOut/max(mag(SfOut), VSMALL) : vector(Zero);

    // the one-step face dump (sweptFaceDump): the pieces of the face update
    const bool dump = sweptFaceDump_ && mesh_.time().timeIndex() == 1;
    struct pieceRow { label cell; word cls; vector aw; vector ac; bool cop; };
    DynamicList<pieceRow> pieceRows(dump ? 16 : 0);

    // The pieces
    scalar threeV = 0;
    vector aPreWet(Zero);
    scalar coverage = 0, coveragePre = 0;
    bool rankFallback = false;
    List<davof::halfSpace> planes;
    for (const label e : ring)
    {
        // a ring entry: a local cell c = e, or a cell of another rank
        // (e <= -2) whose packed planes and state are read instead
        const davof::remoteCell* rc = (e <= -2) ? &sweptRemote_[-2 - e] : nullptr;
        const label c = rc ? rc->gid : e;
        // the cell's half-spaces
        if (rc)
        {
            davof::remoteCellPlanes(*rc, planes);
        }
        else
        {
            const labelList& cf = cells[c];
            planes.setSize(cf.size());
            forAll(cf, j)
            {
                const vector nOut = (own[cf[j]] == c) ? Sf[cf[j]] : -Sf[cf[j]];
                planes[j].n = nOut/max(mag(nOut), VSMALL);
                planes[j].x0 = Cf[cf[j]];
            }
        }
        // the cell's state and plane
        const scalar a0c = rc ? rc->alpha0 : a0[c];
        const bool hasS = rc ? rc->hasSurface : hasSurface_.test(c);
        bool wetPlane = false, wet = true;
        davof::halfSpace plic;
        plic.n = Zero; plic.x0 = Zero;
        if (hasS)
        {
            vector nRaw; point x0;
            planeOf(rc ? rc->S : cellSurfaces_[c], c, h, sweptPlaneTol_, "sweptFaceCells", nRaw, x0);
            plic.n = nRaw/max(mag(nRaw), VSMALL);
            plic.x0 = x0;
            wetPlane = true;
        }
        else if (a0c >= 1 - 1e-8)
        {
            wet = true;
        }
        else if (a0c <= 1e-8)
        {
            wet = false;
        }
        else
        {
            // a partial cell without a surface (a wisp): wet by its alpha
            ++nSweptPartialCells_;
            wet = (a0c >= 0.5);
        }

        if (sweeps)
        {
            const davof::pieceResult r = davof::sweptPieceVolume(T, planes, wetPlane, plic, xf);
            coverage += r.coverage;
            if (r.coverage > 1e-12*sumMagT)
            {
                if (rc) sweptFaceRemote_ = true;
                ++nSweptPieces_;
                if (!r.ok) rankFallback = true;
                if (r.degenerate) ++nSweptDegenerate_;
                if ((wetPlane || wet) && !r.degenerate)
                {
                    threeV += 3*r.volume;
                    sweptPieceClosureMax_ = max(sweptPieceClosureMax_, r.closureResidual);
                }
                nSweptClipTris_ += r.nCut;
            }
        }
        scalar cov = 0;
        vector ac(Zero);
        bool cop = false;
        const vector aw = davof::sweptPieceArea(Tpre, planes, wetPlane, plic, cov, nRef, &ac, &cop);
        coveragePre += cov;
        if (rc && cov > 0) sweptFaceRemote_ = true;
        if (cop) ++nSweptCoplanar_;
        if (wetPlane || wet) aPreWet += aw;
        if (dump && cov > 0)
        {
            word cls = (!rc && c == donor) ? "donor" : (wetPlane ? "interface" : (a0c >= 1 - 1e-8 ? "full" : (a0c <= 1e-8 ? "empty" : "wisp")));
            pieceRows.append({c, cls, (wetPlane || wet) ? aw : vector(Zero), ac, cop});
        }
    }

    // the parts of the boundary beyond the ring: skipped, counted
    if (sumMagT > VSMALL)
    {
        const scalar beyond = max(scalar(0), (sumMagT - coverage)/sumMagT);
        const scalar beyondPre = (sumMagPre > VSMALL) ? max(scalar(0), (sumMagPre - coveragePre)/sumMagPre) : 0;
        const scalar b = sweeps ? max(beyond, beyondPre) : beyondPre;
        if (b > 1e-9)
        {
            ++nSweptBeyond_;
            sweptBeyondMax_ = max(sweptBeyondMax_, b);
        }
    }

    // the face update of sweptFace: the pre-image fan clipped by the
    // donor's plane alone (sweptCellsFaceUpdate donor; computed for the dump
    // as well)
    auto donorWet = [&]() -> vector
    {
        vector nRaw; point x0;
        planeOf(S, donor, h, sweptPlaneTol_, "sweptFaceCells", nRaw, x0);
        vector w(Zero);
        forAll(Tpre, i)
        {
            const FixedList<point, 3>& t = Tpre[i];
            const scalar s0 = (t[0] - x0) & nRaw, s1 = (t[1] - x0) & nRaw, s2 = (t[2] - x0) & nRaw;
            if (s0 <= 0 && s1 <= 0 && s2 <= 0)
            {
                w += 0.5*((t[1] - t[0]) ^ (t[2] - t[0]));
            }
            else if (!(s0 > 0 && s1 > 0 && s2 > 0))
            {
                vector aw; point xw;
                davof::triNegativePolygon(t[0], t[1], t[2], s0, s1, s2, aw, xw);
                w += aw;
            }
        }
        return w;
    };
    const scalar A2 = magSqr(aPreFull);
    auto fractionOf = [&](const vector& w) -> scalar
    {
        return (A2 > VSMALL) ? (w & aPreFull)/A2 : 0;
    };

    if (rankFallback)
    {
        // the one-donor clip of this face (counted)
        ++nSweptRankFallback_;
        sweptFaceFlux(facei, donor, S, dt, phiF, sweeps, UInterp, Up, af, aWet, vLiq);
        if (dump)
        {
            writeFaceDump(facei, donor, phiF, sweeps, fractionOf(aPreWet), af, af, true,
                aPreFull, pieceRows.size());
        }
        return;
    }

    vLiq = sweeps ? threeV/3 : 0;
    // the kinematic consistency: the liquid share of the swept volume
    // times the flux phi_f dt (sweptScaleToPhi)
    if (sweptScaleToPhi_ && sweeps && vTotal > VSMALL)
    {
        vLiq *= mag(phiF)*dt/vTotal;
    }
    const scalar afPieces = fractionOf(aPreWet);
    const vector wDonor = (dump || sweptCellsFaceUpdate_ == "donor") ? donorWet() : vector(Zero);
    if (sweptCellsFaceUpdate_ == "donor")
    {
        aPreWet = wDonor;
    }
    aWet = aPreWet;
    af = fractionOf(aWet);

    if (dump)
    {
        writeFaceDump(facei, donor, phiF, sweeps, afPieces, fractionOf(wDonor), af, false,
            aPreFull, pieceRows.size());
        for (const pieceRow& r : pieceRows)
        {
            dumpPieces_() << facei << ',' << r.cell << ',' << r.cls << ','
                << r.aw.x() << ',' << r.aw.y() << ',' << r.aw.z() << ','
                << r.ac.x() << ',' << r.ac.y() << ',' << r.ac.z() << ','
                << (r.cop ? 1 : 0) << nl;
        }
    }
}


void Foam::advection::isoAdvection::writeFaceDump
(
    const label facei,
    const label donor,
    const scalar phiF,
    const bool sweeps,
    const scalar afPieces,
    const scalar afDonor,
    const scalar afUsed,
    const bool rankFallback,
    const vector& Apre,
    const label nPieces
)
{
    if (!dumpFaces_)
    {
        const fileName dir(mesh_.time().globalPath()/"postProcessing"/"sweptFaceDump");
        mkDir(dir);
        dumpFaces_.reset(new OFstream(dir/"faces.csv"));
        dumpFaces_().precision(16);
        dumpFaces_() << "facei,owner,neighbour,donor,phi,sweeps,afPieces,afDonor,afUsed,"
            << "rankFallback,nPieces,Apre_x,Apre_y,Apre_z,Sf_x,Sf_y,Sf_z,Cf_x,Cf_y,Cf_z" << nl;
        dumpPieces_.reset(new OFstream(dir/"pieces.csv"));
        dumpPieces_().precision(16);
        dumpPieces_() << "facei,cell,class,aWet_x,aWet_y,aWet_z,aCell_x,aCell_y,aCell_z,coplanar" << nl;
    }
    const vector& S = mesh_.faceAreas()[facei];
    const point& C = mesh_.faceCentres()[facei];
    const label nb = mesh_.isInternalFace(facei) ? mesh_.faceNeighbour()[facei] : -1;
    dumpFaces_() << facei << ',' << mesh_.faceOwner()[facei] << ',' << nb << ',' << donor << ','
        << phiF << ',' << (sweeps ? 1 : 0) << ',' << afPieces << ',' << afDonor << ',' << afUsed << ','
        << (rankFallback ? 1 : 0) << ',' << nPieces << ','
        << Apre.x() << ',' << Apre.y() << ',' << Apre.z() << ','
        << S.x() << ',' << S.y() << ',' << S.z() << ','
        << C.x() << ',' << C.y() << ',' << C.z() << nl;
}


void Foam::advection::isoAdvection::packSweptCell
(
    const label celli,
    const pointVectorField& Up,
    DynamicList<scalar>& v
) const
{
    const pointField& P = mesh_.points();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();
    const labelList& own = mesh_.faceOwner();
    v.clear();
    v.append(1);
    davof::appendVector(mesh_.cellCentres()[celli], v);
    v.append(mesh_.V()[celli]);
    v.append(alpha1_.oldTime().primitiveField()[celli]);
    const bool has = hasSurface_.test(celli);
    v.append(has ? 1 : 0);
    davof::packCellSurface(has ? cellSurfaces_[celli] : davof::cellSurface(), v);
    const boundBox bb(P, mesh_.cellPoints()[celli], false);
    davof::appendVector(bb.min(), v);
    davof::appendVector(bb.max(), v);
    davof::appendVector(U_.primitiveField()[celli], v);
    const labelList& cFaces = mesh_.cells()[celli];
    v.append(cFaces.size());
    for (const label f : cFaces)
    {
        davof::appendVector((own[f] == celli) ? Sf[f] : -Sf[f], v);
        davof::appendVector(Cf[f], v);
    }
    // the tetrahedra of the cell-point interpolation (cellPointWeight's
    // order: polyMeshTetDecomposition::cellTetIndices)
    const List<tetIndices> tets(polyMeshTetDecomposition::cellTetIndices(mesh_, celli));
    v.append(tets.size());
    for (const tetIndices& tetIs : tets)
    {
        const triFace tri(tetIs.faceTriIs(mesh_));
        for (direction k = 0; k < 3; ++k) davof::appendVector(P[tri[k]], v);
        for (direction k = 0; k < 3; ++k) davof::appendVector(Up[tri[k]], v);
    }
}


void Foam::advection::isoAdvection::sweptExchange
(
    const boolList& donors,
    const pointVectorField& Up
)
{
    sweptRemote_.clear();
    sweptRemoteRing_.clear();
    nSweptRemoteCells_ = 0;
    sweptExchangeTime_ = 0;
    sweptExchangeParts_ = Zero;
    if (!Pstream::parRun()) return;
    auto now = []() { return std::chrono::steady_clock::now(); };
    auto secs = [](const auto& a, const auto& b)
    {
        return std::chrono::duration<scalar>(b - a).count();
    };
    const auto t0 = now();

    // the wait for the slowest rank (the work before the exchange is not
    // balanced: the interface ranks fit and anchor their surfaces), timed
    // apart from the exchange itself
    UPstream::barrier(UPstream::worldComm);
    const auto t1 = now();

    // The cell-point-cell stencil of the donors with global indices
    // (OpenFOAM's zoneDistribute / zoneCPCStencil, H. Scheufler, J. Roenby,
    // J. Comput. Phys. 383 (2019) 1-23, doi:10.1016/j.jcp.2019.01.009). The
    // stencil of a cell is computed once and cached by zoneCPCStencil, but
    // every update synchronises the coupled boundary points of the whole
    // mesh; it is therefore called only when a donor of some rank has no
    // stencil yet, and then for the donors and two layers of their face
    // neighbours, so that the interface moving by a fraction of a cell per
    // step finds its next donors' stencils in place
    zoneDistribute& zd = zoneDistribute::New(mesh_);
    if (sweptStencilDone_.size() != mesh_.nCells())
    {
        sweptStencilDone_.reset();
        sweptStencilDone_.resize(mesh_.nCells());
    }
    bool anyNew = false;
    forAll(donors, celli)
    {
        if (donors[celli] && !sweptStencilDone_.test(celli))
        {
            anyNew = true;
            break;
        }
    }
    reduce(anyNew, orOp<bool>());
    if (anyNew)
    {
        const labelListList& cc = mesh_.cellCells();
        boolList zone(donors);
        for (label layer = 0; layer < 2; ++layer)
        {
            const boolList prev(zone);
            forAll(prev, celli)
            {
                if (!prev[celli]) continue;
                for (const label nb : cc[celli]) zone[nb] = true;
            }
        }
        zd.updateStencil(zone);
        forAll(zone, celli)
        {
            if (zone[celli]) sweptStencilDone_.set(celli);
        }
    }
    const labelListList& stencil = zd.getStencil();
    const globalIndex& gNum = zd.globalNumbering();
    const auto t2 = now();

    // the remote entries (cells, or boundary faces of another rank, which
    // the owner answers with an empty entry)
    labelHashSet need;
    forAll(donors, celli)
    {
        if (!donors[celli]) continue;
        for (const label g : stencil[celli])
        {
            if (!gNum.isLocal(g)) need.insert(g);
        }
    }
    const labelList needSorted(need.sortedToc());

    // one transfer: the map from the global indices (the request round of
    // the mapDistributeBase constructor; kept while no rank's set of remote
    // cells changes), the packed cells sent per rank
    bool changed = (!sweptMap_ || needSorted != sweptPrevNeed_);
    reduce(changed, orOp<bool>());
    if (changed)
    {
        labelList compact(needSorted);
        List<Map<label>> compactMap;
        sweptMap_.reset(new mapDistributeBase(gNum, compact, compactMap));
        sweptPrevNeed_ = needSorted;
        sweptPrevCompact_ = std::move(compact);
    }
    const mapDistributeBase& map = sweptMap_();
    const labelList& compact = sweptPrevCompact_;
    const auto t3 = now();
    List<scalarList> data(gNum.localSize());
    DynamicList<scalar> buf(512);
    for (const labelList& sub : map.subMap())
    {
        for (const label idx : sub)
        {
            if (!data[idx].empty()) continue;
            if (idx < mesh_.nCells())
            {
                packSweptCell(idx, Up, buf);
                data[idx] = buf;
            }
            else
            {
                data[idx] = scalarList(1, Zero);
            }
        }
    }
    const auto t4 = now();
    map.distribute(data);
    const auto t5 = now();

    // unpack: the remote cells by ascending global index, and per donor the
    // indices of its remote ring cells
    Map<label> slotOf(2*needSorted.size());
    sweptRemote_.setSize(needSorted.size());
    label n = 0;
    forAll(needSorted, i)
    {
        davof::remoteCell rc;
        if (davof::unpackRemoteCell(data[compact[i]], rc))
        {
            rc.gid = needSorted[i];
            slotOf.insert(needSorted[i], n);
            sweptRemote_[n++] = std::move(rc);
        }
    }
    sweptRemote_.setSize(n);
    nSweptRemoteCells_ = n;
    forAll(donors, celli)
    {
        if (!donors[celli]) continue;
        DynamicList<label> ks;
        for (const label g : stencil[celli])
        {
            const auto iter = slotOf.cfind(g);
            if (iter.good()) ks.append(iter.val());
        }
        if (ks.size())
        {
            Foam::sort(ks);
            sweptRemoteRing_.insert(celli, labelList(std::move(ks)));
        }
    }
    const auto t6 = now();
    sweptExchangeTime_ = secs(t0, t6);
    sweptExchangeParts_[0] = secs(t0, t1);   // wait at the barrier
    sweptExchangeParts_[1] = secs(t1, t2);   // stencil
    sweptExchangeParts_[2] = secs(t2, t3);   // map
    sweptExchangeParts_[3] = secs(t3, t4);   // pack
    sweptExchangeParts_[4] = secs(t4, t5);   // transfer
    sweptExchangeParts_[5] = secs(t5, t6);   // unpack
}


Foam::label Foam::advection::isoAdvection::classifyPreImage
(
    const label donor,
    const pointField& q,
    const point& xq
) const
{
    const labelListList& cellPoints = mesh_.cellPoints();
    const labelListList& pointCells = mesh_.pointCells();
    const scalar h = Foam::cbrt(mesh_.V()[donor]);
    // the start-of-step volume fractions
    const scalarField& a0 = alpha1_.oldTime().primitiveField();

    // the class of one point: 0 donor, 1 interface, 2 full, 3 empty, 4 outside
    auto classOf = [&](const point& x) -> label
    {
        if (inCellTol(mesh_, x, donor, 1e-9, h)) return 0;
        for (const label pi : cellPoints[donor])
        {
            for (const label c : pointCells[pi])
            {
                if (c == donor || !inCellTol(mesh_, x, c, 1e-9, h)) continue;
                if (hasSurface_.test(c)) return 1;
                if (a0[c] >= 1 - 1e-8) return 2;
                if (a0[c] <= 1e-8) return 3;
                return 1;   // a partial cell without a surface counts as interface
            }
        }
        return 4;
    };
    // the priority 1 > 2 > 3 > 4 > 0 over the points
    static const label rank[5] = {0, 4, 3, 2, 1};
    label best = classOf(xq);
    forAll(q, i)
    {
        const label c = classOf(q[i]);
        if (rank[c] > rank[best]) best = c;
    }
    return best;
}


// ************************************************************************* //
