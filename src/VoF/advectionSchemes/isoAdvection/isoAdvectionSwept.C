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
    to the one-donor clip, counted; the cap set is read from the final
    clipped polygons and a boundary triangle lying in a cell plane with the
    swept volume outside that cell is a contact, so the fallback does not
    depend on the order of a cell's faces, davofSweptPieces.H, 2026-11-02),
    the largest closure residual of the
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

    - CONSISTENCY UNDER OBLIQUE TRAVEL (2026-11-02, both keys default false,
      the default path unchanged). A region of a cell lies in the swept
      volumes of two faces when the travel is oblique to the mesh (a corner
      of C in those of C|N and N|NE under north-east travel); the receiving
      cell balances only if both faces give that region the same liquid.
      Two rules can differ: the faces outside the candidate set keep the
      algebraic upwind flux alpha_donor phi_f dt, and a wisp's pieces are
      wet or dry by alpha >= 1/2. sweptConsistencyLog (serial) probes every
      algebraic face whose vertices touch a cell with a surface or another
      alpha, cuts its swept volume into the ring pieces and compares the
      consistent rule (a surface clips, any other cell gives alpha_c
      V_piece) with the algebraic flux; it sums the wisp pieces' mismatch of
      the computed faces and writes the per-cell residue of the step
      (sweptConsistencyProbe). sweptConsistent computes every face whose
      donor shares a vertex with an active cell (a surface, an interface
      cell, a wisp; the vertex marks synchronised over the processor
      patches), with or without a surface, and gives a piece in a cell
      without a surface alpha_c V_piece and alpha_c times its pre-image
      area: the piecewise algebraic value, the same in every face that
      sweeps the piece; a face without a surface whose pieces need the rank
      fallback takes the algebraic value (counted).

    - THE VOLUME CORRECTION (key sweptVolumeCorrection, 2026-11-02/03;
      without the key: scale where sweptScaleToPhi is true, none otherwise).
      The traced swept volume V_g = (1/3) sum_T x_T.a_T of the closed list
      differs from the flux volume V_f = |phi_f| dt (phi of the half step,
      setFlow at t + dt/2, the temporal pair of the midpoint trace) under a
      non-uniform velocity. scale multiplies the liquid volume by V_f / V_g
      (it divides by V_g, which goes to zero for a face twisted near a
      rotation axis); pyramid traces only the face points, starts the apex
      of the pre-image fan at the mean of the traced vertices and moves it
      along the pre-image's outward area vector A_pre by
      delta = 3 (V_f - V_g)/|A_pre|, the volume correction of the unsplit
      face-vertex flux calculation (UFVFC; see sweptTriangles); the corrected
      list is clipped and capped as before, and the face update uses the
      corrected fan (its total area vector is unchanged). Logged per step:
      max |V_f - V_g|/V_f before the correction, max |delta|/h, the twisted
      faces (V_g below 0.1 of the sum of the unsigned tetrahedra about x_f),
      the faces without a pre-image area, the largest per-axis displacement
      of a traced vertex over h.

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
#include "volFields.H"
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
    // The same trace through the cache of a mesh point (sweptTraceCache,
    // 2026-11-04): the midpoint x_m and the cell c that held it with the
    // velocity there are kept per point; another face holding the point
    // repeats only the cell search from its own donor (the ring order of
    // that donor, as the per-face trace) and takes the kept velocity when
    // the same cell results, so that the traced point is bit for bit the
    // per-face one; a different cell (a midpoint within the tolerance of a
    // cell face, claimed by two cells) interpolates anew.
    auto traceBackCached = [&](const label pi) -> point
    {
        const point& x = P[pi];
        const vector& u1 = Up[pi];
        const point xm = x - 0.5*dt*u1;
        const label c = cellOf(xm);
        if (traceDone_.test(pi) && traceCell_[pi] == c)
        {
            ++nSweptTraceHits_;
            if (c <= -2) sweptFaceRemote_ = true;
            if (c >= 0 || c <= -2) return x - dt*traceU_[pi];
            ++nSweptEulerTrace_;
            return x - dt*u1;
        }
        ++nSweptTraceMisses_;
        bool found = (c >= 0 || c <= -2);
        vector um(Zero);
        if (c >= 0) um = UInterp.interpolate(xm, c);
        else if (c <= -2) { sweptFaceRemote_ = true; um = davof::interpolateRemote(sweptRemote_[-2 - c], xm); }
        traceCell_[pi] = c; traceU_[pi] = um; traceDone_.set(pi);
        if (!found) ++nSweptEulerTrace_;
        return found ? (x - dt*um) : (x - dt*u1);
    };

    // The pre-image: the vertices by the point velocities, the centre by
    // the cell-point interpolate in the donor (the prism's displacement)
    const point xf = Cf[facei];
    q.setSize(F);
    forAll(f, i)
    {
        // sweptTraceCache (2026-11-04): a mesh point is traced once per
        // step (its trace depends on the point and the velocity only) and
        // reused by every other face that holds it; the cell search of the
        // midpoint starts from the first face's donor, so a midpoint within
        // the tolerance of a cell face can be interpolated in the other
        // cell: the same velocity in exact arithmetic, round-off apart
        const label pi = f[i];
        q[i] = sweptTraceCache_ ? traceBackCached(pi) : traceBack(P[pi], Up[pi]);
    }
    if (sweptVolumeCorrection_ == "pyramid")
    {
        // sweptVolumeCorrection pyramid (2026-11-03): only the face points
        // are traced; the apex of the pre-image fan starts at the mean of
        // the traced vertices (the translated face centre of a parallelogram
        // face under a uniform translation) and is placed by the pyramid
        // below
        xq = Zero;
        forAll(q, i) xq += q[i];
        xq /= scalar(F);
    }
    else
    {
        xq = traceBack(xf, UInterp.interpolate(xf, donor));
    }

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

    // the largest per-axis displacement of a traced vertex over h (the one
    // layer the pieces cover needs it below 1)
    if (sweeps)
    {
        forAll(f, i)
        {
            const vector dsp = P[f[i]] - q[i];
            sweptDispAxisMax_ = max(sweptDispAxisMax_, cmptMax(cmptMag(dsp))/h);
        }
    }

    vTotal = 0;
    if (sweeps)
    {
        // the signed volume about x_f (translation invariant; the local
        // origin avoids the cancellation of the global coordinates), and the
        // unsigned one (the sum of the tetrahedra's magnitudes)
        scalar sixV = 0, sixU = 0;
        forAll(T, k)
        {
            const scalar t6 = (T[k][0] - xf) & ((T[k][1] - xf) ^ (T[k][2] - xf));
            sixV += t6;
            sixU += mag(t6);
        }
        // the flux volume of the step: V_f = |phi_f| dt with phi of the half
        // step (advectorVoF evaluates setFlow at t + dt/2), the temporal
        // pair of the explicit midpoint rule of the trace above; no other
        // temporal scheme is mixed in (2026-11-02)
        const scalar Vf = mag(phiF)*dt;
        if (mag(Vf) > VSMALL)
        {
            sweptVolDefectMax_ = max(sweptVolDefectMax_, mag(Vf - sixV/6)/Vf);
        }
        if (sixV < 0.1*sixU) ++nSweptTwisted_;

        if (sweptVolumeCorrection_ == "pyramid")
        {
            // THE PYRAMID (sweptVolumeCorrection pyramid, 2026-11-02): the
            // apex x_f' of the pre-image fan (the triangles k % 4 == 1) is
            // moved along the pre-image polygon's area vector A_pre,
            // oriented outward of the swept volume, so that the closed list
            // holds V_f. The volume of the closed list is
            // V_g = (1/3) sum_T x_T . a_T (the divergence theorem,
            // int_V (1/3) div x dV = int_S (1/3) x.n dS, the sum used in
            // T. Tolle, D. Gruending, D. Bothe, T. Maric, triSurfaceImmersion:
            // Computing volume fractions and signed distances from
            // triangulated surfaces immersed in unstructured meshes, Comput.
            // Phys. Commun. 273 (2022) 108249, doi:10.1016/j.cpc.2021.108249,
            // arXiv:2101.08511), exact for a
            // thin or twisted swept volume as well; it is linear in the apex,
            // V(x_f' + delta m) = V(x_f') + (1/3) delta m . A_pre, because
            // the fan's area vector sum_i (1/2) (p_i' - c) x (p_{i+1}' - c)
            // does not depend on the apex c. With m = A_pre/|A_pre|,
            // delta = 3 (V_f - V_g)/|A_pre| closes the difference in closed
            // form: the correction of the fluxed volume of the unsplit
            // face-vertex flux calculation (UFVFC; T. Maric, H. Marschall,
            // D. Bothe, An enhanced un-split face-vertex flux-based VoF
            // method, J. Comput. Phys. 371 (2018) 967-993,
            // doi:10.1016/j.jcp.2018.03.048; no arXiv version found),
            // without a volumetric triangulation. Moving the apex along the swept face's normal,
            // and not in any other direction, keeps the pre-image fan from
            // folding over itself (a further self-intersection). The list is
            // kept in the orientation of phi_f (no flip: a twisted list near
            // a rotation axis has V_g < 0, and the pyramid takes it to V_f).
            vector Apre(Zero);
            forAll(f, i)
            {
                const FixedList<point, 3>& t = T[4*i + 1];
                Apre += 0.5*((t[1] - t[0]) ^ (t[2] - t[0]));
            }
            const scalar Amag = mag(Apre);
            if (Amag > SMALL*h*h)
            {
                const scalar delta = 3*(Vf - sixV/6)/Amag;
                xq += (delta/Amag)*Apre;
                forAll(f, i) T[4*i + 1][0] = xq;
                sweptPyrDeltaMax_ = max(sweptPyrDeltaMax_, mag(delta)/h);
                sixV = 0;
                forAll(T, k)
                {
                    sixV += (T[k][0] - xf) & ((T[k][1] - xf) ^ (T[k][2] - xf));
                }
            }
            else
            {
                ++nSweptPyrSkipped_;
            }
        }
        else if (sixV < 0)
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
    sweptFaceConsDelta_ = 0;

    // The consistency of the swept flux (2026-11-02): sweptConsistent gives
    // a piece in a cell without a surface alpha_c V_piece (consRule);
    // sweptConsistencyLog sums the consistent rule's liquid next to the used
    // one (threeVCons) and the wisp pieces' mismatch of the 0/1 rule against
    // alpha_c V_piece (the face's sums, scaled with the flux below)
    const bool consRule = sweptConsistent_;
    const bool consLog = sweptConsistencyLog_;
    scalar threeVCons = 0;
    label nWispF = 0;
    scalar wispSumF = 0, wispAbsF = 0, wispMaxF = 0;

    // the per-class statistics (sweptClassStats, 2026-11-05): the wall
    // clock around the calls, read only when the key is on
    const bool stats = sweptClassStats_;
    auto tick = [&]() -> double
    {
        return stats ? std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count() : 0.0;
    };
    double tStatVol = 0, tStatArea = 0;
    label nStatPieces = 0, nStatVolClips = 0, nStatAreaClips = 0, nStatTris = 0;

    List<FixedList<point, 3>> T;
    pointField q;
    point xq;
    label preSign;
    scalar vTotal;
    const double tTrace0 = tick();
    sweptTriangles(facei, donor, dt, phiF, sweeps, UInterp, Up, T, q, xq, preSign, vTotal);
    const double tStatTrace = tick() - tTrace0;
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
    // the longest triangle edge of the list: the size bound of every polygon
    // a clip can produce (sweptPreTest strict)
    scalar LmaxT = 0;
    forAll(T, k)
    {
        sumMagT += mag(0.5*((T[k][1] - T[k][0]) ^ (T[k][2] - T[k][0])));
        LmaxT = max(LmaxT, max(mag(T[k][1] - T[k][0]), max(mag(T[k][2] - T[k][1]), mag(T[k][0] - T[k][2]))));
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
    const pointField& P = mesh_.points();
    const point xf = Cf[facei];
    // The donor fast path (sweptDonorFastPath, 2026-11-05): the vertex set of
    // the swept volume (as the pre-test's) tested against the donor's planes
    // first; strictly inside every one by more than 1e-10 L_max, the face
    // needs the donor alone (isoAdvection.H)
    bool donorOnly = false;
    DynamicList<point> Vsw(2*F + 2);
    if (sweptDonorFastPath_ || sweptPreTest_ != "none" || sweptStateShortcut_ || sweptPreTestCheck_)
    {
        Vsw.append(xf);
        for (const label pi : mesh_.faces()[facei]) Vsw.append(P[pi]);
        Vsw.append(xq);
        forAll(q, i) Vsw.append(q[i]);
    }
    if (sweptDonorFastPath_ && !sweptPreTestCheck_ && !sweptStateShortcut_ && !sweptFaceDump_)
    {
        donorOnly = true;
        const scalar margin = 1e-10*LmaxT;
        const labelList& cf = cells[donor];
        forAll(cf, j)
        {
            const vector nOut = (own[cf[j]] == donor) ? Sf[cf[j]] : -Sf[cf[j]];
            const vector n = nOut/max(mag(nOut), VSMALL);
            const point& x0 = Cf[cf[j]];
            for (const point& v : Vsw)
            {
                if (((v - x0) & n) >= -margin) { donorOnly = false; break; }
            }
            if (!donorOnly) break;
        }
        if (donorOnly) ++nSweptDonorFast_;
    }
    boundBox sweptBox(T[0][0], T[0][0]);
    forAll(T, k)
    {
        sweptBox.add(T[k][0]); sweptBox.add(T[k][1]); sweptBox.add(T[k][2]);
    }
    sweptBox.inflate(1e-9);
    // the cell's box from the per-step cache (curvedFlux, 2026-11-04; the
    // same boundBox(points, cellPoints) as built here per face before)
    auto touches = [&](const label c) -> bool
    {
        return cellBoxes_[c].overlaps(sweptBox);
    };
    DynamicList<label> ring(32);
    ring.append(donor);
    if (!donorOnly)
    {
        for (const label pi : cellPoints[donor])
        {
            for (const label c : pointCells[pi])
            {
                if (c != donor && !ring.found(c) && touches(c)) ring.append(c);
            }
        }
        // in parallel the ring cells of other ranks (sweptExchange), after the
        // local ones, as -2 - k (k their index in sweptRemote_)
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
    // the donor's centre: the side of a cell plane on which a boundary
    // triangle lying in that plane bounds the swept volume (the contact
    // rule of sweptPieceVolume, davofSweptPieces.H, 2026-11-02)
    const point xDonor = mesh_.cellCentres()[donor];

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

    // the debug face (key sweptDebugPoint, 2026-10-31): every step, on the
    // rank that computes it, the face whose centre lies within 1e-9 of the
    // point writes its pre-image and every ring cell's state, piece and
    // pre-image part to postProcessing/sweptDebugFace/<step>_<rank>.txt
    const bool debugFace = sweptDebug_ && mag(Cf[facei] - sweptDebugPoint_) < 1e-9;
    autoPtr<OFstream> dbg;
    if (debugFace)
    {
        const fileName dir(mesh_.time().globalPath()/"postProcessing"/"sweptDebugFace");
        mkDir(dir);
        dbg.reset(new OFstream(dir/(Foam::name(mesh_.time().timeIndex()) + "_" + Foam::name(Pstream::myProcNo()) + ".txt")));
        dbg().precision(17);
        dbg() << "face " << facei << " Cf " << Cf[facei] << " Sf " << Sf[facei] << " donor " << donor
            << " phi " << phiF << " sweeps " << sweeps << " vTotal " << vTotal << nl
            << "xq " << xq << nl << "q " << q << nl;
    }

    // The pieces. Two passes over the ring (2026-11-04): the first takes
    // every cell's half-spaces and state and runs the pre-test
    // (sweptPreTest) that decides whether the swept volume can hold a piece
    // of the cell at all; the state shortcut (sweptStateShortcut) then
    // settles a face whose touched cells are all full or all empty without
    // a clip; the second pass clips the touched cells of the other faces
    // exactly as before (the same ring order, the same sums).
    scalar threeV = 0;
    vector aPreWet(Zero);
    scalar coverage = 0, coveragePre = 0;
    bool rankFallback = false;
    List<davof::halfSpace> planes;

    // The pre-test (sweptPreTest, 2026-11-04). The vertex set of the swept
    // volume: x_f, the face points, the apex x_f' (where the pyramid moved
    // it) and the traced points; a cell whose bounding box overlaps the
    // swept volume's (the ring) still misses it when one of the cell's
    // planes has the whole vertex set on its outer side. strict: skipped
    // only when every vertex lies beyond 1e-10 L (L the longest triangle
    // edge, the bound of every polygon size the clips meet, the tolerance
    // of the contact and coplanar tests): such a cell's clips remove every
    // triangle entirely (all signed distances positive, no triangle within
    // the coplanar tolerance), so the result is the same bit for bit.
    // touching (the default): skipped also when the vertices lie on the
    // plane (signed distance >= 0, the cells sharing only a face, an edge
    // or a point with the swept volume: 18 to 27 ring cells against 1.4
    // pieces per face), unless a pre-image triangle lies in that plane (the
    // coplanar rule of the face update keeps it). In exact arithmetic such
    // a cell holds no volume and no area; the clips of the unskipped code
    // leave zero-area polygons whose round-off area (1e-20 of the face)
    // entered the wetted area vector, so the face fractions move by
    // round-off (verified 2026-11-04 on the 20-step copies). none: no test.
    const bool preTest = (sweptPreTest_ != "none") && !donorOnly;
    const bool preStrict = (sweptPreTest_ == "strict");
    // a ring cell's planes (flat) and state
    struct ringState
    {
        label c; const davof::remoteCell* rc; scalar a0c; bool hasS, wetPlane, wet, wisp, touched;
        davof::halfSpace plic; label p0, np;
        label sepPlane = -1; scalar sepMin = 0;   // the separating plane of a skipped cell (the check)
    };
    DynamicList<ringState> states(ring.size());
    DynamicList<davof::halfSpace> ringPlanes(6*ring.size());

    // the vertices of the swept volume held by a touched ring cell (the
    // coverage condition of the state shortcut, 2026-11-04): with every
    // vertex inside a touched cell the swept volume, within the convex hull
    // of its vertices, lies in the union of the touched cells when that
    // union is convex (the vertex ring of a hex cell, a 3 x 3 x 3 block;
    // on a mesh whose ring is not convex the condition is not sufficient
    // and the shortcut should stay off), so the pieces cover it entirely
    // and the shortcut's V_f and A_pre are the pieces' sums up to round-off.
    // A face whose pre-image leaves the mesh (the inflow boundary) or the
    // ring fails it and keeps the pieces with their coverage accounting.
    List<bool> coveredV(Vsw.size(), false);
    const double tPre0 = tick();
    for (const label e : ring)
    {
        // a ring entry: a local cell c = e, or a cell of another rank
        // (e <= -2) whose packed planes and state are read instead
        ringState st;
        st.rc = (e <= -2) ? &sweptRemote_[-2 - e] : nullptr;
        st.c = st.rc ? st.rc->gid : e;
        const davof::remoteCell* rc = st.rc;
        const label c = st.c;
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
        // the donor fast path: no cell planes (their clips would be no-ops)
        st.p0 = ringPlanes.size(); st.np = donorOnly ? 0 : planes.size();
        if (!donorOnly) for (const davof::halfSpace& hs : planes) ringPlanes.append(hs);
        // the cell's state and plane
        const scalar a0c = rc ? rc->alpha0 : a0[c];
        const bool hasS = rc ? rc->hasSurface : hasSurface_.test(c);
        bool wetPlane = false, wet = true, wisp = false;
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
        else if (a0c >= 1 - sweptWispTol_)
        {
            wet = true;
        }
        else if (a0c <= sweptWispTol_)
        {
            wet = false;
        }
        else
        {
            // a partial cell without a surface (a wisp): wet by its alpha
            // (sweptConsistent: alpha_c V_piece, below)
            ++nSweptPartialCells_;
            wet = (a0c >= 0.5);
            wisp = true;
            // sweptWispSurface neighbour (2026-11-02): the plane of the
            // wisp's closest face neighbour with a surface clips its piece
            // and its part of the face update, as the second ring gives a
            // face without a surface its plane (a plane cutting a small
            // corner of the cell, |m_c| below the interface threshold);
            // anchored (2026-11-03): the plane with the neighbour's normal
            // positioned by the wisp's alpha (curvedFlux)
            const bool wispHas =
                rc
              ? rc->wispSurface
              : (hasWispSurface_.size() == mesh_.nCells() && hasWispSurface_.test(c));
            if (sweptWispSurface_ != "none" && wispHas)
            {
                vector nRaw; point x0;
                planeOf(rc ? rc->S : wispSurfaces_[c], c, h, sweptPlaneTol_, "sweptFaceCells", nRaw, x0);
                plic.n = nRaw/max(mag(nRaw), VSMALL);
                plic.x0 = x0;
                wetPlane = true;
                wisp = false;
                ++nSweptWispClipped_;
            }
        }
        st.a0c = a0c; st.hasS = hasS; st.wetPlane = wetPlane; st.wet = wet; st.wisp = wisp; st.plic = plic;
        // the pre-test
        st.touched = true;
        if (preTest)
        {
            forAll(planes, i)
            {
                const davof::halfSpace& hs = planes[i];
                scalar smin = GREAT;
                for (const point& v : Vsw) smin = min(smin, (v - hs.x0) & hs.n);
                bool separated;
                if (preStrict)
                {
                    separated = (smin > 1e-10*LmaxT);
                }
                else
                {
                    separated = (smin >= 0);
                    if (separated)
                    {
                        // a pre-image triangle lying in the plane keeps the
                        // face update's coplanar rule: no skip
                        forAll(Tpre, j)
                        {
                            const FixedList<point, 3>& t = Tpre[j];
                            const scalar L = max(mag(t[1] - t[0]), mag(t[2] - t[0]));
                            scalar sMax = 0;
                            for (label k = 0; k < 3; ++k) sMax = max(sMax, mag((t[k] - hs.x0) & hs.n));
                            if (L > 0 && sMax <= 1e-10*L) { separated = false; break; }
                        }
                    }
                }
                if (separated) { st.touched = false; st.sepPlane = i; st.sepMin = smin; break; }
            }
        }
        if (st.touched) ++nSweptPreClipped_; else ++nSweptPreSkipped_;
        if (st.touched && (sweptStateShortcut_ || sweptPreTestCheck_))
        {
            forAll(Vsw, vi)
            {
                if (coveredV[vi]) continue;
                bool in = true;
                forAll(planes, i)
                {
                    if (((Vsw[vi] - planes[i].x0) & planes[i].n) > 0) { in = false; break; }
                }
                if (in) coveredV[vi] = true;
            }
        }
        states.append(st);
    }
    const double tStatPre = tick() - tPre0;
    bool covered = (sweptStateShortcut_ || sweptPreTestCheck_);
    for (const bool b : coveredV) covered = covered && b;

    // The state shortcut (sweptStateShortcut, 2026-11-04): every touched
    // cell full (alpha^n exactly 1, no surface) or every one empty
    // (exactly 0): the pieces would all be wet (or all dry), summing to the
    // swept volume (or to nothing), and the pre-image fan wholly wet (or
    // dry): V_liq = V_f (the scaling, or V_g) and the wetted area A_pre, or
    // zero. Exact 0 and 1 rather than the wisp tolerance: with
    // sweptConsistent a piece of a cell at alpha = 1 - 1e-10 carries
    // alpha V_piece, which the shortcut would round to V_piece; the faces
    // the tolerance would take in addition are counted apart. The empty
    // shortcut gives the unclipped sums bit for bit (nothing is added);
    // the full one differs by round-off (the sum of the pieces against
    // V_g, the scaled V_liq against |phi_f| dt, the sum of the wetted parts
    // against A_pre).
    bool shortFull = false, shortEmpty = false;
    if (sweptStateShortcut_ && !covered) ++nSweptShortNotInside_;
    if (sweptStateShortcut_ && covered)
    {
        bool allFull = true, allEmpty = true, allFullTol = true, allEmptyTol = true, anyRemote = false;
        for (const ringState& st : states)
        {
            if (!st.touched) continue;
            if (st.rc) anyRemote = true;
            if (st.wetPlane) { allFull = allEmpty = allFullTol = allEmptyTol = false; break; }
            if (st.a0c != 1) allFull = false;
            if (st.a0c != 0) allEmpty = false;
            if (st.a0c < 1 - sweptWispTol_) allFullTol = false;
            if (st.a0c > sweptWispTol_) allEmptyTol = false;
        }
        shortFull = allFull; shortEmpty = allEmpty;
        if (allFullTol && !allFull) ++nSweptShortFullTol_;
        if (allEmptyTol && !allEmpty) ++nSweptShortEmptyTol_;
        if (shortFull || shortEmpty)
        {
            if (shortFull) ++nSweptShortFull_; else ++nSweptShortEmpty_;
            if (anyRemote) sweptFaceRemote_ = true;
            threeV = shortFull ? 3*vTotal : 0;
            aPreWet = shortFull ? aPreFull : vector(Zero);
            coverage = sumMagT; coveragePre = sumMagPre;
        }
        else
        {
            ++nSweptShortPieces_;
        }
    }

    // sweptPreTestCheck (2026-11-04, a debug key): the skipped cells are
    // clipped as well and a cell with a piece (coverage above 1e-12 of the
    // boundary, as the pieces count) is a violation of the pre-test; a
    // shortcut face runs the pieces too and their sums are compared with
    // the shortcut's (the pieces' values are kept, the counters report)
    const bool check = sweptPreTestCheck_;
    const bool shortcut = shortFull || shortEmpty;
    const scalar threeVsc = threeV; const vector aPreWetSc = aPreWet;
    if (check && shortcut) { threeV = 0; aPreWet = Zero; coverage = 0; coveragePre = 0; }
    for (const ringState& st : states)
    {
        if (!st.touched)
        {
            if (check)
            {
                const SubList<davof::halfSpace> planes(ringPlanes, st.np, st.p0);
                // a violation: a piece of nonzero volume, a nonzero wetted
                // area, or a covered part that is not a degenerate slab (the
                // zero-thickness piece of triangles lying in a cell plane,
                // counted apart: volume 0, nothing taken from it)
                bool viol = false, slab = false; scalar covV = 0, covA = 0;
                davof::pieceResult r;
                if (sweeps)
                {
                    r = davof::sweptPieceVolume(T, planes, st.wetPlane, st.plic, xf, &xDonor, sweptCapSet_ == "crossing");
                    covV = r.coverage;
                    // the pieces' own gate: a piece counts above 1e-12 of
                    // the boundary; below it nothing is added
                    if (r.coverage > 1e-12*sumMagT) { if (r.degenerate) slab = true; else viol = true; }
                }
                scalar cov = 0;
                const vector aw = davof::sweptPieceArea(Tpre, planes, st.wetPlane, st.plic, cov, nRef);
                covA = cov;
                if (cov > 1e-12*sumMagPre || magSqr(aw) > 0) viol = true;
                if (slab && !viol) ++nSweptPreSlab_;
                if (viol)
                {
                    ++nSweptPreViol_;
                    if (nSweptPreViol_ <= 40)
                    {
                        const davof::halfSpace& sp = planes[max(st.sepPlane, label(0))];
                        Pout<< "sweptPreTestCheck: step " << mesh_.time().timeIndex()
                            << " face " << facei << (facei >= mesh_.nInternalFaces() ? " (boundary)" : "")
                            << " donor " << donor << " cell " << st.c << (st.rc ? " (remote)" : "")
                            << " volume coverage " << covV/max(sumMagT, VSMALL)
                            << " piece volume/V_g " << (vTotal > VSMALL ? r.volume/vTotal : 0)
                            << " degenerate " << r.degenerate << " caps " << r.nCaps << " contacts " << r.nContact
                            << " area coverage " << covA/max(sumMagPre, VSMALL) << " |aw|/|A_pre| " << mag(aw)/max(mag(aPreFull), VSMALL)
                            << " alpha " << st.a0c << " surface " << st.hasS
                            << " separating plane " << st.sepPlane << " n " << sp.n << " smin/h " << st.sepMin/h
                            << " donor side " << ((xDonor - sp.x0) & sp.n)/h << endl;
                    }
                }
            }
            continue;
        }
        if (shortcut && !check) continue;
        const davof::remoteCell* rc = st.rc;
        const label c = st.c;
        const scalar a0c = st.a0c;
        const bool hasS = st.hasS;
        const bool wetPlane = st.wetPlane, wet = st.wet, wisp = st.wisp;
        const davof::halfSpace& plic = st.plic;
        const SubList<davof::halfSpace> planes(ringPlanes, st.np, st.p0);

        scalar dbgVol = 0, dbgCovV = 0;
        label dbgCaps = -1, dbgContact = 0;
        bool dbgOk = true, dbgDeg = false;
        if (sweeps)
        {
            const double tv0 = tick();
            const davof::pieceResult r =
                davof::sweptPieceVolume
                (
                    T, planes, wetPlane, plic, xf, &xDonor,
                    sweptCapSet_ == "crossing"
                );
            tStatVol += tick() - tv0;
            nStatVolClips += r.nClip; nStatTris += r.nTri;
            if (r.coverage > 1e-12*sumMagT) ++nStatPieces;
            dbgVol = r.volume;
            dbgCovV = r.coverage;
            dbgCaps = r.nCaps; dbgContact = r.nContact; dbgOk = r.ok; dbgDeg = r.degenerate;
            coverage += r.coverage;
            if (r.coverage > 1e-12*sumMagT)
            {
                if (rc) sweptFaceRemote_ = true;
                ++nSweptPieces_;
                if (!r.ok) rankFallback = true;
                if (r.degenerate) ++nSweptDegenerate_;
                if (consRule && !wetPlane)
                {
                    // a cell without a surface: alpha_c times the piece (the
                    // whole piece: no PLIC clip without a plane)
                    if (!r.degenerate)
                    {
                        threeV += 3*a0c*r.volume;
                        sweptPieceClosureMax_ = max(sweptPieceClosureMax_, r.closureResidual);
                    }
                }
                else if ((wetPlane || wet) && !r.degenerate)
                {
                    threeV += 3*r.volume;
                    sweptPieceClosureMax_ = max(sweptPieceClosureMax_, r.closureResidual);
                }
                if (consLog && !r.degenerate)
                {
                    threeVCons += 3*(wetPlane ? scalar(1) : a0c)*r.volume;
                    if (wisp)
                    {
                        // the used weight (0/1, or alpha_c under the
                        // consistent rule) minus alpha_c, times the piece
                        const scalar wUsed = consRule ? a0c : (wet ? scalar(1) : scalar(0));
                        const scalar mis = (wUsed - a0c)*r.volume;
                        ++nWispF;
                        wispSumF += mis;
                        wispAbsF += mag(mis);
                        wispMaxF = max(wispMaxF, mag(mis));
                    }
                }
                nSweptClipTris_ += r.nCut;
            }
        }
        scalar cov = 0;
        vector ac(Zero);
        bool cop = false;
        const double ta0 = tick();
        const vector aw = davof::sweptPieceArea(Tpre, planes, wetPlane, plic, cov, nRef, &ac, &cop, stats ? &nStatAreaClips : nullptr);
        tStatArea += tick() - ta0;
        coveragePre += cov;
        if (rc && cov > 0) sweptFaceRemote_ = true;
        if (cop) ++nSweptCoplanar_;
        if (consRule && !wetPlane) aPreWet += a0c*aw;
        else if (wetPlane || wet) aPreWet += aw;
        if (debugFace && (cov > 0 || dbgCovV > 0 || dbgContact > 0))
        {
            dbg() << "cell " << c << " remote " << (rc ? 1 : 0) << " alpha0 " << a0c
                << " surface " << hasS << " wet " << wet << " plic.n " << plic.n
                << " plic.x0 " << plic.x0 << " pieceVolume " << dbgVol
                << " pieceCoverage " << dbgCovV << " nCaps " << dbgCaps
                << " ok " << dbgOk << " degenerate " << dbgDeg
                << " contacts " << dbgContact << " preCoverage " << cov
                << " aWet " << aw << " aCell " << ac << " coplanar " << cop
                << " planes";
            for (const davof::halfSpace& hs : planes) dbg() << ' ' << hs.n << ' ' << hs.x0;
            dbg() << nl;
        }
        if (dump && cov > 0)
        {
            word cls = (!rc && c == donor) ? "donor" : (wetPlane ? "interface" : (a0c >= 1 - 1e-8 ? "full" : (a0c <= 1e-8 ? "empty" : "wisp")));
            pieceRows.append({c, cls, (wetPlane || wet) ? aw : vector(Zero), ac, cop});
        }
    }

    // the per-class statistics (sweptClassStats): the class of the face from
    // the pre-test's touched cells, the sums of this face
    if (stats)
    {
        label nTouched = 0;
        for (const ringState& st : states) if (st.touched) ++nTouched;
        const label k = (facei >= mesh_.nInternalFaces()) ? 2 : (nTouched <= 1 ? 0 : 1);
        ++nClsFaces_[k]; nClsTouched_[k] += nTouched; nClsPieces_[k] += nStatPieces;
        nClsVolClips_[k] += nStatVolClips; nClsAreaClips_[k] += nStatAreaClips; nClsTris_[k] += nStatTris;
        tClsVol_[k] += tStatVol; tClsArea_[k] += tStatArea; tClsTrace_[k] += tStatTrace; tClsPre_[k] += tStatPre;
    }

    // the shortcut against the pieces (sweptPreTestCheck)
    if (check && shortcut)
    {
        bool viol = false;
        if (sweeps && vTotal > VSMALL)
        {
            if (shortFull && mag(threeV/3 - vTotal) > 1e-10*vTotal) viol = true;
            if (shortEmpty && mag(threeV) > 1e-12*vTotal) viol = true;
        }
        if (mag(aPreWet - aPreWetSc) > 1e-10*mag(aPreFull)) viol = true;
        if (viol)
        {
            ++nSweptShortViol_;
            if (nSweptShortViol_ <= 20)
            {
                Pout<< "sweptPreTestCheck: step " << mesh_.time().timeIndex() << " face " << facei
                    << (facei >= mesh_.nInternalFaces() ? " (boundary)" : "") << " donor " << donor << " sweeps " << sweeps
                    << (shortFull ? " full" : " empty") << " shortcut: pieces V_liq/V_g "
                    << (vTotal > VSMALL ? threeV/3/vTotal : 0) << " against " << (vTotal > VSMALL ? threeVsc/3/vTotal : 0)
                    << ", coverage " << coverage/max(sumMagT, VSMALL) << ", wetted area " << aPreWet << " against " << aPreWetSc
                    << ", ring " << ring.size() << " touched:";
                for (const ringState& st : states)
                {
                    if (!st.touched) continue;
                    const SubList<davof::halfSpace> planes(ringPlanes, st.np, st.p0);
                    davof::pieceResult r;
                    if (sweeps) r = davof::sweptPieceVolume(T, planes, st.wetPlane, st.plic, xf, &xDonor, sweptCapSet_ == "crossing");
                    Pout<< " [cell " << st.c << (st.rc ? "r" : "") << " alpha " << st.a0c
                        << " cov " << r.coverage/max(sumMagT, VSMALL) << " V/V_g " << (vTotal > VSMALL ? r.volume/vTotal : 0)
                        << " deg " << r.degenerate << " ok " << r.ok << " caps " << r.nCaps << " contacts " << r.nContact << "]";
                }
                Pout<< endl;
            }
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
        if (sweptNoDonorSurface_)
        {
            // no surface (sweptConsistent): the donor's alpha over the
            // pre-image
            return a0[donor]*aPreFull;
        }
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
        sweptFaceFallback_ = true;
        if (sweptNoDonorSurface_)
        {
            // a face without a surface (sweptConsistent): no plane to clip
            // by; the algebraic value, the donor's alpha (counted)
            ++nConsFallback_;
            af = a0[donor];
            aWet = a0[donor]*aPreFull;
            vLiq = sweeps ? a0[donor]*mag(phiF)*dt : 0;
            return;
        }
        sweptFaceFlux(facei, donor, S, dt, phiF, sweeps, UInterp, Up, af, aWet, vLiq);
        if (debugFace)
        {
            dbg() << "result (rank fallback, the one-donor clip) af " << af << " vLiq " << vLiq << nl;
        }
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
        // the full shortcut: the flux volume itself, not V_g (|phi_f| dt/V_g)
        if (shortFull) vLiq = mag(phiF)*dt;
    }
    if (consLog && sweeps)
    {
        // the consistent rule's liquid, scaled as the flux, and the face's
        // wisp pieces (sweptConsistencyLog)
        const scalar scale =
            (sweptScaleToPhi_ && vTotal > VSMALL) ? mag(phiF)*dt/vTotal : scalar(1);
        sweptFaceConsDelta_ = vLiq - scale*threeVCons/3;
        nConsWispPieces_ += nWispF;
        consWispMisSum_ += scale*wispSumF;
        consWispMisAbs_ += scale*wispAbsF;
        consWispMisMax_ = max(consWispMisMax_, scale*wispMaxF);
    }
    const scalar afPieces = fractionOf(aPreWet);
    const vector wDonor = (dump || sweptCellsFaceUpdate_ == "donor") ? donorWet() : vector(Zero);
    if (sweptCellsFaceUpdate_ == "donor")
    {
        aPreWet = wDonor;
    }
    aWet = aPreWet;
    af = fractionOf(aWet);
    if (debugFace)
    {
        dbg() << "result af " << af << " vLiq " << vLiq << " threeV " << threeV
            << " coverage " << coverage << " coveragePre " << coveragePre
            << " rankFallback " << rankFallback << nl;
    }

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


void Foam::advection::isoAdvection::sweptConsistencyProbe
(
    const scalar dt,
    const interpolationCellPoint<vector>& UInterp,
    const pointVectorField& Up
)
{
    // The diagnostic of sweptConsistencyLog (2026-11-02, serial). Every
    // internal face whose flux curvedFlux did not set keeps the algebraic
    // upwind flux alpha_donor phi_f dt. Its swept volume (the prism f - d
    // under a uniform velocity) is cut into the pieces of the donor's vertex
    // ring as in sweptFaceCellsFlux, and its liquid by the consistent rule
    // (a cell with a surface clipped by its plane, any other cell alpha_c
    // V_piece; scaled by |phi_f| dt / V_swept with sweptScaleToPhi) is
    // compared with the algebraic flux. Only faces with a cell of another
    // alpha or with a surface among the cells sharing a vertex of the face
    // are probed (the swept volume of a face reaches only those cells while
    // the displacement per step is below h along each axis; for the others
    // every piece carries alpha_donor and the two values agree). The
    // counters that sweptTriangles touches are restored.
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    const label nInt = mesh_.nInternalFaces();
    const faceList& faces = mesh_.faces();
    const cellList& cells = mesh_.cells();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();
    const vectorField& C = mesh_.cellCentres();
    const scalarField& V = mesh_.V();
    const pointField& P = mesh_.points();
    const labelListList& cellPoints = mesh_.cellPoints();
    const labelListList& pointCells = mesh_.pointCells();
    const scalarField& a0 = alpha1_.oldTime().primitiveField();
    const scalarField& phiIn = phi_.primitiveField();
    const scalarField& dVfIn = dVf_.primitiveField();
    const scalar tol = sweptConsistentTol_;

    const label s1 = nSweptFaces_, s2 = nSweptFlips_, s3 = nSweptOutOfDonor_;
    const label s4 = nSweptEulerTrace_;
    const scalar s5 = sweptClosureMax_;
    const bool s6 = sweptFaceRemote_;

    nConsProbed_ = nConsReach_ = nConsReachInterface_ = nConsReachWisp_ = 0;
    nConsReachOpposite_ = nConsProbeFallback_ = 0;
    consAlgMisSum_ = consAlgMisAbs_ = consAlgMisMax_ = 0;

    // the state of a cell: 0 a surface, 1 a wisp, 2 full, 3 empty
    auto stateOf = [&](const label c) -> label
    {
        if (hasSurface_.test(c)) return 0;
        if (a0[c] >= 1 - tol) return 2;
        if (a0[c] <= tol) return 3;
        return 1;
    };

    List<FixedList<point, 3>> T;
    pointField q;
    point xq;
    label preSign;
    scalar vTotal;
    List<davof::halfSpace> planes;
    DynamicList<label> ring(32);
    for (label facei = 0; facei < nInt; ++facei)
    {
        if (alphafSetFace_.test(facei)) continue;
        const scalar phiF = phiIn[facei];
        if (phiF == 0) continue;
        const label donor = (phiF > 0) ? own[facei] : nei[facei];
        const label other = (phiF > 0) ? nei[facei] : own[facei];

        // the prefilter: a cell sharing a vertex of the face with a surface
        // or another alpha
        bool differs = false;
        for (const label pi : faces[facei])
        {
            for (const label c : pointCells[pi])
            {
                if (c == donor) continue;
                if (hasSurface_.test(c) || a0[c] != a0[donor]) { differs = true; break; }
            }
            if (differs) break;
        }
        if (!differs) continue;

        const vector d = UInterp.interpolate(Cf[facei], donor)*dt;
        const bool sweeps = mag(d & Sf[facei]) > 1e-10*mag(d)*mag(Sf[facei]);
        if (!sweeps) continue;

        sweptTriangles(facei, donor, dt, phiF, sweeps, UInterp, Up, T, q, xq, preSign, vTotal);
        scalar sumMagT = 0;
        forAll(T, k)
        {
            sumMagT += mag(0.5*((T[k][1] - T[k][0]) ^ (T[k][2] - T[k][0])));
        }
        boundBox sweptBox(T[0][0], T[0][0]);
        forAll(T, k)
        {
            sweptBox.add(T[k][0]); sweptBox.add(T[k][1]); sweptBox.add(T[k][2]);
        }
        sweptBox.inflate(1e-9);
        ring.clear();
        ring.append(donor);
        for (const label pi : cellPoints[donor])
        {
            for (const label c : pointCells[pi])
            {
                if (c != donor && !ring.found(c) && boundBox(P, cellPoints[c], false).overlaps(sweptBox))
                {
                    ring.append(c);
                }
            }
        }

        const label sDonor = stateOf(donor);
        scalar threeV = 0;
        bool fallback = false;
        label reach = -1;   // 0 surface, 1 wisp, 2 a pure cell of another state
        for (const label c : ring)
        {
            const labelList& cf = cells[c];
            planes.setSize(cf.size());
            forAll(cf, j)
            {
                const vector nOut = (own[cf[j]] == c) ? Sf[cf[j]] : -Sf[cf[j]];
                planes[j].n = nOut/max(mag(nOut), VSMALL);
                planes[j].x0 = Cf[cf[j]];
            }
            const bool hasS = hasSurface_.test(c);
            davof::halfSpace plic;
            plic.n = Zero; plic.x0 = Zero;
            if (hasS)
            {
                vector nRaw; point x0;
                planeOf(cellSurfaces_[c], c, Foam::cbrt(V[donor]), sweptPlaneTol_, "sweptFaceCells", nRaw, x0);
                plic.n = nRaw/max(mag(nRaw), VSMALL);
                plic.x0 = x0;
            }
            const davof::pieceResult r =
                davof::sweptPieceVolume
                (
                    T, planes, hasS, plic, Cf[facei], &C[donor],
                    sweptCapSet_ == "crossing"
                );
            if (r.coverage <= 1e-12*sumMagT) continue;
            if (!r.ok) fallback = true;
            if (r.degenerate) continue;
            threeV += 3*(hasS ? scalar(1) : a0[c])*r.volume;
            if (c == donor) continue;
            // does the swept volume reach c (its whole piece, not its
            // liquid), and is c in another state than the donor?
            scalar vPiece = r.volume;
            if (hasS)
            {
                vPiece = davof::sweptPieceVolume
                (
                    T, planes, false, plic, Cf[facei], &C[donor], sweptCapSet_ == "crossing"
                ).volume;
            }
            if (vPiece <= 1e-9*max(vTotal, VSMALL)) continue;
            const label sc = stateOf(c);
            label k = -1;
            if (sc == 0) k = 0;
            else if (sc == 1) k = 1;
            else if (sc != sDonor) k = 2;
            if (k >= 0 && (reach < 0 || k < reach)) reach = k;
        }
        if (fallback)
        {
            ++nConsProbeFallback_;
            continue;
        }
        scalar vGeo = threeV/3;
        if (sweptScaleToPhi_ && vTotal > VSMALL)
        {
            vGeo *= mag(phiF)*dt/vTotal;
        }
        // the consistent rule's liquid minus the algebraic flux, out of the
        // donor; the residue (the used change of a cell's liquid minus the
        // consistent rule's) is +mis in the donor and -mis in the other cell
        const scalar mis = vGeo - mag(dVfIn[facei]);
        ++nConsProbed_;
        if (reach >= 0)
        {
            ++nConsReach_;
            if (reach == 0) ++nConsReachInterface_;
            else if (reach == 1) ++nConsReachWisp_;
            else ++nConsReachOpposite_;
        }
        consAlgMisSum_ += mis;
        consAlgMisAbs_ += mag(mis);
        consAlgMisMax_ = max(consAlgMisMax_, mag(mis));
        consResidueStep_[donor] += mis;
        consResidueStep_[other] -= mis;
    }

    nSweptFaces_ = s1; nSweptFlips_ = s2; nSweptOutOfDonor_ = s3;
    nSweptEulerTrace_ = s4; sweptClosureMax_ = s5; sweptFaceRemote_ = s6;

    // The residues of the step: the cells, the sum of |R|, the largest
    // |R|/V; the cumulative field sweptResidue (R/V summed over the steps)
    if (!consResiduePtr_)
    {
        consResiduePtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("sweptResidue", alpha1_.group()),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, Zero)
            )
        );
    }
    scalarField& Rcum = consResiduePtr_->primitiveFieldRef();
    const fileName dir(mesh_.time().globalPath()/"postProcessing"/"sweptConsistency");
    if (!consSummary_)
    {
        mkDir(dir);
        consSummary_.reset(new OFstream(dir/"summary.dat"));
        consSummary_().precision(10);
        consSummary_() << "# time step nProbed nReach nReachInterface nReachWisp"
            << " nReachPure nProbeFallback algMisSum algMisAbs algMisMax"
            << " nWispPieces wispMisSum wispMisAbs wispMisMax nResidueCells"
            << " residueAbs residueMaxRel nCandAdded nNoSurface nConsFallback"
            << " (volumes in m^3; residueMaxRel = max |R|/V)" << nl;
    }
    OFstream res(dir/("residue_" + Foam::name(mesh_.time().timeIndex()) + ".csv"));
    res.precision(10);
    res << "cell,x,y,z,alpha0,R" << nl;
    label nRes = 0;
    scalar resAbs = 0, resMax = 0;
    forAll(consResidueStep_, c)
    {
        const scalar R = consResidueStep_[c];
        if (mag(R) <= 1e-15*V[c]) continue;
        ++nRes;
        resAbs += mag(R);
        resMax = max(resMax, mag(R)/V[c]);
        Rcum[c] += R/V[c];
        res << c << ',' << C[c].x() << ',' << C[c].y() << ',' << C[c].z() << ','
            << a0[c] << ',' << R/V[c] << nl;
    }
    consSummary_() << mesh_.time().value() << ' ' << mesh_.time().timeIndex() << ' '
        << nConsProbed_ << ' ' << nConsReach_ << ' ' << nConsReachInterface_ << ' '
        << nConsReachWisp_ << ' ' << nConsReachOpposite_ << ' ' << nConsProbeFallback_ << ' '
        << consAlgMisSum_ << ' ' << consAlgMisAbs_ << ' ' << consAlgMisMax_ << ' '
        << nConsWispPieces_ << ' ' << consWispMisSum_ << ' ' << consWispMisAbs_ << ' '
        << consWispMisMax_ << ' ' << nRes << ' ' << resAbs << ' ' << resMax << ' '
        << nConsCandAdded_ << ' ' << nConsNoSurface_ << ' ' << nConsFallback_ << endl;
    Info<< "isoAdvection: swept consistency: algebraic faces probed " << nConsProbed_
        << " (reaching another state " << nConsReach_ << ": a surface "
        << nConsReachInterface_ << ", a wisp " << nConsReachWisp_ << ", a pure cell "
        << nConsReachOpposite_ << "; fallbacks " << nConsProbeFallback_
        << "), mismatch sum " << consAlgMisSum_ << " abs " << consAlgMisAbs_
        << " max " << consAlgMisMax_ << "; wisp pieces " << nConsWispPieces_
        << ", mismatch sum " << consWispMisSum_ << " abs " << consWispMisAbs_
        << " max " << consWispMisMax_ << "; residue cells " << nRes
        << ", abs " << resAbs << ", max |R|/V " << resMax << endl;
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
    // the surface flag: 1 its own surface, 2 a wisp's neighbour surface
    // (sweptWispSurface neighbour, 2026-11-02), 0 none
    const bool has = hasSurface_.test(celli);
    const bool wispHas =
        !has && hasWispSurface_.size() == mesh_.nCells() && hasWispSurface_.test(celli);
    v.append(has ? 1 : (wispHas ? 2 : 0));
    davof::packCellSurface
    (
        has ? cellSurfaces_[celli] : (wispHas ? wispSurfaces_[celli] : davof::cellSurface()),
        v
    );
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
    // the cells sent to the OTHER ranks only: the map's own-rank entry is
    // the identity over every local cell and boundary face (packing it
    // packed the whole mesh every step, 0.7 s of a 2.2 s advect call at
    // N = 128 on 8 ranks, fork 50b8bef)
    forAll(map.subMap(), proci)
    {
        if (proci == Pstream::myProcNo()) continue;
        for (const label idx : map.subMap()[proci])
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
