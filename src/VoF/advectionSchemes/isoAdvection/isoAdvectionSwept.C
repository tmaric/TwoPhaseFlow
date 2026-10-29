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
    (fluxIntegration sweptFace, 2026-10-28; the working note
    leia/docs/davof/davof-swept-volume, whose Python prototype swept_volume.py
    this routine follows).

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
    it is negative (counted). Every triangle is clipped against the donor
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
    the cap term vanishes identically (the hydrostatic form). dVf carries the
    sign of phi_f, as the prism flux does, and both cells of a face see the
    same volume (conservation by antisymmetry). The existing bounding and
    clipping of alpha are untouched.

    SCOPE AND LIMITS (stated here and in etc/README-davof.md):

    - PLANES ONLY. The formula needs x.n constant on the cap. The donor's
      surface must be a plane: a cellSurface without the paraboloid flag, or
      one whose curvature tensor is at round-off, max|b_ij| h <= sweptPlaneTol
      (default 1e-6; on the translated plane the ring tensor is 1e-10 /m at
      h = 2e-4 m, b h = 2e-14). A true paraboloid or a sphere surface is a
      FatalError, never a silent fallback: use clippingSurface plane, or
      fluxIntegration prism with the paraboloid.

    - ONE DONOR CELL. The whole swept volume is clipped by the donor's plane.
      The faces whose pre-image vertices (or centre) leave the donor cell,
      by the face-plane test with the tolerance 1e-9 h, are counted per step
      ("out of donor"); the multi-cell cut is not implemented.

    - TRAJECTORIES. Backward midpoint rule (explicit RK2, second order in
      dt): x' = x - dt U(x - dt/2 U(x)). The vertex velocity U(p_i) is the
      volPointInterpolation value at the vertex, the face centre's is the
      cell-point interpolate of U at x_f in the donor (the prism's
      displacement d = U(x_f) dt, so for a uniform U the pre-image is the
      prism's f - d and the swept volume IS the prism). The midpoint velocity
      is the cell-point interpolate in the cell that holds the midpoint (the
      donor, a cell sharing a vertex with it, else the mesh's octree search);
      a midpoint outside the mesh, or across a processor boundary, takes the
      Euler step x - dt U(x) (exact for a uniform U as well).

    - DIAGNOSTICS per step in the Info line: the faces computed, the sign
      flips, the faces out of the donor, the largest closure residual
      |sum_T a_T| / sum_T |a_T| (round-off for a closed list by construction)
      and the number of clipped triangles.

SourceFiles
    isoAdvectionSwept.C

\*---------------------------------------------------------------------------*/

#include "isoAdvection.H"
#include "davofSimplexGeometry.H"
#include "interpolationCellPoint.H"
#include "pointFields.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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
    const face& f = mesh_.faces()[facei];
    const pointField& P = mesh_.points();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();
    const labelList& own = mesh_.faceOwner();
    const cellList& cells = mesh_.cells();
    const label F = f.size();
    const scalar h = Foam::cbrt(mesh_.V()[donor]);

    af = 0;
    aWet = Zero;
    vLiq = 0;

    // The donor's plane: (x - x0).nRaw = 0, liquid where the value is
    // negative (the surface's own sign convention: e3 points out of the
    // liquid). A paraboloid with the slopes a1, a2 and a zero tensor is the
    // plane with the normal e3 - a1 e1 - a2 e2 through o + c e3.
    if (S.sphere)
    {
        FatalErrorInFunction
            << "fluxIntegration sweptFace: the surface of cell " << donor
            << " is a sphere (clippingSurface exactSphere). The swept-face"
            << " flux is implemented for a plane only (x.n is constant on the"
            << " cap); use clippingSurface plane or fluxIntegration prism."
            << exit(FatalError);
    }
    const scalar bMax = max(mag(S.b11), max(mag(S.b12), mag(S.b22)));
    if (S.paraboloid && bMax*h > sweptPlaneTol_)
    {
        FatalErrorInFunction
            << "fluxIntegration sweptFace: the surface of cell " << donor
            << " is a paraboloid (curvature tensor b11 " << S.b11 << " b12 "
            << S.b12 << " b22 " << S.b22 << " /m, max|b| h = " << bMax*h
            << " above sweptPlaneTol " << sweptPlaneTol_ << "). The"
            << " swept-face flux is implemented for a plane only (x.n is"
            << " constant on the cap); use clippingSurface plane, or"
            << " fluxIntegration prism with the paraboloid."
            << exit(FatalError);
    }
    const vector nRaw = S.paraboloid ? (S.e3 - S.a1*S.e1 - S.a2*S.e2) : S.e3;
    const point x0 = S.o + S.c*S.e3;
    auto sdist = [&](const point& x) -> scalar
    {
        return (x - x0) & nRaw;
    };

    // The point-in-cell test by the face planes (convex cells), with a
    // tolerance in units of the cell size
    auto inCell = [&](const point& x, const label c, const scalar tolRel) -> bool
    {
        for (const label cf : cells[c])
        {
            const vector nOut = (own[cf] == c) ? Sf[cf] : -Sf[cf];
            if (((x - Cf[cf]) & nOut) > tolRel*h*mag(nOut))
            {
                return false;
            }
        }
        return true;
    };

    // The cell that holds a point of the trajectory: the donor, else the
    // cells sharing a vertex with the donor (a midpoint lies within Co h/2
    // of the donor's boundary, inside that ring for Co < 1), else the mesh's
    // octree search (CELL_TETS; the tree and the tet base points are built
    // once per step in curvedFlux, on every rank; the face-plane search of
    // polyMesh::findCell scans every cell and is not used); -1 outside the
    // mesh or across a processor boundary
    const labelListList& cellPoints = mesh_.cellPoints();
    const labelListList& pointCells = mesh_.pointCells();
    auto cellOf = [&](const point& x) -> label
    {
        if (inCell(x, donor, 1e-12)) return donor;
        for (const label pi : cellPoints[donor])
        {
            for (const label c : pointCells[pi])
            {
                if (c != donor && inCell(x, c, 1e-12)) return c;
            }
        }
        return mesh_.findCell(x, polyMesh::CELL_TETS);
    };

    // The velocity at a point of the trajectory: the cell-point interpolate
    // in the cell that holds the point; nothing where no local cell holds it
    auto velocityNear = [&](const point& x, bool& found) -> vector
    {
        const label c = cellOf(x);
        found = (c >= 0);
        return found ? UInterp.interpolate(x, c) : vector(Zero);
    };

    // The backward midpoint trace from x with the velocity u1 at x
    auto traceBack = [&](const point& x, const vector& u1) -> point
    {
        const point xm = x - 0.5*dt*u1;
        bool found = false;
        const vector um = velocityNear(xm, found);
        return found ? (x - dt*um) : (x - dt*u1);
    };

    // The pre-image: the vertices by the point velocities, the centre by
    // the cell-point interpolate in the donor (the prism's displacement)
    const point xf = Cf[facei];
    pointField q(F);
    forAll(f, i)
    {
        q[i] = traceBack(P[f[i]], Up[f[i]]);
    }
    const point xq = traceBack(xf, UInterp.interpolate(xf, donor));
    if (qOut) *qOut = q;
    if (xqOut) *xqOut = xq;

    // The out-of-donor count: any traced point outside the donor, of the
    // faces that sweep a volume (a face parallel to the velocity carries no
    // flux, and its traced vertices leave the donor along its side faces by
    // construction; the studies of 2026-10-28 counted those faces too)
    if (sweeps)
    {
        bool outOfDonor = !inCell(xq, donor, 1e-9);
        forAll(q, i)
        {
            if (!inCell(q[i], donor, 1e-9)) outOfDonor = true;
        }
        if (outOfDonor) ++nSweptOutOfDonor_;
    }

    // The triangle list, outward for phi_f > 0 (the face fan along S_f, the
    // pre-image fan reversed, the lateral quads split along one diagonal)
    List<FixedList<point, 3>> T(4*F);
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
    label preSign = -1;
    auto flipAll = [&]()
    {
        forAll(T, k) Swap(T[k][1], T[k][2]);
        preSign = -preSign;
    };
    if (phiF < 0) flipAll();

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
        }
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

    // The face fraction of the pre-image in the face's orientation
    aWet = preSign*aPreWet;
    const vector Apre = preSign*aPreFull;
    const scalar A2 = magSqr(Apre);
    af = (A2 > VSMALL) ? (aWet & Apre)/A2 : 0;

    ++nSweptFaces_;
}


Foam::label Foam::advection::isoAdvection::classifyPreImage
(
    const label donor,
    const pointField& q,
    const point& xq
) const
{
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();
    const labelList& own = mesh_.faceOwner();
    const cellList& cells = mesh_.cells();
    const labelListList& cellPoints = mesh_.cellPoints();
    const labelListList& pointCells = mesh_.pointCells();
    const scalar h = Foam::cbrt(mesh_.V()[donor]);
    // the start-of-step volume fractions
    const scalarField& a0 = alpha1_.oldTime().primitiveField();

    auto inCell = [&](const point& x, const label c) -> bool
    {
        for (const label cf : cells[c])
        {
            const vector nOut = (own[cf] == c) ? Sf[cf] : -Sf[cf];
            if (((x - Cf[cf]) & nOut) > 1e-9*h*mag(nOut)) return false;
        }
        return true;
    };
    // the class of one point: 0 donor, 1 interface, 2 full, 3 empty, 4 outside
    auto classOf = [&](const point& x) -> label
    {
        if (inCell(x, donor)) return 0;
        for (const label pi : cellPoints[donor])
        {
            for (const label c : pointCells[pi])
            {
                if (c == donor || !inCell(x, c)) continue;
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
