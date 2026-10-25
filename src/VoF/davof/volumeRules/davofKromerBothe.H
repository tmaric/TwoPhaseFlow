/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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

Namespace
    Foam::davof

Description
    The liquid volume of a polyhedron truncated by a paraboloid, and the
    liquid area of a planar polygon cut by a quadric, by the face-based
    method of

        J. Kromer, D. Bothe, Third-order accurate initialization of volume
        fractions on unstructured meshes with arbitrary polyhedral cells,
        Journal of Computational Physics 475 (2023) 111840,
        doi:10.1016/j.jcp.2022.111840, arXiv:2111.01073,

    which builds on

        J. Kromer, D. Bothe, Highly accurate computation of volume fractions
        using differential geometry, Journal of Computational Physics 396
        (2019) 761-784, doi:10.1016/j.jcp.2019.07.005, arXiv:1805.03136.

    Their third order is the error of the osculating paraboloid that replaces
    the hypersurface; the DAVOF clipping surface (davofCurvedGeometry.H) IS a
    paraboloid in the frame (o; e1, e2, e3) of the cell's normal,

        phi(x) = w - c - a1 u - a2 v - 1/2 (b11 u^2 + 2 b12 u v + b22 v^2),
        liquid: phi < 0, i.e. w < h(t) = c + a.t + 1/2 t^T b t,  t = (u, v),

    so the only approximation left is the one-dimensional quadrature along
    the traces. The construction (their Section 3, in this notation):

    1. The divergence theorem with the field x - o (divergence 3) over the
       liquid part P^- of the polyhedron P,

           |P^-| = 1/3 [ sum_k ((x_k - o).n_k) A_k^wet
                       + int_{Gamma cap P} (x - o).n_Gamma dA ],

       the planar faces k (outward unit normal n_k, wetted area A_k^wet; the
       integrand is constant on a plane) and the paraboloid patch Gamma
       inside P, on which, with x - o = (t, h(t)) and n_Gamma dA = (-grad h,
       1) dt, the integrand is h - t.grad h = c - 1/2 t^T b t =: g(t), the
       slopes dropping out.

    2. Green's theorem in the (t1, t2) plane with H1 = int g dt1,

           H1(t) = c t1 - (b11 t1^3/6 + b12 t1^2 t2/2 + b22 t1 t2^2/2),

       turns the patch integral into the line integral  oint_{dS} H1 dt2
       around the projection S of the patch, counterclockwise about e3, whose
       boundary is the projection of the traces of Gamma on the faces of P
       inside the face polygons (their Eqs. 19-21 and 39).

    3. The wetted area of a face by Green's theorem as well, 1/2 oint
       (xi1 dxi2 - xi2 dxi1) around the liquid region in the face plane
       (coordinates xi): the liquid pieces of the polygon edges, exact, and
       the same trace arcs with the liquid on the left (their Eqs. 41-44,
       polygonal part and conic cap, here evaluated in one integral).

    The trace of a quadric on a plane is a conic. It is classified and
    parametrised as in their Eqs. 26-34 (rotation to the principal axes;
    ellipse, hyperbola, parabola, line) ONLY to order the crossings of the
    polygon boundary along the conic and to split long arcs; every arc is
    then integrated as the graph over its chord with the exact root of the
    quadratic along the chord's normal (the arc of a convex conic between two
    of its points is such a graph for a turning below 180 degrees, and arcs
    are split to at most 90 degrees of the ellipse parameter), so that the
    nearly degenerate conics of faces almost parallel to the paraboloid's
    axis, where the principal-form parametrisation loses digits, are
    integrated without loss. The quadrature along an arc is a 12-point
    Gauss-Legendre rule refined adaptively; nodes and weights by Newton on
    the Legendre polynomials (W. H. Press, S. A. Teukolsky, W. T. Vetterling,
    B. P. Flannery, Numerical Recipes, 3rd ed., Cambridge University Press,
    2007, Sect. 4.6, ISBN 978-0-521-88068-8). A fully enclosed trace (an
    ellipse inside a face, their Section 3.6) is integrated as four arcs.

    Non-planar mesh faces are integrated as the fan triangles about the face
    centre (planar); the interior fan edges cancel in the line integrals and
    carry the same face constant. The routines return false where the
    configuration is not resolved (a degenerate conic with crossings, a
    failed root along an arc, an area outside [0, A]); the caller then falls
    back to the split quadrature of davofCurvedGeometry.H and counts.

SourceFiles
    (header only)

\*---------------------------------------------------------------------------*/

#ifndef davofKromerBothe_H
#define davofKromerBothe_H

#include "davofCurvedGeometry.H"
#include <cmath>

namespace Foam
{
namespace davof
{

// * * * * * * * * * * * * * * Debug knobs * * * * * * * * * * * * * * * * //

//- The least number of sub-arcs every trace arc is split into (1: only the
//  quarter-turn rule of the ellipse); a convergence check of the arc
//  quadrature (key kromerBotheMinSplit)
inline label& kbMinSplit()
{
    static label v = 1;
    return v;
}

//- Whether the top and bottom faces of a prism and the faces of a cell are
//  integrated as fan triangles (true, the default; warped faces) or as whole
//  polygons (planar faces only; a consistency check of the fan)
inline bool& kbFanFacesOn()
{
    static bool v = true;
    return v;
}


// * * * * * * * * * * * * * * The conic trace  * * * * * * * * * * * * * * //

//- A conic Q11 x^2 + 2 Q12 x y + Q22 y^2 + 2 L1 x + 2 L2 y + K = 0 in the
//  plane coordinates xi (scaled by Ls inside), in its principal form after
//  the rotation xi = R u, R = [[cs, -sn], [sn, cs]]:
//  B1 u1^2 + B2 u2^2 + 2 l1 u1 + 2 l2 u2 + K = 0 (Kromer and Bothe 2023,
//  Eqs. 26-34). Used to order the boundary crossings along the conic and to
//  place arc split points, not to integrate.
struct conic2D
{
    enum Kind { EMPTY = 0, ELLIPSE, HYPERBOLA, PARABOLA, LINE, DEGENERATE };

    Kind kind = DEGENERATE;
    scalar Ls = 1;
    scalar cs = 1, sn = 0;
    scalar B1 = 0, B2 = 0, l1 = 0, l2 = 0, K = 0;
    scalar uc1 = 0, uc2 = 0, b = 0, r1 = 0, r2 = 0;
    bool swapped = false;
    scalar p01 = 0, p02 = 0, d1 = 0, d2 = 0;

    // the quadric's coefficients in the unscaled xi (for q and grad q)
    scalar Q11 = 0, Q12 = 0, Q22 = 0, L1 = 0, L2 = 0, K0 = 0;

    scalar q(const scalar x, const scalar y) const
    {
        return Q11*x*x + 2*Q12*x*y + Q22*y*y + 2*L1*x + 2*L2*y + K0;
    }
    void gradq(const scalar x, const scalar y, scalar& gx, scalar& gy) const
    {
        gx = 2*(Q11*x + Q12*y + L1);
        gy = 2*(Q12*x + Q22*y + L2);
    }

    void toU(scalar x1, scalar x2, scalar& u1, scalar& u2) const
    {
        x1 /= Ls; x2 /= Ls;
        u1 = cs*x1 + sn*x2;
        u2 = -sn*x1 + cs*x2;
    }
    void toXi(const scalar u1, const scalar u2, scalar& x1, scalar& x2) const
    {
        x1 = Ls*(cs*u1 - sn*u2);
        x2 = Ls*(sn*u1 + cs*u2);
    }

    label nBranches() const { return kind == HYPERBOLA ? 2 : 1; }
    bool cyclic() const { return kind == ELLIPSE; }

    //- The point of branch br at the parameter s (false: no such point)
    bool at(const label br, const scalar s, scalar& x1, scalar& x2) const
    {
        scalar u1 = 0, u2 = 0;
        switch (kind)
        {
            case ELLIPSE:
                u1 = uc1 + r1*std::cos(s);
                u2 = uc2 + r2*std::sin(s);
                break;
            case HYPERBOLA:
            {
                const scalar sg = (br == 0) ? 1 : -1;
                if (!swapped)
                {
                    const scalar g = (b - B2*s*s)/B1;
                    if (g < 0) return false;
                    u1 = uc1 + sg*std::sqrt(g); u2 = uc2 + s;
                }
                else
                {
                    const scalar g = (b - B1*s*s)/B2;
                    if (g < 0) return false;
                    u2 = uc2 + sg*std::sqrt(g); u1 = uc1 + s;
                }
                break;
            }
            case PARABOLA:
                if (!swapped) { u2 = s; u1 = -(B2*s*s + 2*l2*s + K)/(2*l1); }
                else          { u1 = s; u2 = -(B1*s*s + 2*l1*s + K)/(2*l2); }
                break;
            case LINE:
                u1 = p01 + s*d1; u2 = p02 + s*d2;
                break;
            default:
                return false;
        }
        toXi(u1, u2, x1, x2);
        // polish onto the conic along its gradient (two Newton steps): the
        // principal form loses digits for eccentric conics
        for (label it = 0; it < 2; ++it)
        {
            scalar gx, gy;
            gradq(x1, x2, gx, gy);
            const scalar g2 = gx*gx + gy*gy;
            if (!(g2 > 0)) break;
            const scalar f = q(x1, x2)/g2;
            x1 -= f*gx; x2 -= f*gy;
        }
        return true;
    }

    //- The (branch, parameter) of a point of the conic
    void paramOf(const scalar x1, const scalar x2, label& br, scalar& s) const
    {
        scalar u1, u2;
        toU(x1, x2, u1, u2);
        br = 0;
        switch (kind)
        {
            case ELLIPSE:
                s = std::atan2((u2 - uc2)/r2, (u1 - uc1)/r1);
                break;
            case HYPERBOLA:
                if (!swapped) { br = (u1 - uc1 >= 0) ? 0 : 1; s = u2 - uc2; }
                else          { br = (u2 - uc2 >= 0) ? 0 : 1; s = u1 - uc1; }
                break;
            case PARABOLA:
                s = swapped ? u1 : u2;
                break;
            case LINE:
                s = (u1 - p01)*d1 + (u2 - p02)*d2;
                break;
            default:
                s = 0;
        }
    }
};


//- The conic of the quadric q restricted to the plane: q(xi) = phi(x0 + xi1
//  f1 + xi2 f2) from six exact evaluations (phi quadratic), classified with
//  the coefficients scaled to a face of size Ls.
template<class PhiAt>
inline conic2D faceConic
(
    const point& x0, const vector& f1, const vector& f2, const scalar Ls,
    PhiAt&& phi
)
{
    conic2D C;
    C.Ls = Ls;
    const scalar s = Ls;
    const scalar q0 = phi(x0);
    const scalar qp1 = phi(x0 + s*f1), qm1 = phi(x0 - s*f1);
    const scalar qp2 = phi(x0 + s*f2), qm2 = phi(x0 - s*f2);
    const scalar q11 = phi(x0 + s*f1 + s*f2);
    // in the scaled coordinates xi' = xi/Ls
    const scalar A11 = 0.5*(qp1 + qm1 - 2*q0);
    const scalar A22 = 0.5*(qp2 + qm2 - 2*q0);
    const scalar a1 = 0.25*(qp1 - qm1);
    const scalar a2 = 0.25*(qp2 - qm2);
    const scalar A12 = 0.5*(q11 - A11 - A22 - 2*a1 - 2*a2 - q0);
    C.Q11 = A11/(Ls*Ls); C.Q12 = A12/(Ls*Ls); C.Q22 = A22/(Ls*Ls);
    C.L1 = a1/Ls; C.L2 = a2/Ls; C.K0 = q0;

    // principal axes
    scalar psi = 0;
    if (mag(A12) > 0)
    {
        psi = 0.5*std::atan2(2*A12, A11 - A22);
    }
    C.cs = std::cos(psi); C.sn = std::sin(psi);
    const scalar cs = C.cs, sn = C.sn;
    C.B1 = A11*cs*cs + 2*A12*cs*sn + A22*sn*sn;
    C.B2 = A11*sn*sn - 2*A12*cs*sn + A22*cs*cs;
    C.l1 = a1*cs + a2*sn;
    C.l2 = -a1*sn + a2*cs;
    C.K = q0;

    const scalar scaleAll = max(max(max(mag(C.B1), mag(C.B2)), max(mag(C.l1), mag(C.l2))), mag(C.K));
    if (!(scaleAll > 0))
    {
        C.kind = conic2D::DEGENERATE;
        return C;
    }
    const scalar epsB = 1e-13*scaleAll;
    const bool z1 = mag(C.B1) <= epsB, z2 = mag(C.B2) <= epsB;
    if (z1 && z2)
    {
        const scalar l2n = C.l1*C.l1 + C.l2*C.l2;
        if (l2n > sqr(1e-13*scaleAll))
        {
            const scalar ln = std::sqrt(l2n);
            C.kind = conic2D::LINE;
            C.d1 = -C.l2/ln; C.d2 = C.l1/ln;
            C.p01 = -C.K*C.l1/(2*l2n); C.p02 = -C.K*C.l2/(2*l2n);
        }
        else
        {
            C.kind = conic2D::DEGENERATE;
        }
        return C;
    }
    if (z1)
    {
        if (mag(C.l1) > epsB) { C.kind = conic2D::PARABOLA; C.swapped = false; }
        else C.kind = conic2D::DEGENERATE;
        return C;
    }
    if (z2)
    {
        if (mag(C.l2) > epsB) { C.kind = conic2D::PARABOLA; C.swapped = true; }
        else C.kind = conic2D::DEGENERATE;
        return C;
    }
    C.uc1 = -C.l1/C.B1; C.uc2 = -C.l2/C.B2;
    C.b = C.l1*C.l1/C.B1 + C.l2*C.l2/C.B2 - C.K;
    if (C.B1*C.B2 > 0)
    {
        if (!(C.b/C.B1 > 0))
        {
            C.kind = conic2D::EMPTY;
            return C;
        }
        C.kind = conic2D::ELLIPSE;
        C.r1 = std::sqrt(C.b/C.B1); C.r2 = std::sqrt(C.b/C.B2);
        return C;
    }
    if (mag(C.b) <= 1e-13*scaleAll)
    {
        C.kind = conic2D::DEGENERATE;
        return C;
    }
    C.kind = conic2D::HYPERBOLA;
    C.swapped = !(C.b/C.B1 > 0);
    return C;
}


// * * * * * * * * * * * * The polygon under the surface  * * * * * * * * * //

//- One convex planar polygon face of a polyhedron, counterclockwise about
//  its outward unit normal n, with the in-plane frame (x0; f1, f2),
//  f1 x f2 = n, and the size Ls.
struct kbFace
{
    List<point> p;
    vector n;
    point x0;
    vector f1, f2;
    scalar Ls;
    List<scalar> xi1, xi2;   // the vertices in plane coordinates

    void setFrame()
    {
        x0 = Zero;
        forAll(p, i) x0 += p[i];
        x0 /= p.size();
        n = Zero;
        for (label i = 1; i + 1 < p.size(); ++i)
        {
            n += 0.5*((p[i] - p[0]) ^ (p[i + 1] - p[0]));
        }
        const scalar A = mag(n);
        n = (A > 0) ? n/A : vector(0, 0, 1);
        Ls = max(std::sqrt(A), SMALL);
        vector d = p[0] - x0;
        d -= (d & n)*n;
        if (mag(d) < 1e-12*Ls) { d = p[1] - x0; d -= (d & n)*n; }
        f1 = d/max(mag(d), VSMALL);
        f2 = n ^ f1;
        xi1.setSize(p.size()); xi2.setSize(p.size());
        forAll(p, i)
        {
            xi1[i] = (p[i] - x0) & f1;
            xi2[i] = (p[i] - x0) & f2;
        }
    }

    //- Reverse the orientation (the normal flips)
    void reverse()
    {
        for (label i = 0, j = p.size() - 1; i < j; ++i, --j) Swap(p[i], p[j]);
    }

    point at(const scalar a, const scalar b) const { return x0 + a*f1 + b*f2; }

    //- Inside the convex polygon (counterclockwise), with a tolerance
    bool inside(const scalar a, const scalar b, const scalar tol) const
    {
        forAll(p, i)
        {
            const label j = (i + 1) % p.size();
            const scalar cr = (xi1[j] - xi1[i])*(b - xi2[i]) - (xi2[j] - xi2[i])*(a - xi1[i]);
            if (cr < -tol) return false;
        }
        return true;
    }

    scalar area() const
    {
        scalar A = 0;
        forAll(p, i)
        {
            const label j = (i + 1) % p.size();
            A += 0.5*(xi1[i]*xi2[j] - xi2[i]*xi1[j]);
        }
        return A;
    }
};


//- The accumulators of one polyhedron: the face sum of the divergence
//  theorem and the line integral of H1 around the projected patch; the
//  wetted area of the last face.
struct kbSums
{
    scalar faceSum = 0;      // sum_k ((x_k - o).n_k) A_k^wet
    scalar lineSum = 0;      // oint H1 dt2
    scalar lastArea = 0;     // A^wet of the last face integrated
    label nArcs = 0;
    label nLeaves = 0;
};


//- The reasons for which the rule leaves a face or a polyhedron unresolved
//  (and the volume rule falls back), counted per step and printed by the
//  rule's report(): a degenerate or empty conic with crossings, a conic
//  point that does not exist, a chord-normal root that does not exist, a
//  vanishing gradient on the arc, a wetted area outside [0, A], a volume
//  outside [0, V].
enum kbFail { KB_CONIC = 0, KB_AT, KB_ROOT, KB_GRAD, KB_AREA, KB_VOLUME, KB_NFAIL };
inline label& kbFailCount(const label r)
{
    static label c[KB_NFAIL] = {};
    return c[r];
}


//- The trace arcs and the wetted area of one face under the surface S;
//  accumulates into sums. Returns false where the configuration is not
//  resolved.
inline bool kbIntegrateFace
(
    const kbFace& F,
    const cellSurface& S,
    const scalar tolArea,      // absolute tolerance of the area integrals
    const scalar tolLine,      // absolute tolerance of the H1 line integrals
    const bool withLine,       // accumulate the H1 line integral (volume) too
    kbSums& sums
)
{
    auto phi = [&](const point& x) { return S.value(x); };
    const label nV = F.p.size();
    const scalar A = F.area();
    const scalar Ls = F.Ls;
    sums.lastArea = 0;

    // the vertex values and the boundary crossings: the exact roots along
    // the edges and the vertices on the surface
    List<scalar> fv(nV);
    scalar fmax = 0;
    forAll(F.p, i) { fv[i] = phi(F.p[i]); fmax = max(fmax, mag(fv[i])); }
    const scalar epsV = 1e-14*max(fmax, SMALL);

    struct crossing { scalar a, b; label br; scalar s; };
    DynamicList<crossing> cr(8);
    DynamicList<scalar> roots(4);
    scalar edgeArea = 0;   // the liquid pieces of the edges, 1/2 (P x Q)
    bool anyCrossing = false;
    forAll(F.p, i)
    {
        const label j = (i + 1) % nV;
        roots.clear();
        segmentRoots(F.p[i], F.p[j], phi, roots);
        for (label a = 1; a < roots.size(); ++a)
        {
            for (label k = a; k > 0 && roots[k] < roots[k - 1]; --k) Swap(roots[k], roots[k - 1]);
        }
        const scalar dx = F.xi1[j] - F.xi1[i], dy = F.xi2[j] - F.xi2[i];
        scalar r0 = 0;
        for (label k = 0; k <= roots.size(); ++k)
        {
            const scalar r1 = (k < roots.size()) ? roots[k] : scalar(1);
            if (r1 > r0)
            {
                const scalar rm = 0.5*(r0 + r1);
                const point pm = F.p[i] + rm*(F.p[j] - F.p[i]);
                if (phi(pm) < 0)
                {
                    const scalar ax = F.xi1[i] + r0*dx, ay = F.xi2[i] + r0*dy;
                    const scalar bx = F.xi1[i] + r1*dx, by = F.xi2[i] + r1*dy;
                    edgeArea += 0.5*(ax*by - ay*bx);
                }
            }
            if (k < roots.size())
            {
                crossing c;
                c.a = F.xi1[i] + roots[k]*dx; c.b = F.xi2[i] + roots[k]*dy;
                c.br = 0; c.s = 0;
                cr.append(c);
                anyCrossing = true;
            }
            r0 = r1;
        }
        if (mag(fv[i]) <= epsV)
        {
            crossing c;
            c.a = F.xi1[i]; c.b = F.xi2[i]; c.br = 0; c.s = 0;
            cr.append(c);
            anyCrossing = true;
        }
    }

    const conic2D C = faceConic(F.x0, F.f1, F.f2, Ls, phi);

    // the arcs inside the polygon: (branch, s_a, s_b, endpoints)
    struct arc { label br; scalar sa, sb; scalar ax, ay, bx, by; bool exactEnds; };
    DynamicList<arc> arcs(8);

    if (!anyCrossing)
    {
        // no crossing: the polygon is on one side, unless an ellipse lies
        // fully inside it
        bool enclosed = false;
        if (C.kind == conic2D::ELLIPSE)
        {
            scalar x1, x2;
            if (C.at(0, 0, x1, x2) && F.inside(x1, x2, -1e-12*Ls))
            {
                enclosed = true;
            }
        }
        if (!enclosed)
        {
            const point pc = F.at(0, 0);
            scalar fc = phi(pc);
            // the centroid value decides; a centroid on the surface without
            // crossings is a tangency: use the vertices
            if (mag(fc) <= epsV) fc = fv[0];
            sums.lastArea = (fc < 0) ? A : 0;
            sums.faceSum += ((F.x0 - S.o) & F.n)*sums.lastArea;
            return true;
        }
        for (label k = 0; k < 4; ++k)
        {
            arc ar;
            ar.br = 0; ar.sa = k*0.5*constant::mathematical::pi; ar.sb = (k + 1)*0.5*constant::mathematical::pi;
            ar.exactEnds = false;
            arcs.append(ar);
        }
    }
    else
    {
        if (C.kind == conic2D::DEGENERATE || C.kind == conic2D::EMPTY)
        {
            ++kbFailCount(KB_CONIC);
            return false;
        }
        forAll(cr, i) C.paramOf(cr[i].a, cr[i].b, cr[i].br, cr[i].s);
        // sort by (branch, s)
        for (label i = 1; i < cr.size(); ++i)
        {
            for (label k = i; k > 0 && (cr[k].br < cr[k - 1].br || (cr[k].br == cr[k - 1].br && cr[k].s < cr[k - 1].s)); --k)
            {
                Swap(cr[k], cr[k - 1]);
            }
        }
        const scalar twoPi = 2*constant::mathematical::pi;
        for (label br = 0; br < C.nBranches(); ++br)
        {
            label i0 = -1, i1 = -1;
            forAll(cr, i)
            {
                if (cr[i].br == br) { if (i0 < 0) i0 = i; i1 = i; }
            }
            if (i0 < 0) continue;
            const label n = i1 - i0 + 1;
            if (C.cyclic())
            {
                if (n == 1)
                {
                    // a single touching point: the whole ellipse if inside
                    scalar x1, x2;
                    if (C.at(0, cr[i0].s + constant::mathematical::pi, x1, x2) && F.inside(x1, x2, -1e-12*Ls))
                    {
                        for (label k = 0; k < 4; ++k)
                        {
                            arc ar;
                            ar.br = 0; ar.sa = cr[i0].s + k*0.25*twoPi; ar.sb = cr[i0].s + (k + 1)*0.25*twoPi;
                            ar.exactEnds = false;
                            arcs.append(ar);
                        }
                    }
                    continue;
                }
                for (label k = 0; k < n; ++k)
                {
                    const crossing& ca = cr[i0 + k];
                    const crossing& cb = cr[i0 + (k + 1) % n];
                    scalar ds = cb.s - ca.s;
                    if (k == n - 1) ds += twoPi;
                    if (ds < 0) ds += twoPi;
                    arc ar;
                    ar.br = 0; ar.sa = ca.s; ar.sb = ca.s + ds;
                    ar.ax = ca.a; ar.ay = ca.b; ar.bx = cb.a; ar.by = cb.b;
                    ar.exactEnds = true;
                    arcs.append(ar);
                }
            }
            else
            {
                for (label k = 0; k + 1 < n; ++k)
                {
                    const crossing& ca = cr[i0 + k];
                    const crossing& cb = cr[i0 + k + 1];
                    arc ar;
                    ar.br = br; ar.sa = ca.s; ar.sb = cb.s;
                    ar.ax = ca.a; ar.ay = ca.b; ar.bx = cb.a; ar.by = cb.b;
                    ar.exactEnds = true;
                    arcs.append(ar);
                }
            }
        }
    }

    // the arcs inside the polygon, split to at most a quarter turn
    const List<Pair<scalar>>& g = gl12();
    scalar arcArea = 0, arcLine = 0;
    forAll(arcs, ia)
    {
        const arc& ar = arcs[ia];
        if (!(ar.sb > ar.sa)) continue;
        // inside?
        {
            scalar mx, my;
            if (!C.at(ar.br, 0.5*(ar.sa + ar.sb), mx, my)) continue;
            if (!F.inside(mx, my, 1e-12*Ls)) continue;
        }
        label m = max(label(1), kbMinSplit());
        if (C.kind == conic2D::ELLIPSE)
        {
            m = max(m, label(std::ceil((ar.sb - ar.sa)/(0.5*constant::mathematical::pi) - 1e-9)));
        }
        for (label k = 0; k < m; ++k)
        {
            const scalar s0 = ar.sa + (ar.sb - ar.sa)*k/m;
            const scalar s1 = ar.sa + (ar.sb - ar.sa)*(k + 1)/m;
            scalar ax, ay, bx, by;
            if (k == 0 && ar.exactEnds) { ax = ar.ax; ay = ar.ay; }
            else if (!C.at(ar.br, s0, ax, ay)) { ++kbFailCount(KB_AT); return false; }
            if (k == m - 1 && ar.exactEnds) { bx = ar.bx; by = ar.by; }
            else if (!C.at(ar.br, s1, bx, by)) { ++kbFailCount(KB_AT); return false; }
            // the side of the chord the arc lies on, from its midpoint
            scalar mx, my;
            if (!C.at(ar.br, 0.5*(s0 + s1), mx, my)) { ++kbFailCount(KB_AT); return false; }

            const scalar ex = bx - ax, ey = by - ay;
            const scalar Le = std::sqrt(ex*ex + ey*ey);
            if (Le <= 1e-14*Ls) continue;
            const scalar nx = -ey/Le, ny = ex/Le;      // left normal of the chord
            const scalar rM = (mx - ax)*nx + (my - ay)*ny;
            const scalar sigma = (rM >= 0) ? 1 : -1;
            // a flat arc (a line trace, or a parabola of round-off curvature
            // on a plane) has a midpoint on its chord to round-off: its side
            // is not known, and the root along the chord's normal is the one
            // nearest to the chord; selecting it by the side took the far
            // root of the near-degenerate quadratic (|r| ~ 1e9 Ls) in up to
            // 70000 prisms per step on the translated plane (2026-10-25)
            const bool sideKnown = mag(rM) > 1e-9*Ls;

            // the point of the arc over the chord parameter l and its derivative
            auto onArc = [&](const scalar l, scalar& x, scalar& y, scalar& dx, scalar& dy) -> bool
            {
                const scalar cx = ax + l*ex, cy = ay + l*ey;
                const scalar al = C.Q11*nx*nx + 2*C.Q12*nx*ny + C.Q22*ny*ny;
                const scalar be = 2*((C.Q11*cx + C.Q12*cy + C.L1)*nx + (C.Q12*cx + C.Q22*cy + C.L2)*ny);
                const scalar ga = C.q(cx, cy);
                scalar r;
                if (mag(al) <= 1e-14*(mag(be)/Ls + mag(ga)/(Ls*Ls)))
                {
                    if (!(mag(be) > 0)) { ++kbFailCount(KB_ROOT); return false; }
                    r = -ga/be;
                }
                else
                {
                    const scalar disc = be*be - 4*al*ga;
                    if (disc < 0) { ++kbFailCount(KB_ROOT); return false; }
                    const scalar sq = std::sqrt(disc);
                    const scalar qq = -0.5*(be + ((be >= 0) ? sq : -sq));
                    const scalar ra = (qq != 0) ? ga/qq : -be/(2*al);
                    const scalar rb = (al != 0) ? qq/al : ra;
                    const bool sa = (ra*sigma >= 0), sb = (rb*sigma >= 0);
                    if (!sideKnown) r = (mag(ra) <= mag(rb)) ? ra : rb;
                    else if (sa && sb) r = (mag(ra) <= mag(rb)) ? ra : rb;
                    else if (sa) r = ra;
                    else if (sb) r = rb;
                    else r = (mag(ra) <= mag(rb)) ? ra : rb;
                }
                x = cx + r*nx; y = cy + r*ny;
                scalar gx, gy;
                C.gradq(x, y, gx, gy);
                const scalar gn = gx*nx + gy*ny;
                if (!(mag(gn) > 0)) { ++kbFailCount(KB_GRAD); return false; }
                const scalar rp = -(gx*ex + gy*ey)/gn;
                dx = ex + rp*nx; dy = ey + rp*ny;
                return true;
            };

            // orientation of the arc for the area (liquid on the left) and
            // for the line integral (the projected patch on the left), from
            // the midpoint
            scalar mX, mY, mdx, mdy;
            if (!onArc(0.5, mX, mY, mdx, mdy)) return false;
            scalar gx, gy;
            C.gradq(mX, mY, gx, gy);
            const scalar leftX = -mdy, leftY = mdx;
            const scalar sArea = ((leftX*(-gx) + leftY*(-gy)) >= 0) ? 1 : -1;
            scalar sLine = 1;
            if (withLine)
            {
                const point X = F.at(mX, mY);
                const vector T = mdx*F.f1 + mdy*F.f2;
                const scalar t1 = (X - S.o) & S.e1, t2 = (X - S.o) & S.e2;
                const scalar h1 = S.a1 + S.b11*t1 + S.b12*t2;
                const scalar h2 = S.a2 + S.b12*t1 + S.b22*t2;
                const scalar n1 = F.n & S.e1, n2 = F.n & S.e2, n3 = F.n & S.e3;
                const scalar m1 = -(n1 + h1*n3), m2 = -(n2 + h2*n3);
                const scalar tau1 = T & S.e1, tau2 = T & S.e2;
                const scalar nu1 = -tau2, nu2 = tau1;
                sLine = ((nu1*m1 + nu2*m2) >= 0) ? 1 : -1;
            }

            // the integrands over l in [0, 1]
            auto integrands = [&](const scalar l, scalar& fA, scalar& fL) -> bool
            {
                scalar x, y, dx, dy;
                if (!onArc(l, x, y, dx, dy)) return false;
                fA = 0.5*(x*dy - y*dx);
                fL = 0;
                if (withLine)
                {
                    const point X = F.at(x, y);
                    const scalar t1 = (X - S.o) & S.e1, t2 = (X - S.o) & S.e2;
                    const scalar H1 = S.c*t1 - (S.b11*t1*t1*t1/6 + 0.5*S.b12*t1*t1*t2 + 0.5*S.b22*t1*t2*t2);
                    const scalar dt2 = dx*(F.f1 & S.e2) + dy*(F.f2 & S.e2);
                    fL = H1*dt2;
                }
                return true;
            };
            auto gaussL = [&](const scalar l0, const scalar l1, scalar& IA, scalar& IL) -> bool
            {
                IA = 0; IL = 0;
                for (const Pair<scalar>& gw : g)
                {
                    scalar fA, fL;
                    if (!integrands(l0 + gw.first()*(l1 - l0), fA, fL)) return false;
                    IA += gw.second()*(l1 - l0)*fA;
                    IL += gw.second()*(l1 - l0)*fL;
                }
                return true;
            };
            bool ok = true;
            std::function<void(scalar, scalar, scalar, scalar, label, scalar&, scalar&)> rec =
                [&](scalar l0, scalar l1, scalar pA, scalar pL, label depth, scalar& oA, scalar& oL)
                {
                    const scalar lm = 0.5*(l0 + l1);
                    scalar a0, L0, a1, L1;
                    if (!gaussL(l0, lm, a0, L0) || !gaussL(lm, l1, a1, L1)) { ok = false; oA = pA; oL = pL; return; }
                    ++sums.nLeaves;
                    const scalar cA = a0 + a1, cL = L0 + L1;
                    if ((mag(cA - pA) <= tolArea && mag(cL - pL) <= tolLine) || depth >= 20)
                    {
                        oA = cA; oL = cL; return;
                    }
                    scalar rA0, rL0, rA1, rL1;
                    rec(l0, lm, a0, L0, depth + 1, rA0, rL0);
                    rec(lm, l1, a1, L1, depth + 1, rA1, rL1);
                    oA = rA0 + rA1; oL = rL0 + rL1;
                };
            scalar pA, pL;
            if (!gaussL(0, 1, pA, pL)) return false;
            scalar IA, IL;
            rec(0, 1, pA, pL, 0, IA, IL);
            if (!ok) return false;
            arcArea += sArea*IA;
            arcLine += sLine*IL;
            ++sums.nArcs;
        }
    }

    const scalar Awet = edgeArea + arcArea;
    if (Awet < -1e-9*A || Awet > (1 + 1e-9)*A)
    {
        ++kbFailCount(KB_AREA);
        // the first faces of a run that fail the area check, for the
        // diagnosis of the rule on flat traces (2026-10-25)
        static label nPrinted = 0;
        if (nPrinted < 40)
        {
            ++nPrinted;
            Pout<< "kromerBothe area check failed: kind " << label(C.kind)
                << " nV " << nV << " A " << A << " Ls " << Ls
                << " Awet/A " << Awet/A << " edge/A " << edgeArea/A
                << " arc/A " << arcArea/A << " crossings " << cr.size()
                << " arcs " << arcs.size() << " nArcs " << sums.nArcs
                << " scaled B1 " << C.B1 << " B2 " << C.B2 << " l1 " << C.l1
                << " l2 " << C.l2 << " K " << C.K
                << " surface b11 " << S.b11 << " b12 " << S.b12 << " b22 " << S.b22
                << " c " << S.c << " a1 " << S.a1 << " a2 " << S.a2
                << " fv " << fv << nl;
        }
        return false;
    }
    sums.lastArea = min(max(Awet, scalar(0)), A);
    sums.faceSum += ((F.x0 - S.o) & F.n)*sums.lastArea;
    if (withLine) sums.lineSum += arcLine;
    return true;
}


//- The fan triangles of a polygon about its centroid, oriented so that
//  their normals point away from the interior point xin, appended to faces.
inline void kbFanFaces(const List<point>& poly, const point& xin, DynamicList<kbFace>& faces)
{
    if (!kbFanFacesOn())
    {
        kbFace F;
        F.p = poly;
        F.setFrame();
        if (((F.x0 - xin) & F.n) < 0)
        {
            F.reverse();
            F.setFrame();
        }
        if (mag(F.n) > 0 && F.Ls > SMALL) faces.append(F);
        return;
    }
    point xc = Zero;
    forAll(poly, i) xc += poly[i];
    xc /= poly.size();
    forAll(poly, i)
    {
        const label j = (i + 1) % poly.size();
        kbFace F;
        F.p.setSize(3);
        F.p[0] = xc; F.p[1] = poly[i]; F.p[2] = poly[j];
        F.setFrame();
        if (((F.x0 - xin) & F.n) < 0)
        {
            F.reverse();
            F.setFrame();
        }
        if (mag(F.n) > 0 && F.Ls > SMALL) faces.append(F);
    }
}


//- The liquid volume of the polyhedron given by its faces (outward
//  oriented) under the paraboloid S. Returns false where a face is not
//  resolved.
inline bool kbPolyhedronLiquidVolume
(
    const DynamicList<kbFace>& faces,
    const cellSurface& S,
    const scalar tolAbs,         // of the volume
    scalar& vLiq,
    label& nLeaves
)
{
    kbSums sums;
    scalar Ltot = 0;
    forAll(faces, k) Ltot = max(Ltot, faces[k].Ls);
    // the face terms are multiplied by |x - o| ~ Ltot and the H1 terms are
    // volumes already: tolerances of the area and the line integrals
    const scalar tolArea = tolAbs/max(Ltot, SMALL)/max(faces.size(), label(1));
    const scalar tolLine = tolAbs/max(faces.size(), label(1));
    forAll(faces, k)
    {
        if (!kbIntegrateFace(faces[k], S, tolArea, tolLine, true, sums))
        {
            nLeaves += sums.nLeaves;
            return false;
        }
    }
    nLeaves += sums.nLeaves;
    vLiq = (sums.faceSum + sums.lineSum)/3.0;
    return true;
}


//- The faces of cellI as fan triangles, outward (the owner's face order,
//  reversed for the neighbour).
inline void kbCellFaces(const fvMesh& mesh, const label cellI, DynamicList<kbFace>& faces)
{
    const pointField& P = mesh.points();
    const faceList& fcs = mesh.faces();
    const labelList& own = mesh.faceOwner();
    const vectorField& Cf = mesh.faceCentres();
    const labelList& cFaces = mesh.cells()[cellI];
    faces.clear();
    forAll(cFaces, cf)
    {
        const label faceI = cFaces[cf];
        const face& f = fcs[faceI];
        const bool isOwner = (own[faceI] == cellI);
        const point& xf = Cf[faceI];
        if (!kbFanFacesOn())
        {
            kbFace F;
            F.p.setSize(f.size());
            forAll(f, ip) F.p[ip] = isOwner ? P[f[ip]] : P[f[f.size() - 1 - ip]];
            F.setFrame();
            if (F.Ls > SMALL) faces.append(F);
            continue;
        }
        forAll(f, ip)
        {
            const label ip1 = f.nextLabel(ip);
            kbFace F;
            F.p.setSize(3);
            F.p[0] = xf;
            if (isOwner) { F.p[1] = P[f[ip]]; F.p[2] = P[ip1]; }
            else         { F.p[1] = P[ip1];   F.p[2] = P[f[ip]]; }
            F.setFrame();
            if (F.Ls > SMALL) faces.append(F);
        }
    }
}


//- The faces of the prism between the face faceI displaced by -d and the
//  face: the fan triangles of both and the side parallelograms, outward.
inline void kbPrismFaces(const fvMesh& mesh, const label faceI, const vector& d, DynamicList<kbFace>& faces)
{
    const pointField& P = mesh.points();
    const face& f = mesh.faces()[faceI];
    const point xin = mesh.faceCentres()[faceI] - 0.5*d;
    faces.clear();
    List<point> top(f.size()), bottom(f.size());
    forAll(f, i) { top[i] = P[f[i]]; bottom[i] = P[f[i]] - d; }
    kbFanFaces(top, xin, faces);
    kbFanFaces(bottom, xin, faces);
    forAll(f, i)
    {
        const label j = f.fcIndex(i);
        kbFace F;
        F.p.setSize(4);
        F.p[0] = top[i]; F.p[1] = top[j]; F.p[2] = bottom[j]; F.p[3] = bottom[i];
        F.setFrame();
        if (((F.x0 - xin) & F.n) < 0) { F.reverse(); F.setFrame(); }
        if (F.Ls > SMALL) faces.append(F);
    }
}


//- The liquid area fraction of the face faceI displaced by d under S by the
//  Green integral (fan triangles). Returns false where not resolved.
inline bool kbFaceLiquidFraction
(
    const fvMesh& mesh,
    const label faceI,
    const vector& d,
    const cellSurface& S,
    const scalar tolRel,
    scalar& alphaf,
    label& nLeaves
)
{
    const pointField& P = mesh.points();
    const face& f = mesh.faces()[faceI];
    const vector& Sf = mesh.faceAreas()[faceI];
    const scalar Af = mag(Sf);
    if (!(Af > 0)) { alphaf = 0; return true; }
    const point xf = mesh.faceCentres()[faceI] + d;
    kbSums sums;
    scalar Awet = 0, Atot = 0;
    forAll(f, ip)
    {
        const label ip1 = f.nextLabel(ip);
        kbFace F;
        F.p.setSize(3);
        F.p[0] = xf; F.p[1] = P[f[ip]] + d; F.p[2] = P[ip1] + d;
        F.setFrame();
        if (!(F.Ls > SMALL)) continue;
        const scalar At = F.area();
        if (!kbIntegrateFace(F, S, tolRel*Af, GREAT, false, sums))
        {
            nLeaves += sums.nLeaves;
            return false;
        }
        Awet += sums.lastArea;
        Atot += At;
    }
    nLeaves += sums.nLeaves;
    alphaf = (Atot > 0) ? Awet/Atot : 0;
    return true;
}

} // End namespace davof
} // End namespace Foam

#endif

// ************************************************************************* //
