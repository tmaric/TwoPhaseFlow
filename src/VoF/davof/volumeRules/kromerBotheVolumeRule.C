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

\*---------------------------------------------------------------------------*/

#include "kromerBotheVolumeRule.H"
#include "davofKromerBothe.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace davofVolumeRules
{
    defineTypeNameAndDebug(kromerBothe, 0);
    addToRunTimeSelectionTable(davofVolumeRule, kromerBothe, dictionary);
}
}


namespace
{
    // the self-test (defined below)
    void kbSelfTest(const Foam::scalar tolRel, const Foam::label maxDepth);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofVolumeRules::kromerBothe::kromerBothe
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    davofVolumeRule(mesh, dict),
    debug_(dict.getOrDefault<bool>("volumeRuleDebug", false)),
    nDebugPrinted_(0),
    dbgMaxPrism_(0), dbgSumPrism_(0), dbgMaxCell_(0), dbgSumCell_(0),
    dbgNPrism_(0), dbgNCell_(0)
{
    // the debug knobs of davofKromerBothe.H
    davof::kbMinSplit() = dict.getOrDefault<label>("kromerBotheMinSplit", 1);
    davof::kbFanFacesOn() = dict.getOrDefault<bool>("kromerBotheFan", true);
    if (debug_)
    {
        Info<< "kromerBothe: minSplit " << davof::kbMinSplit()
            << ", fan faces " << davof::kbFanFacesOn() << endl;
        kbSelfTest(tolRel_, maxDepth_);
    }
}


namespace
{

// The self-test of the rule (volumeRuleDebug): the liquid volume of an
// axis-aligned box under the axis-aligned paraboloid z = c + b (x^2 + y^2)/2
// (e3 = z, the liquid below), against a reference that is exact in x (the
// clamped quadratic integrated piecewise between its exact roots) and
// adaptive Gauss-Legendre in y, and against the split quadrature. Thin boxes
// in x reproduce the equatorial prisms of the sphere, whose traces on the
// side faces are parabolas and on the top and bottom circles.
void kbSelfTest(const Foam::scalar tolRel, const Foam::label maxDepth)
{
    using namespace Foam;
    const scalar h = 1.0/32, dx = 0.2*h, b = -7.0, c = 0.0;
    davof::cellSurface S;
    S.paraboloid = true;
    S.o = point::zero;
    S.e1 = vector(1, 0, 0); S.e2 = vector(0, 1, 0); S.e3 = vector(0, 0, 1);
    S.b11 = b; S.b22 = b; S.b12 = 0; S.c = c; S.a1 = 0; S.a2 = 0;
    auto phi = [&](const point& x) { return S.value(x); };
    const List<Pair<scalar>>& g = davof::gl12();

    // the reference
    auto inner = [&](const scalar y, const scalar x0, const scalar x1, const scalar z0) -> scalar
    {
        // clamp(gq(x), 0, h) with gq(x) = (c - z0 + b y^2/2) + b x^2/2
        const scalar k = c - z0 + 0.5*b*y*y;
        DynamicList<scalar> br(6);
        br.append(x0); br.append(x1);
        for (const scalar lvl : {scalar(0), h})
        {
            // b x^2/2 = lvl - k
            const scalar r2 = 2*(lvl - k)/b;
            if (r2 > 0)
            {
                const scalar r = std::sqrt(r2);
                if (r > x0 && r < x1) br.append(r);
                if (-r > x0 && -r < x1) br.append(-r);
            }
        }
        for (label i = 1; i < br.size(); ++i)
            for (label j = i; j > 0 && br[j] < br[j - 1]; --j) Swap(br[j], br[j - 1]);
        scalar sum = 0;
        for (label i = 0; i + 1 < br.size(); ++i)
        {
            const scalar xa = br[i], xb = br[i + 1];
            if (!(xb > xa)) continue;
            const scalar xm = 0.5*(xa + xb);
            const scalar gm = k + 0.5*b*xm*xm;
            if (gm >= h) sum += h*(xb - xa);
            else if (gm > 0) sum += k*(xb - xa) + b*(xb*xb*xb - xa*xa*xa)/6;
        }
        return sum;
    };
    auto reference = [&](const scalar x0, const scalar x1, const scalar y0, const scalar y1, const scalar z0) -> scalar
    {
        std::function<scalar(scalar, scalar, scalar, label)> rec =
            [&](scalar ya, scalar yb, scalar parent, label depth) -> scalar
            {
                auto gauss = [&](scalar a, scalar bq) {
                    scalar s = 0;
                    for (const Pair<scalar>& gw : g) s += gw.second()*(bq - a)*inner(a + gw.first()*(bq - a), x0, x1, z0);
                    return s;
                };
                const scalar ym = 0.5*(ya + yb);
                const scalar c0 = gauss(ya, ym), c1 = gauss(ym, yb);
                if ((depth >= 4 && mag(c0 + c1 - parent) <= 1e-15*h*h*h) || depth >= 40) return c0 + c1;
                return rec(ya, ym, c0, depth + 1) + rec(ym, yb, c1, depth + 1);
            };
        scalar p = 0;
        for (const Pair<scalar>& gw : g) p += gw.second()*(y1 - y0)*inner(y0 + gw.first()*(y1 - y0), x0, x1, z0);
        return rec(y0, y1, p, 0);
    };

    Info<< "kromerBothe self-test: thin boxes (dx = 0.2 h, h = 1/32) under z = "
        << c << " + " << b << " (x^2 + y^2)/2" << nl;
    label nBox = 0;
    for (const scalar x0 : {0.0, 0.02, 0.05, 0.1})
    {
        for (const scalar y0 : {0.0, 0.03, 0.1})
        {
            for (const scalar zShift : {0.5, 0.2, 0.8})
            {
                const scalar x1 = x0 + dx, y1 = y0 + h;
                // the surface height at the box centre, the box placed so that
                // the surface sits at zShift of its height there
                const scalar zc = c + 0.5*b*(sqr(0.5*(x0 + x1)) + sqr(0.5*(y0 + y1)));
                const scalar z0 = zc - zShift*h, z1 = z0 + h;
                const scalar pv = dx*h*h;
                const scalar vRef = reference(x0, x1, y0, y1, z0);
                // the box as a Kromer-Bothe polyhedron (fan faces about an
                // interior point)
                const point xin(0.5*(x0 + x1), 0.5*(y0 + y1), 0.5*(z0 + z1));
                DynamicList<davof::kbFace> faces(12);
                auto quad = [&](point a, point bq, point cq, point d) {
                    List<point> poly(4); poly[0] = a; poly[1] = bq; poly[2] = cq; poly[3] = d;
                    davof::kbFanFaces(poly, xin, faces);
                };
                quad(point(x0,y0,z0), point(x1,y0,z0), point(x1,y1,z0), point(x0,y1,z0));
                quad(point(x0,y0,z1), point(x1,y0,z1), point(x1,y1,z1), point(x0,y1,z1));
                quad(point(x0,y0,z0), point(x1,y0,z0), point(x1,y0,z1), point(x0,y0,z1));
                quad(point(x0,y1,z0), point(x1,y1,z0), point(x1,y1,z1), point(x0,y1,z1));
                quad(point(x0,y0,z0), point(x0,y1,z0), point(x0,y1,z1), point(x0,y0,z1));
                quad(point(x1,y0,z0), point(x1,y1,z0), point(x1,y1,z1), point(x1,y0,z1));
                scalar vKb = 0;
                label nLeavesKb = 0;
                const bool ok = davof::kbPolyhedronLiquidVolume(faces, S, tolRel*pv, vKb, nLeavesKb);
                // the split quadrature: the prism of the face x = x1 displaced by -d
                const vector d(dx, 0, 0);
                const point fc(x1, 0.5*(y0 + y1), 0.5*(z0 + z1));
                const point q0(x1,y0,z0), q1(x1,y1,z0), q2(x1,y1,z1), q3(x1,y0,z1);
                label nls = 0;
                scalar vSplit = 0;
                vSplit += davof::prismLiquidVolume(fc, q0, q1, d, phi, tolRel*pv, maxDepth, nls);
                vSplit += davof::prismLiquidVolume(fc, q1, q2, d, phi, tolRel*pv, maxDepth, nls);
                vSplit += davof::prismLiquidVolume(fc, q2, q3, d, phi, tolRel*pv, maxDepth, nls);
                vSplit += davof::prismLiquidVolume(fc, q3, q0, d, phi, tolRel*pv, maxDepth, nls);
                ++nBox;
                Info<< "  box " << nBox << " x0 " << x0 << " y0 " << y0 << " zShift " << zShift
                    << ": ref " << vRef/pv << " kb " << vKb/pv << (ok ? "" : " (NOT RESOLVED)")
                    << " split " << vSplit/pv
                    << "  (kb - ref)/pv " << (vKb - vRef)/pv
                    << "  (split - ref)/pv " << (vSplit - vRef)/pv << nl;
                if (mag(vKb - vRef) > 1e-9*pv)
                {
                    forAll(faces, k)
                    {
                        davof::kbSums sums;
                        davof::kbIntegrateFace(faces[k], S, 1e-16*pv, 1e-16*pv, true, sums);
                        Info<< "    face " << k << " n " << faces[k].n << " wet " << sums.lastArea
                            << " of " << faces[k].area() << " faceTerm " << sums.faceSum
                            << " line " << sums.lineSum << " arcs " << sums.nArcs << nl;
                    }
                }
            }
        }
    }

    // Exact per-face references for the axis-aligned frame (e3 = z, e1 = x,
    // e2 = y) with the general tensor: the liquid is z < hz(x, y) = c + (b11
    // x^2 + 2 b12 x y + b22 y^2)/2, a quadratic in x for fixed y and in y
    // for fixed x, so the clamped heights integrate exactly between their
    // roots along one direction and adaptively along the other
    struct AxisRef
    {
        scalar c, b11, b12, b22;
        scalar hz(scalar x, scalar y) const { return c + 0.5*(b11*x*x + 2*b12*x*y + b22*y*y); }
        // int_{xa}^{xb} clamp(hz(x, y) - zlo, 0, H) dx, exact
        scalar clampX(scalar y, scalar xa, scalar xb, scalar zlo, scalar H) const
        {
            // hz - zlo = A x^2 + B x + C
            const scalar A = 0.5*b11, B = b12*y, C = c + 0.5*b22*y*y - zlo;
            DynamicList<scalar> br(8);
            br.append(xa); br.append(xb);
            for (const scalar lvl : {scalar(0), H})
            {
                if (mag(A) > VSMALL)
                {
                    const scalar disc = B*B - 4*A*(C - lvl);
                    if (disc > 0)
                    {
                        const scalar sq = std::sqrt(disc);
                        for (const scalar r : {(-B - sq)/(2*A), (-B + sq)/(2*A)})
                            if (r > xa && r < xb) br.append(r);
                    }
                }
                else if (mag(B) > VSMALL)
                {
                    const scalar r = (lvl - C)/B;
                    if (r > xa && r < xb) br.append(r);
                }
            }
            for (label i = 1; i < br.size(); ++i)
                for (label j = i; j > 0 && br[j] < br[j - 1]; --j) Swap(br[j], br[j - 1]);
            scalar sum = 0;
            for (label i = 0; i + 1 < br.size(); ++i)
            {
                const scalar p = br[i], q = br[i + 1];
                if (!(q > p)) continue;
                const scalar m = 0.5*(p + q);
                const scalar gm = A*m*m + B*m + C;
                if (gm >= H) sum += H*(q - p);
                else if (gm > 0) sum += A*(q*q*q - p*p*p)/3 + B*(q*q - p*p)/2 + C*(q - p);
            }
            return sum;
        }
        // int_{xa}^{xb} 1[zlo < hz < zhi] w(x) dx with w(x) = c - (b11 x^2 + 2 b12 x y + b22 y^2)/2 (the patch integrand), exact
        scalar bandX(scalar y, scalar xa, scalar xb, scalar zlo, scalar zhi) const
        {
            const scalar A = 0.5*b11, B = b12*y, C = c + 0.5*b22*y*y;
            DynamicList<scalar> br(8);
            br.append(xa); br.append(xb);
            for (const scalar lvl : {zlo, zhi})
            {
                if (mag(A) > VSMALL)
                {
                    const scalar disc = B*B - 4*A*(C - lvl);
                    if (disc > 0)
                    {
                        const scalar sq = std::sqrt(disc);
                        for (const scalar r : {(-B - sq)/(2*A), (-B + sq)/(2*A)})
                            if (r > xa && r < xb) br.append(r);
                    }
                }
                else if (mag(B) > VSMALL)
                {
                    const scalar r = (lvl - C)/B;
                    if (r > xa && r < xb) br.append(r);
                }
            }
            for (label i = 1; i < br.size(); ++i)
                for (label j = i; j > 0 && br[j] < br[j - 1]; --j) Swap(br[j], br[j - 1]);
            scalar sum = 0;
            for (label i = 0; i + 1 < br.size(); ++i)
            {
                const scalar p = br[i], q = br[i + 1];
                if (!(q > p)) continue;
                const scalar m = 0.5*(p + q);
                const scalar hm = A*m*m + B*m + C;
                if (hm > zlo && hm < zhi)
                {
                    // w = c - A x^2 - B x - (C - c) ... w = c - 0.5 b11 x^2 - b12 x y - 0.5 b22 y^2
                    sum += c*(q - p) - A*(q*q*q - p*p*p)/3 - B*(q*q - p*p)/2 - 0.5*b22*y*y*(q - p);
                }
            }
            return sum;
        }
        // the adaptive outer integral over y of a function f(y)
        scalar outer(const std::function<scalar(scalar)>& f, scalar ya, scalar yb, scalar tol) const
        {
            const List<Pair<scalar>>& g = davof::gl12();
            auto gauss = [&](scalar a, scalar bq) {
                scalar s = 0;
                for (const Pair<scalar>& gw : g) s += gw.second()*(bq - a)*f(a + gw.first()*(bq - a));
                return s;
            };
            std::function<scalar(scalar, scalar, scalar, label)> rec =
                [&](scalar a, scalar bq, scalar parent, label depth) -> scalar
                {
                    const scalar m = 0.5*(a + bq);
                    const scalar c0 = gauss(a, m), c1 = gauss(m, bq);
                    if ((depth >= 4 && mag(c0 + c1 - parent) <= tol) || depth >= 40) return c0 + c1;
                    return rec(a, m, c0, depth + 1) + rec(m, bq, c1, depth + 1);
                };
            return rec(ya, yb, gauss(ya, yb), 0);
        }
    };

    // The second series: the frame tilted as at the sphere's equatorial cells
    // (e3 with a small component along the thin direction x) and a
    // non-diagonal tensor; the split quadrature at a tightened tolerance and
    // depth (1e-15 of the prism, depth 20) is the comparison next to the
    // split quadrature at the run's settings, and for the axis-aligned frame
    // (a = b = 0) the exact per-face references above
    Info<< "kromerBothe self-test, tilted frame: e3 = (sin a, -cos a sin b, cos a cos b), b12 = 0.5" << nl;
    for (const scalar aTilt : {0.0, 0.1, 0.3})
    {
        for (const scalar bTilt : {0.0, 0.33})
        {
            davof::cellSurface T;
            T.paraboloid = true;
            T.o = point::zero;
            T.e3 = vector(std::sin(aTilt), -std::cos(aTilt)*std::sin(bTilt), std::cos(aTilt)*std::cos(bTilt));
            davof::frameOf(T.e3, T.e1, T.e2);
            T.b11 = b; T.b22 = b; T.b12 = 0.5; T.c = 0; T.a1 = 0; T.a2 = 0;
            auto phiT = [&](const point& x) { return T.value(x); };
            for (const scalar t1 : {0.0, 0.02})
            {
                for (const scalar zShift : {0.5, 0.2, 0.8})
                {
                    const scalar t2 = 0.03;
                    // a point of the surface and the box about it
                    const scalar hs = T.c + 0.5*(T.b11*t1*t1 + 2*T.b12*t1*t2 + T.b22*t2*t2);
                    const point ps = T.o + t1*T.e1 + t2*T.e2 + hs*T.e3;
                    const point xin = ps + (0.5 - zShift)*h*vector(0, 0, 1);
                    const scalar x0 = xin.x() - 0.5*dx, x1 = xin.x() + 0.5*dx;
                    const scalar y0 = xin.y() - 0.5*h, y1 = xin.y() + 0.5*h;
                    const scalar z0 = xin.z() - 0.5*h, z1 = xin.z() + 0.5*h;
                    const scalar pv = dx*h*h;
                    DynamicList<davof::kbFace> faces(12);
                    auto quad = [&](point a, point bq, point cq, point d) {
                        List<point> poly(4); poly[0] = a; poly[1] = bq; poly[2] = cq; poly[3] = d;
                        davof::kbFanFaces(poly, xin, faces);
                    };
                    quad(point(x0,y0,z0), point(x1,y0,z0), point(x1,y1,z0), point(x0,y1,z0));
                    quad(point(x0,y0,z1), point(x1,y0,z1), point(x1,y1,z1), point(x0,y1,z1));
                    quad(point(x0,y0,z0), point(x1,y0,z0), point(x1,y0,z1), point(x0,y0,z1));
                    quad(point(x0,y1,z0), point(x1,y1,z0), point(x1,y1,z1), point(x0,y1,z1));
                    quad(point(x0,y0,z0), point(x0,y1,z0), point(x0,y1,z1), point(x0,y0,z1));
                    quad(point(x1,y0,z0), point(x1,y1,z0), point(x1,y1,z1), point(x1,y0,z1));
                    scalar vKb = 0;
                    label nLeavesKb = 0;
                    const bool ok = davof::kbPolyhedronLiquidVolume(faces, T, tolRel*pv, vKb, nLeavesKb);
                    const vector d(dx, 0, 0);
                    const point fc(x1, 0.5*(y0 + y1), 0.5*(z0 + z1));
                    const point q0(x1,y0,z0), q1(x1,y1,z0), q2(x1,y1,z1), q3(x1,y0,z1);
                    auto splitAt = [&](const scalar tolA, const label depth) -> scalar {
                        label nls = 0;
                        scalar v = 0;
                        v += davof::prismLiquidVolume(fc, q0, q1, d, phiT, tolA, depth, nls);
                        v += davof::prismLiquidVolume(fc, q1, q2, d, phiT, tolA, depth, nls);
                        v += davof::prismLiquidVolume(fc, q2, q3, d, phiT, tolA, depth, nls);
                        v += davof::prismLiquidVolume(fc, q3, q0, d, phiT, tolA, depth, nls);
                        return v;
                    };
                    const scalar vSplit = splitAt(tolRel*pv, maxDepth);
                    const scalar vTight = splitAt(1e-15*pv, 20);
                    Info<< "  tilt a " << aTilt << " b " << bTilt << " t1 " << t1 << " zShift " << zShift
                        << ": kb " << vKb/pv << (ok ? "" : " (NOT RESOLVED)") << " split " << vSplit/pv
                        << " splitTight " << vTight/pv
                        << "  (kb - tight)/pv " << (vKb - vTight)/pv
                        << "  (split - tight)/pv " << (vSplit - vTight)/pv << nl;
                    if (mag(vKb - vTight) > 1e-9*pv)
                    {
                        // the Kromer-Bothe terms per box face (the fan
                        // triangles summed): wetted area, face term, line
                        scalar wetF[6] = {0, 0, 0, 0, 0, 0}, lineF[6] = {0, 0, 0, 0, 0, 0}, termF[6] = {0, 0, 0, 0, 0, 0};
                        forAll(faces, k)
                        {
                            davof::kbSums sums;
                            davof::kbIntegrateFace(faces[k], T, 1e-16*pv, 1e-16*pv, true, sums);
                            const label bf = k/4;   // four fan triangles per box face, in the order quad() was called
                            wetF[bf] += sums.lastArea; lineF[bf] += sums.lineSum; termF[bf] += sums.faceSum;
                        }
                        const char* names[6] = {"z0 (n -z)", "z1 (n +z)", "y0 (n -y)", "y1 (n +y)", "x0 (n -x)", "x1 (n +x)"};
                        scalar lineTot = 0, termTot = 0;
                        for (label bf = 0; bf < 6; ++bf) { lineTot += lineF[bf]; termTot += termF[bf]; }
                        if (aTilt == 0 && bTilt == 0)
                        {
                            // exact references: e3 = z, e1 = x, e2 = y
                            AxisRef R{T.c, T.b11, T.b12, T.b22};
                            const scalar tolA = 1e-16*h*h;
                            // z-faces: the area where hz > z_c, as the clamped height with a unit band
                            auto areaZ = [&](scalar zc) {
                                return R.outer([&](scalar y) {
                                    // length of {x : hz(x, y) > zc} = clamp with H -> the indicator: use a thin band of height 1 scaled
                                    const scalar eps = 1e-9*h;
                                    return R.clampX(y, x0, x1, zc, eps)/eps;
                                }, y0, y1, tolA);
                            };
                            const scalar refZ0 = areaZ(z0), refZ1 = areaZ(z1);
                            const scalar refY0 = R.clampX(y0, x0, x1, z0, h), refY1 = R.clampX(y1, x0, x1, z0, h);
                            // x-faces: clamp along y at fixed x (hz quadratic in y): reuse clampX with the roles swapped
                            AxisRef Rs{T.c, T.b22, T.b12, T.b11};
                            const scalar refX0 = Rs.clampX(x0, y0, y1, z0, h), refX1 = Rs.clampX(x1, y0, y1, z0, h);
                            const scalar refLine = R.outer([&](scalar y) { return R.bandX(y, x0, x1, z0, z1); }, y0, y1, 1e-16*pv);
                            const scalar refWet[6] = {refZ0, refZ1, refY0, refY1, refX0, refX1};
                            // the face terms from the exact areas: (x_k - o).n_k A_k
                            const scalar dist[6] = {-z0, z1, -y0, y1, -x0, x1};
                            scalar refTerm = 0;
                            for (label bf = 0; bf < 6; ++bf)
                            {
                                refTerm += dist[bf]*refWet[bf];
                                Info<< "    " << names[bf] << ": wet kb " << wetF[bf] << " exact " << refWet[bf]
                                    << " diff/area " << (wetF[bf] - refWet[bf])/(bf < 2 ? dx*h : (bf < 4 ? dx*h : h*h))
                                    << "  line kb " << lineF[bf] << nl;
                            }
                            Info<< "    sums: faceTerm kb " << termTot << " exact " << refTerm
                                << "  line kb " << lineTot << " exact (int_S g) " << refLine
                                << "  V kb " << (termTot + lineTot)/3 << " exact " << (refTerm + refLine)/3
                                << " tight " << vTight << nl;
                            // the split quadrature's integrand against the exact liquid
                            // length of the segment [p - d, p] on a grid of the face
                            scalar maxDiff = 0; point worst = point::zero;
                            for (label i = 0; i < 60; ++i)
                            {
                                for (label j = 0; j < 60; ++j)
                                {
                                    const scalar y = y0 + (i + 0.5)/60*(y1 - y0), z = z0 + (j + 0.5)/60*(z1 - z0);
                                    const point p(x1, y, z);
                                    const scalar fSplit = davof::segmentLiquidFraction(p - d, p, phiT)*dx;
                                    const scalar fExact = R.clampX(y, x0, x1, z, 1e30)/1e30*0 + R.clampX(y, x0, x1, z, h*1e-9)/(h*1e-9);
                                    if (mag(fSplit - fExact) > maxDiff) { maxDiff = mag(fSplit - fExact); worst = p; }
                                }
                            }
                            Info<< "    split integrand vs exact segment length on a 60x60 grid: max diff/dx "
                                << maxDiff/dx << " at " << worst << nl;
                            // the sweep direction of the split quadrature's triangles reversed
                            label nlr = 0;
                            scalar vRev = 0;
                            vRev += davof::prismLiquidVolume(q0, q1, fc, d, phiT, 1e-15*pv, 20, nlr);
                            vRev += davof::prismLiquidVolume(q1, q2, fc, d, phiT, 1e-15*pv, 20, nlr);
                            vRev += davof::prismLiquidVolume(q2, q3, fc, d, phiT, 1e-15*pv, 20, nlr);
                            vRev += davof::prismLiquidVolume(q3, q0, fc, d, phiT, 1e-15*pv, 20, nlr);
                            // the split quadrature on the two triangles of the face's diagonal
                            label nld = 0;
                            scalar vDiag = 0;
                            vDiag += davof::prismLiquidVolume(q0, q1, q2, d, phiT, 1e-15*pv, 20, nld);
                            vDiag += davof::prismLiquidVolume(q0, q2, q3, d, phiT, 1e-15*pv, 20, nld);
                            Info<< "    split tight: fan " << vTight/pv << " (sweep reversed " << vRev/pv
                                << ", leaves " << nlr << "), diagonal triangles " << vDiag/pv
                                << " (leaves " << nld << "), exact " << (refTerm + refLine)/3/pv << nl;
                        }
                        else
                        {
                            for (label bf = 0; bf < 6; ++bf)
                                Info<< "    " << names[bf] << ": wet kb " << wetF[bf] << " line kb " << lineF[bf] << nl;
                        }
                    }
                }
            }
        }
    }
    Info<< endl;
}


// The faces of a polyhedron and their wetted areas under S, for the debug print
void printFaces
(
    const Foam::DynamicList<Foam::davof::kbFace>& faces,
    const Foam::davof::cellSurface& S,
    const Foam::scalar tolAbs
)
{
    using namespace Foam;
    forAll(faces, k)
    {
        const davof::kbFace& F = faces[k];
        davof::kbSums sums;
        const bool ok = davof::kbIntegrateFace(F, S, tolAbs, tolAbs, true, sums);
        Pout<< "    face " << k << " n " << F.n << " x0 " << F.x0 << " Ls " << F.Ls
            << " area " << F.area() << " wet " << sums.lastArea
            << " faceTerm " << sums.faceSum << " line " << sums.lineSum
            << " arcs " << sums.nArcs << (ok ? "" : " NOT RESOLVED") << nl;
        forAll(F.p, i)
        {
            Pout<< "      p " << F.p[i] << " phi " << S.value(F.p[i]) << nl;
        }
    }
}

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::davofVolumeRules::kromerBothe::cellLiquidVolume
(
    const label celli,
    const davof::cellSurface& S,
    const scalar tolRel
) const
{
    if (S.sphere)
    {
        ++nFallback_;
        return splitCellLiquidVolume(celli, S, tolRel);
    }
    DynamicList<davof::kbFace> faces(24);
    davof::kbCellFaces(mesh_, celli, faces);
    scalar v = 0;
    const scalar Vc = mesh_.V()[celli];
    if
    (
        davof::kbPolyhedronLiquidVolume(faces, S, tolRel*Vc, v, nLeaves_)
     && v >= -1e-9*Vc && v <= (1 + 1e-9)*Vc
    )
    {
        if (debug_)
        {
            const scalar vs = splitCellLiquidVolume(celli, S, tolRel);
            const scalar rel = mag(v - vs)/Vc;
            dbgMaxCell_ = max(dbgMaxCell_, rel);
            dbgSumCell_ += rel;
            ++dbgNCell_;
        }
        return min(max(v, scalar(0)), Vc);
    }
    ++nFallback_;
    return splitCellLiquidVolume(celli, S, tolRel);
}


Foam::scalar Foam::davofVolumeRules::kromerBothe::prismLiquidVolume
(
    const label facei,
    const vector& d,
    const davof::cellSurface& S,
    const scalar tolRel
) const
{
    const vector& Sf = mesh_.faceAreas()[facei];
    const scalar prismVol = mag(d & Sf);
    // a face parallel to the displacement sweeps no volume (relative test:
    // the flat prism of a round-off product is not a closed polyhedron and
    // the face sum does not cancel)
    if (prismVol <= 1e-10*mag(d)*mag(Sf)) return 0;
    if (S.sphere)
    {
        ++nFallback_;
        return splitPrismLiquidVolume(facei, d, S, tolRel);
    }
    DynamicList<davof::kbFace> faces(12);
    davof::kbPrismFaces(mesh_, facei, d, faces);
    scalar v = 0;
    if
    (
        davof::kbPolyhedronLiquidVolume(faces, S, tolRel*prismVol, v, nLeaves_)
     && v >= -1e-9*prismVol && v <= (1 + 1e-9)*prismVol
    )
    {
        v = min(max(v, scalar(0)), prismVol);
        if (debug_)
        {
            const scalar vs = splitPrismLiquidVolume(facei, d, S, tolRel);
            const scalar rel = mag(v - vs)/prismVol;
            dbgMaxPrism_ = max(dbgMaxPrism_, rel);
            dbgSumPrism_ += rel;
            ++dbgNPrism_;
            if (nDebugPrinted_ < 12 && mag(v - vs) > 1e-8*prismVol)
            {
                ++nDebugPrinted_;
                Pout<< "kromerBothe prism mismatch: face " << facei << " d " << d
                    << " prismVol " << prismVol << " kb " << v << " split " << vs
                    << " diff/prismVol " << (v - vs)/prismVol << nl;
                // a third value: the brute-force midpoint grid over the face
                // (a parallelogram spanned from its first vertex) with the
                // exact liquid fraction of every segment [p - d, p], at
                // three grid sizes and Richardson-extrapolated (the
                // integrand has kinks on the traces, so the grid converges
                // at second order at best)
                const face& f = mesh_.faces()[facei];
                if (f.size() == 4)
                {
                    const pointField& P = mesh_.points();
                    const point& p0 = P[f[0]];
                    const vector ea = P[f[1]] - p0, eb = P[f[3]] - p0;
                    auto phi = [&](const point& x) { return S.value(x); };
                    scalar vb[3];
                    const label ns[3] = {200, 400, 800};
                    for (label k = 0; k < 3; ++k)
                    {
                        const label n = ns[k];
                        scalar sum = 0;
                        for (label i = 0; i < n; ++i)
                        {
                            for (label j = 0; j < n; ++j)
                            {
                                const point p = p0 + ((i + 0.5)/n)*ea + ((j + 0.5)/n)*eb;
                                sum += davof::segmentLiquidFraction(p - d, p, phi);
                            }
                        }
                        vb[k] = prismVol*sum/(scalar(n)*n);
                    }
                    const scalar rich = vb[2] + (vb[2] - vb[1])/3.0;
                    Pout<< "  brute force (midpoint grid, exact segments) n=200 " << vb[0]
                        << " n=400 " << vb[1] << " n=800 " << vb[2]
                        << " Richardson " << rich << " (kb - bf)/prismVol " << (v - rich)/prismVol
                        << " (split - bf)/prismVol " << (vs - rich)/prismVol << nl;
                }
                printFaces(faces, S, tolRel*prismVol);
            }
        }
        return v;
    }
    ++nFallback_;
    return splitPrismLiquidVolume(facei, d, S, tolRel);
}


void Foam::davofVolumeRules::kromerBothe::report() const
{
    if (!debug_) return;
    scalar mp = dbgMaxPrism_, sp = dbgSumPrism_, mc = dbgMaxCell_, sc = dbgSumCell_;
    label np = dbgNPrism_, nc = dbgNCell_;
    if (Pstream::parRun())
    {
        reduce(mp, maxOp<scalar>()); reduce(sp, sumOp<scalar>());
        reduce(mc, maxOp<scalar>()); reduce(sc, sumOp<scalar>());
        reduce(np, sumOp<label>()); reduce(nc, sumOp<label>());
    }
    Info<< "kromerBothe debug: |kb - split| over " << np << " prisms max "
        << mp << " mean " << sp/max(np, label(1)) << " of the prism volume; over "
        << nc << " cell volumes max " << mc << " mean " << sc/max(nc, label(1))
        << " of the cell volume" << endl;
    dbgMaxPrism_ = dbgSumPrism_ = dbgMaxCell_ = dbgSumCell_ = 0;
    dbgNPrism_ = dbgNCell_ = 0;
    nDebugPrinted_ = 0;
}


Foam::scalar Foam::davofVolumeRules::kromerBothe::faceLiquidFraction
(
    const label facei,
    const vector& d,
    const davof::cellSurface& S,
    label& nFallback
) const
{
    scalar af = 0;
    if (!S.sphere && davof::kbFaceLiquidFraction(mesh_, facei, d, S, tolRel_, af, nLeaves_))
    {
        if (debug_ && nDebugPrinted_ < 12)
        {
            label nfb = 0;
            const scalar aq = davofVolumeRule::faceLiquidFraction(facei, d, S, nfb);
            if (mag(af - aq) > 1e-8)
            {
                ++nDebugPrinted_;
                Pout<< "kromerBothe face mismatch: face " << facei << " d " << d
                    << " kb " << af << " quadratic cut " << aq << nl;
            }
        }
        return af;
    }
    ++nFallback_;
    return davofVolumeRule::faceLiquidFraction(facei, d, S, nFallback);
}


// ************************************************************************* //
