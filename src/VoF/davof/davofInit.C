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

\*---------------------------------------------------------------------------*/

#include "davofInit.H"
#include "davofSimplexGeometry.H"
#include "davofQuadraticFaceGeometry.H"
#include "davofSphereGeometry.H"
#include "boundBox.H"
#include "profiling.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::davofInit::psi(const point& x) const
{
    if (signedDistance_)
    {
        return -signedDistance(func_, x);
    }
    return -func_.value(x);
}


Foam::scalar Foam::davofInit::cellFromValues
(
    const label celli,
    const scalarField& psiC,
    const scalarField& psiF,
    const scalarField& psiP
)
{
    const pointField& points = mesh_.points();
    const faceList& faces = mesh_.faces();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& C = mesh_.cellCentres();
    const labelList& cFaces = mesh_.cells()[celli];

    const point& xc = C[celli];
    const scalar phic = psiC[celli];

    // One-signed cells carry no interface: alpha is 0 or 1, no tets.
    scalar phiMin = phic, phiMax = phic;
    forAll(cFaces, cf)
    {
        const label facei = cFaces[cf];
        phiMin = min(phiMin, psiF[facei]);
        phiMax = max(phiMax, psiF[facei]);
        const face& f = faces[facei];
        forAll(f, ip)
        {
            phiMin = min(phiMin, psiP[f[ip]]);
            phiMax = max(phiMax, psiP[f[ip]]);
        }
    }
    mTet_[celli] = Zero;
    if (phiMax <= 0) return 1;
    if (phiMin >  0) return 0;

    scalar liquidVol = 0, totalVol = 0;
    vector mt(Zero);

    // The Detrixhe-Aslam decomposition: one tet per face edge (cell centre,
    // face centre, edge end points); the loop order of leia's davofState.
    forAll(cFaces, cf)
    {
        const label facei = cFaces[cf];
        const face& f = faces[facei];
        const point& xf = Cf[facei];
        const scalar phif = psiF[facei];

        forAll(f, ip)
        {
            const label ip1 = f.nextLabel(ip);
            const point& p0 = points[f[ip]];
            const point& p1 = points[ip1];
            const scalar phi0 = psiP[f[ip]];
            const scalar phi1 = psiP[ip1];

            const scalar vol = mag(((xf - xc) ^ (p0 - xc)) & (p1 - xc))/6.0;
            if (vol <= VSMALL) continue;

            liquidVol += davof::tetNegativeFraction(phic, phif, phi0, phi1)*vol;
            totalVol  += vol;

            const point x[4] = {xc, xf, p0, p1};
            const scalar ph[4] = {phic, phif, phi0, phi1};
            vector av;
            point cen;
            scalar a;
            if (davof::tetZeroSet(x, ph, av, cen, a) > 0)
            {
                mt += av;
            }
        }
    }

    const scalar raw =
        (totalVol > VSMALL) ? liquidVol/totalVol : (phic < 0 ? 1 : 0);
    const scalar clamped = min(max(raw, scalar(0)), scalar(1));
    maxAlphaClip_ = max(maxAlphaClip_, mag(raw - clamped));

    // leia's tet patch vector points from psi <= 0 (liquid) to psi > 0;
    // TwoPhaseFlow's normal points into the liquid.
    mTet_[celli] = -mt;
    return clamped;
}


Foam::scalar Foam::davofInit::faceFraction
(
    const label facei,
    const scalarField& psiF,
    const scalarField& psiP
) const
{
    const pointField& points = mesh_.points();
    const face& f = mesh_.faces()[facei];
    const point& xf = mesh_.faceCentres()[facei];
    const scalar phif = psiF[facei];

    scalar liqA = 0, totA = 0;
    forAll(f, ip)
    {
        const label ip1 = f.nextLabel(ip);
        const point& p0 = points[f[ip]];
        const point& p1 = points[ip1];
        const scalar triA = 0.5*mag((p0 - xf) ^ (p1 - xf));
        if (triA <= VSMALL) continue;
        liqA +=
            davof::triNegativeFraction(phif, psiP[f[ip]], psiP[ip1])*triA;
        totA += triA;
    }
    return (totA > VSMALL) ? liqA/totA : (phif < 0 ? scalar(1) : scalar(0));
}


void Foam::davofInit::detrixheAslamState
(
    volScalarField& alpha,
    scalarField& alphafAll,
    scalarField& psiF,
    scalarField& psiP
)
{
    addProfilingInFunction(geometricVoF);
    const vectorField& C = mesh_.cellCentres();
    const vectorField& Cf = mesh_.faceCentres();
    const pointField& P = mesh_.points();

    scalarField psiC(C.size());
    forAll(C, i) psiC[i] = psi(C[i]);
    psiF.setSize(Cf.size());
    forAll(Cf, i) psiF[i] = psi(Cf[i]);
    psiP.setSize(P.size());
    forAll(P, i) psiP[i] = psi(P[i]);

    maxAlphaClip_ = 0;
    mTet_.setSize(mesh_.nCells());
    forAll(C, celli)
    {
        alpha[celli] = cellFromValues(celli, psiC, psiF, psiP);
    }

    alphafAll.setSize(mesh_.nFaces());
    forAll(alphafAll, facei)
    {
        alphafAll[facei] = faceFraction(facei, psiF, psiP);
    }
    alpha.correctBoundaryConditions();
}


void Foam::davofInit::fillFaceField
(
    const scalarField& alphafAll,
    surfaceScalarField& alphaf
) const
{
    scalarField& in = alphaf.primitiveFieldRef();
    forAll(in, facei)
    {
        in[facei] = alphafAll[facei];
    }
    surfaceScalarField::Boundary& bf = alphaf.boundaryFieldRef();
    forAll(bf, patchi)
    {
        const label start = mesh_.boundaryMesh()[patchi].start();
        fvsPatchScalarField& pf = bf[patchi];
        forAll(pf, i)
        {
            pf[i] = alphafAll[start + i];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofInit::davofInit
(
    const fvMesh& mesh,
    const implicitFunction& func,
    const dictionary& dict
)
:
    mesh_(mesh),
    func_(func),
    dict_(dict),
    signedDistance_(dict.getOrDefault<bool>("signedDistance", false)),
    maxAlphaClip_(0),
    nFaceFallback_(0),
    mTet_(mesh.nCells(), Zero)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::davofInit::detrixheAslam
(
    volScalarField& alpha,
    surfaceScalarField& alphaf
)
{
    scalarField alphafAll, psiF, psiP;
    detrixheAslamState(alpha, alphafAll, psiF, psiP);
    nFaceFallback_ = 0;
    fillFaceField(alphafAll, alphaf);
}


void Foam::davofInit::quadraticFaces
(
    volScalarField& alpha,
    surfaceScalarField& alphaf
)
{
    addProfilingInFunction(geometricVoF);
    scalarField alphafAll, psiF, psiP;
    detrixheAslamState(alpha, alphafAll, psiF, psiP);

    const pointField& P = mesh_.points();
    const faceList& faces = mesh_.faces();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();

    label nFallback = 0;
    forAll(alphafAll, facei)
    {
        const scalar A = mag(Sf[facei]);
        if (A <= VSMALL) continue;
        const face& f = faces[facei];
        const point& xf = Cf[facei];
        const scalar phif = psiF[facei];

        scalar liqA = 0, totA = 0;
        bool fb = false;
        forAll(f, ip)
        {
            const label ip1 = f.nextLabel(ip);
            const point& p0 = P[f[ip]];
            const point& p1 = P[ip1];
            const scalar triA = 0.5*mag((p0 - xf) ^ (p1 - xf));
            if (triA <= VSMALL) continue;
            // vertices xf, p0, p1; midpoints of (xf,p0), (p0,p1), (p1,xf)
            const scalar ph[6] =
            {
                phif,
                psiP[f[ip]],
                psiP[ip1],
                psi(0.5*(xf + p0)),
                psi(0.5*(p0 + p1)),
                psi(0.5*(p1 + xf))
            };
            bool tfb = false;
            const scalar fr =
                davof::triQuadraticNegativeFraction(xf, p0, p1, ph, tfb);
            if (tfb) { fb = true; break; }
            liqA += fr*triA;
            totA += triA;
        }
        if (fb) { ++nFallback; continue; }
        if (totA > VSMALL)
        {
            alphafAll[facei] = min(max(liqA/totA, scalar(0)), scalar(1));
        }
    }
    nFaceFallback_ = returnReduce(nFallback, sumOp<label>());

    fillFaceField(alphafAll, alphaf);
}


void Foam::davofInit::exactSphere
(
    volScalarField& alpha,
    surfaceScalarField& alphaf
)
{
    addProfilingInFunction(geometricVoF);
    // The Detrixhe-Aslam state first: the tet reference mTet_ (so that the
    // consistency check measures exact - chord, O(h) by construction) and
    // the fallbacks.
    scalarField alphafAll, psiF, psiP;
    detrixheAslamState(alpha, alphafAll, psiF, psiP);
    nFaceFallback_ = 0;

    const point c = dict_.get<point>("origin");
    const scalar R = dict_.get<scalar>("radius");

    const pointField& P = mesh_.points();
    const faceList& faces = mesh_.faces();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();
    const scalarField& V = mesh_.V();

    // Axis-aligned hexahedra: every face normal along one axis, and the
    // cell's bounding box has the cell's volume.
    forAll(Sf, facei)
    {
        const scalar A = mag(Sf[facei]);
        if (A <= VSMALL) continue;
        const vector n = Sf[facei]/A;
        if (cmptMax(cmptMag(n)) < 1 - 1e-10)
        {
            FatalErrorInFunction
                << "davofMethod exactSphere needs axis-aligned hexahedral "
                << "cells; face " << facei << " has the unit normal " << n
                << exit(FatalError);
        }
    }

    // Faces: the planar polygon inside the ball. A face with every vertex
    // inside the ball is inside (convex faces), fraction exactly 1; a polygon
    // disjoint from the circle sums its sectors to a round-off zero. Both
    // are snapped (1e-12), so the bulk carries exactly equal face fractions
    // and the area normal is zeroed there instead of counting round-off
    // wisps.
    const scalar snapTol = 1e-12;
    forAll(alphafAll, facei)
    {
        const scalar A = mag(Sf[facei]);
        if (A <= VSMALL) { alphafAll[facei] = 0; continue; }
        const face& f = faces[facei];
        bool allInside = true;
        forAll(f, i)
        {
            if (magSqr(P[f[i]] - c) > R*R) { allInside = false; break; }
        }
        if (allInside) { alphafAll[facei] = 1; continue; }
        scalar a = davof::faceBallArea(f, P, Cf[facei], Sf[facei]/A, c, R)/A;
        if (a < snapTol) a = 0;
        else if (a > 1 - snapTol) a = 1;
        alphafAll[facei] = min(max(a, scalar(0)), scalar(1));
    }

    // Cells: axis-aligned boxes.
    const labelListList& cellPoints = mesh_.cellPoints();
    forAll(alpha, celli)
    {
        const pointField cp(P, cellPoints[celli]);
        const boundBox bb(cp, false);
        const vector span = bb.span();
        const scalar Vbb = span.x()*span.y()*span.z();
        if (mag(Vbb - V[celli]) > 1e-8*V[celli])
        {
            FatalErrorInFunction
                << "davofMethod exactSphere needs axis-aligned hexahedral "
                << "cells; cell " << celli << " has volume " << V[celli]
                << " but bounding-box volume " << Vbb << exit(FatalError);
        }
        // Quick decisions: the box entirely inside or entirely outside.
        scalar farthest = 0, nearest2 = 0;
        for (label i = 0; i < 8; ++i)
        {
            const point corner
            (
                (i & 1) ? bb.max().x() : bb.min().x(),
                (i & 2) ? bb.max().y() : bb.min().y(),
                (i & 4) ? bb.max().z() : bb.min().z()
            );
            farthest = max(farthest, mag(corner - c));
        }
        for (direction k = 0; k < 3; ++k)
        {
            const scalar q = min(max(c[k], bb.min()[k]), bb.max()[k]) - c[k];
            nearest2 += q*q;
        }
        if (farthest <= R)
        {
            alpha[celli] = 1;
        }
        else if (nearest2 >= R*R)
        {
            alpha[celli] = 0;
        }
        else
        {
            const scalar a =
                davof::boxBallVolume(bb.min(), bb.max(), c, R)/V[celli];
            alpha[celli] = min(max(a, scalar(0)), scalar(1));
        }
    }
    maxAlphaClip_ = 0;
    alpha.correctBoundaryConditions();
    fillFaceField(alphafAll, alphaf);
}


Foam::scalar Foam::davofInit::wispThreshold
(
    const fvMesh& mesh,
    const label celli,
    const scalar wispTol
)
{
    return wispTol*Foam::pow(mesh.V()[celli], 2.0/3.0);
}


Foam::scalar Foam::davofInit::projectedDistance
(
    const implicitFunction& func,
    const point& x,
    point& foot
)
{
    // Newton along the gradient onto the zero set.
    auto toSurface = [&](point& y)
    {
        for (label it = 0; it < 50; ++it)
        {
            const scalar v = func.value(y);
            const vector g = func.grad(y);
            const scalar g2 = magSqr(g);
            if (g2 <= VSMALL) break;
            const vector dy = (v/g2)*g;
            y -= dy;
            if (mag(dy) <= 1e-15*(mag(y) + SMALL)) break;
        }
    };

    foot = x;
    toSurface(foot);

    // The foot point: move within the tangent plane towards x, re-project.
    for (label it = 0; it < 50; ++it)
    {
        const vector g = func.grad(foot);
        const scalar gm = mag(g);
        if (gm <= VSMALL) break;
        const vector n = g/gm;
        const vector d = x - foot;
        const vector t = d - (d & n)*n;
        if (mag(t) <= 1e-14*(mag(d) + mag(x) + SMALL)) break;
        foot += t;
        toSurface(foot);
    }
    return mag(x - foot);
}


Foam::scalar Foam::davofInit::signedDistance
(
    const implicitFunction& func,
    const point& x
)
{
    point foot;
    const scalar d = projectedDistance(func, x, foot);
    return (func.value(x) < 0) ? -d : d;
}


void Foam::davofInit::areaNormal
(
    const fvMesh& mesh,
    const surfaceScalarField& alphaf,
    vectorField& m,
    scalarField& alphafAll
)
{
    addProfilingInFunction(geometricVoF);
    const vectorField& Sf = mesh.faceAreas();
    const labelUList& own = mesh.faceOwner();
    const labelUList& nei = mesh.faceNeighbour();
    const label nInt = mesh.nInternalFaces();

    m.setSize(mesh.nCells());
    m = Zero;
    alphafAll.setSize(mesh.nFaces());
    alphafAll = 0;

    // A dry face (alpha_f = 0, the bulk of the gas) contributes 0 S_f =
    // (+-0, +-0, +-0), and x + (+-0) = x, x - (+-0) = x for every x
    // including +0 (the sum starts at +0, and +0 + -0 = +0, +0 - +0 = +0):
    // skipping it leaves every m bit-identical and saves the three
    // products and the six scattered additions on most faces.
    for (label facei = 0; facei < nInt; ++facei)
    {
        const scalar af = alphaf[facei];
        alphafAll[facei] = af;
        if (af == 0) continue;
        const vector contrib = af*Sf[facei];
        m[own[facei]] += contrib;
        m[nei[facei]] -= contrib;
    }
    const surfaceScalarField::Boundary& bf = alphaf.boundaryField();
    forAll(bf, patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        const fvsPatchScalarField& pf = bf[patchi];
        forAll(pf, i)
        {
            const label facei = pp.start() + i;
            alphafAll[facei] = pf[i];
            if (pf[i] == 0) continue;
            m[own[facei]] += pf[i]*Sf[facei];
        }
    }
}


Foam::scalar Foam::davofInit::maxConsistency
(
    const surfaceScalarField& alphaf,
    const scalar wispTol
) const
{
    vectorField m;
    scalarField alphafAll;
    areaNormal(mesh_, alphaf, m, alphafAll);
    const scalarField& V = mesh_.V();
    scalar worst = 0;
    forAll(m, celli)
    {
        if (mag(m[celli]) > wispThreshold(mesh_, celli, wispTol))
        {
            worst = max
            (
                worst,
                mag(m[celli] - mTet_[celli])/Foam::pow(V[celli], 2.0/3.0)
            );
        }
    }
    return returnReduce(worst, maxOp<scalar>());
}


// ************************************************************************* //
