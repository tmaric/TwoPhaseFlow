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

#include "ringFitNormal.H"
#include "davof.H"
#include "zoneDistribute.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofNormalModels
{
    defineTypeNameAndDebug(ringFit, 0);
    addToRunTimeSelectionTable(davofNormalModel, ringFit, dictionary);
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::davofNormalModels::ringFit::frameOf
(
    const vector& e3,
    vector& e1,
    vector& e2
)
{
    // The coordinate axis least aligned with e3, made orthogonal to it.
    direction k = 0;
    for (direction i = 1; i < 3; ++i)
    {
        if (mag(e3[i]) < mag(e3[k])) k = i;
    }
    vector a(Zero);
    a[k] = 1;
    e1 = a - (a & e3)*e3;
    e1 /= max(mag(e1), VSMALL);
    e2 = e3 ^ e1;
}


bool Foam::davofNormalModels::ringFit::choleskySolve
(
    scalar* A,
    scalar* b,
    const label n
)
{
    for (label j = 0; j < n; ++j)
    {
        scalar d = A[j*n + j];
        for (label k = 0; k < j; ++k) { d -= A[j*n + k]*A[j*n + k]; }
        if (d <= SMALL) { return false; }
        d = Foam::sqrt(d);
        A[j*n + j] = d;
        for (label i = j + 1; i < n; ++i)
        {
            scalar s = A[i*n + j];
            for (label k = 0; k < j; ++k) { s -= A[i*n + k]*A[j*n + k]; }
            A[i*n + j] = s/d;
        }
    }
    for (label i = 0; i < n; ++i)
    {
        scalar s = b[i];
        for (label k = 0; k < i; ++k) { s -= A[i*n + k]*b[k]; }
        b[i] = s/A[i*n + i];
    }
    for (label i = n - 1; i >= 0; --i)
    {
        scalar s = b[i];
        for (label k = i + 1; k < n; ++k) { s -= A[k*n + i]*b[k]; }
        b[i] = s/A[i*n + i];
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofNormalModels::ringFit::ringFit
(
    const dictionary& dict,
    const reconstruction::davof& scheme
)
:
    ringArea(dict, scheme),
    weight_(dict.getOrDefault<word>("weight", "area")),
    alignTol_(dict.getOrDefault<scalar>("alignTol", 0.2)),
    minCells_(dict.getOrDefault<label>("minCells", 6))
{
    if (weight_ != "area" && weight_ != "alpha")
    {
        FatalIOErrorInFunction(dict)
            << "normalModel ringFit: unknown weight '" << weight_
            << "'. Valid: area, alpha." << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::davofNormalModels::ringFit::modelNormals
(
    const scalar,
    vectorField& n
)
{
    fillInterfaceNormals();

    zoneDistribute& exchange = zoneDistribute::New(mesh_);
    const boolList& zone = scheme_.interfaceCell();
    exchange.setUpCommforZone(zone, true);

    const volVectorField& centre = scheme_.centre();
    const volVectorField& normal = scheme_.normal();
    const volScalarField& alpha = scheme_.alpha1();

    const Map<vector> mapM
    (
        exchange.getDatafromOtherProc(zone, mInterface_)
    );
    const Map<vector> mapX(exchange.getDatafromOtherProc(zone, centre));
    const Map<vector> mapN(exchange.getDatafromOtherProc(zone, normal));
    const Map<scalar> mapA(exchange.getDatafromOtherProc(zone, alpha));

    const labelListList& stencil = exchange.getStencil();
    const DynamicField<label>& labels = scheme_.interfaceLabels();
    const vectorField& m = scheme_.areaNormal().primitiveField();
    const vectorField& C = mesh_.cellCentres();
    const scalarField& V = mesh_.V();

    // The ring data of one cell
    DynamicList<vector> xs(32);
    DynamicList<vector> ns(32);
    DynamicList<scalar> ws(32);

    nFallback_ = 0;

    forAll(labels, i)
    {
        const label celli = labels[i];
        const scalar mMag = mag(m[celli]);

        // The ring: the area sum (the fallback) and the usable planes
        vector sum(Zero);
        xs.clear(); ns.clear(); ws.clear();
        scalar wSum = 0;
        for (const label gblIdx : stencil[celli])
        {
            const vector mk = exchange.getValue(mInterface_, mapM, gblIdx);
            const scalar mkMag = mag(mk);
            if (!(mkMag > 0)) continue;       // not an interface cell
            sum += mk;

            const vector nk = exchange.getValue(normal, mapN, gblIdx);
            const scalar nkMag = mag(nk);
            if (!(nkMag > 0)) continue;       // no positioned plane
            const vector xk = exchange.getValue(centre, mapX, gblIdx);

            scalar w = mkMag;
            if (weight_ == "alpha")
            {
                const scalar ak = exchange.getValue(alpha, mapA, gblIdx);
                w = max(ak*(1 - ak), scalar(0));
            }
            if (!(w > 0)) continue;

            xs.append(xk);
            ns.append(nk/nkMag);
            ws.append(w);
            wSum += w;
        }
        n[i] = sum;

        if (!(mMag > 0) || xs.size() < minCells_ || !(wSum > 0))
        {
            ++nFallback_;
            continue;
        }

        // The frame of the raw normal at the cell's own plane centroid
        const vector e3 = m[celli]/mMag;
        vector e1, e2;
        frameOf(e3, e1, e2);
        const point o = (mag(normal[celli]) > 0) ? centre[celli] : C[celli];
        const scalar h = Foam::cbrt(V[celli]);

        // Normal equations A q = b over q = (a, b, c, d, e, f)
        scalar A[36] = {0};
        scalar b[6] = {0};
        label nUsable = 0;

        forAll(xs, k)
        {
            const vector& nk = ns[k];
            const scalar a3 = nk & e3;
            if (a3 <= alignTol_) continue;
            ++nUsable;

            const vector r = xs[k] - o;
            const scalar U = (r & e1)/h;
            const scalar Vv = (r & e2)/h;
            const scalar W = (r & e3)/h;
            const scalar w = ws[k]/wSum;
            const scalar s1 = -(nk & e1)/a3;
            const scalar s2 = -(nk & e2)/a3;

            // position row: a + b U + c V + d U^2 + e U V + f V^2 = W
            {
                const scalar val[6] = {1, U, Vv, U*U, U*Vv, Vv*Vv};
                for (label p = 0; p < 6; ++p)
                {
                    b[p] += w*val[p]*W;
                    for (label q = 0; q < 6; ++q)
                    {
                        A[p*6 + q] += w*val[p]*val[q];
                    }
                }
            }
            // slope row 1: b + 2 d U + e V = s1
            {
                const label idx[3] = {1, 3, 4};
                const scalar val[3] = {1, 2*U, Vv};
                for (label p = 0; p < 3; ++p)
                {
                    b[idx[p]] += w*val[p]*s1;
                    for (label q = 0; q < 3; ++q)
                    {
                        A[idx[p]*6 + idx[q]] += w*val[p]*val[q];
                    }
                }
            }
            // slope row 2: c + e U + 2 f V = s2
            {
                const label idx[3] = {2, 4, 5};
                const scalar val[3] = {1, U, 2*Vv};
                for (label p = 0; p < 3; ++p)
                {
                    b[idx[p]] += w*val[p]*s2;
                    for (label q = 0; q < 3; ++q)
                    {
                        A[idx[p]*6 + idx[q]] += w*val[p]*val[q];
                    }
                }
            }
        }

        if (nUsable < minCells_ || !choleskySolve(A, b, 6))
        {
            ++nFallback_;
            continue;
        }

        // The fitted surface's normal at the origin: (-b, -c, 1)
        const vector nFit = -b[1]*e1 - b[2]*e2 + e3;
        if (mag(nFit) > 0 && mag(nFit) < ROOTVGREAT)
        {
            n[i] = nFit;
        }
        else
        {
            ++nFallback_;
        }
    }
}


// ************************************************************************* //
