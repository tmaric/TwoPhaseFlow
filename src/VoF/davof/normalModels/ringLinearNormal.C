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

#include "ringLinearNormal.H"
#include "davof.H"
#include "zoneDistribute.H"
#include "addToRunTimeSelectionTable.H"
#include "profiling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofNormalModels
{
    defineTypeNameAndDebug(ringLinear, 0);
    addToRunTimeSelectionTable(davofNormalModel, ringLinear, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofNormalModels::ringLinear::ringLinear
(
    const dictionary& dict,
    const reconstruction::davof& scheme
)
:
    davofNormalModel(typeName, dict, scheme, "all"),
    passes_(dict.getOrDefault<label>("ringLinearPasses", 1)),
    correction_(dict.getOrDefault<Switch>("ringLinearCorrection", Switch(true))),
    diag_(dict.getOrDefault<Switch>("ringLinearDiag", Switch(false)))
{
    if (passes_ != 1 && passes_ != 2)
    {
        FatalIOErrorInFunction(dict)
            << "normalModel ringLinear: ringLinearPasses " << passes_
            << " is not 1 or 2." << exit(FatalIOError);
    }
    if (scheme.interfaceCellRule() != "alphaAndRingSum")
    {
        WarningInFunction
            << "normalModel ringLinear with interfaceCellRule "
            << scheme.interfaceCellRule()
            << ": the ring sum is computed by the rule alphaAndRingSum only"
            << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::davofNormalModels::ringLinear::positionPlanes(const vectorField& n)
{
    reconstruction::davof& scheme = const_cast<reconstruction::davof&>(scheme_);
    const DynamicField<label>& labels = scheme_.interfaceLabels();
    forAll(labels, i)
    {
        scheme.positionPlane(labels[i], n[i]);
    }
    const_cast<volVectorField&>(scheme_.normal()).correctBoundaryConditions();
    const_cast<volVectorField&>(scheme_.centre()).correctBoundaryConditions();
}


void Foam::davofNormalModels::ringLinear::fit(vectorField& n, const label pass)
{
    zoneDistribute& exchange = zoneDistribute::New(mesh_);
    const boolList& zone = scheme_.interfaceCell();
    exchange.setUpCommforZone(zone, true);

    const volVectorField& centre = scheme_.centre();
    const volVectorField& normal = scheme_.normal();
    const volVectorField& mField = scheme_.areaNormal();

    const Map<vector> mapM(exchange.getDatafromOtherProc(zone, mField));
    const Map<vector> mapX(exchange.getDatafromOtherProc(zone, centre));
    const Map<vector> mapN(exchange.getDatafromOtherProc(zone, normal));

    const labelListList& stencil = exchange.getStencil();
    const DynamicField<label>& labels = scheme_.interfaceLabels();
    const vectorField& C = mesh_.cellCentres();

    const bool is2D = (mesh_.nGeometricD() == 2);
    vector e(Zero);
    if (is2D)
    {
        const Vector<label>& gd = mesh_.geometricD();
        for (direction d = 0; d < 3; ++d)
        {
            if (gd[d] == -1) e[d] = 1;
        }
    }
    const label minHold = is2D ? 2 : 3;

    DynamicList<vector> ms(32);
    DynamicList<vector> xs(32);
    DynamicList<scalar> ws(32);
    DynamicList<scalar> q1s(32);
    DynamicList<scalar> q2s(32);

    label nYZero = 0;       // y = 0 by the solvability condition
    label nNegWeight = 0;   // a negative weight for some j
    scalar maxWy = 0;       // max over cells of max_j |W (q_j - qbar).y|

    vectorField nNew(labels.size(), Zero);

    forAll(labels, i)
    {
        const label celli = labels[i];
        const point xk = (mag(normal[celli]) > 0) ? centre[celli] : C[celli];

        ms.clear(); xs.clear(); ws.clear();
        scalar W = 0;
        vector msum(Zero);
        for (const label gblIdx : stencil[celli])
        {
            const vector nk = exchange.getValue(normal, mapN, gblIdx);
            if (!(mag(nk) > 0)) continue;       // holds no plane
            const vector mj = exchange.getValue(mField, mapM, gblIdx);
            const vector xj = exchange.getValue(centre, mapX, gblIdx);
            ms.append(mj);
            xs.append(xj);
            ws.append(mag(mj));
            W += mag(mj);
            msum += mj;
        }

        if (!(W > 0))
        {
            continue;   // no weight: no normal (nNew stays zero)
        }

        const vector abar = msum/W;
        const scalar abarMag = mag(abar);
        if (!(abarMag > 0))
        {
            continue;
        }
        vector nbar = abar/abarMag;

        // The frame and the in-surface coordinates of the ring's centroids
        vector t1, t2;
        if (is2D)
        {
            nbar -= (nbar & e)*e;
            const scalar nbMag = mag(nbar);
            if (nbMag > 0)
            {
                nbar /= nbMag;
            }
            t1 = nbar ^ e;
            const scalar t1Mag = mag(t1);
            if (t1Mag > 0)
            {
                t1 /= t1Mag;
            }
            t2 = Zero;
        }
        else
        {
            davofNormalModel::frameOf(nbar, t1, t2);
        }

        q1s.clear(); q2s.clear();
        scalar q1bar = 0, q2bar = 0;
        forAll(xs, j)
        {
            const vector d = xs[j] - xk;
            const scalar q1 = (d & t1);
            const scalar q2 = is2D ? scalar(0) : (d & t2);
            q1s.append(q1);
            q2s.append(q2);
            q1bar += ws[j]*q1;
            q2bar += ws[j]*q2;
        }
        q1bar /= W;
        q2bar /= W;

        scalar c11 = 0, c12 = 0, c22 = 0;
        forAll(q1s, j)
        {
            const scalar d1 = q1s[j] - q1bar;
            const scalar d2 = q2s[j] - q2bar;
            c11 += ws[j]*d1*d1;
            c12 += ws[j]*d1*d2;
            c22 += ws[j]*d2*d2;
        }

        scalar y1 = 0, y2 = 0;
        bool solvable = (ms.size() >= minHold);
        if (correction_)
        {
            if (solvable)
            {
                if (is2D)
                {
                    if (c11 != 0) { y1 = q1bar/c11; } else { solvable = false; }
                }
                else
                {
                    const scalar det = c11*c22 - c12*c12;
                    if (det != 0)
                    {
                        y1 = (c22*q1bar - c12*q2bar)/det;
                        y2 = (c11*q2bar - c12*q1bar)/det;
                    }
                    else
                    {
                        solvable = false;
                    }
                }
            }
            if (!solvable) ++nYZero;
        }

        vector a(Zero);
        bool neg = false;
        forAll(ms, j)
        {
            const scalar wy =
                (y1*(q1s[j] - q1bar) + y2*(q2s[j] - q2bar));
            const scalar wgt = 1.0/W - wy;
            a += ms[j]*wgt;
            if (wgt < 0) neg = true;
            maxWy = max(maxWy, mag(W*wy));
        }
        if (neg) ++nNegWeight;

        const scalar aMag = mag(a);
        nNew[i] = (aMag > 0) ? a/aMag : vector(Zero);
    }

    n = nNew;

    if (diag_)
    {
        reduce(nYZero, sumOp<label>());
        reduce(nNegWeight, sumOp<label>());
        reduce(maxWy, maxOp<scalar>());
        Info<< "davof: ringLinear pass " << pass << ": y = 0 by solvability in "
            << nYZero << " interface cells, negative weight in " << nNegWeight
            << ", max |W (q - qbar).y| " << maxWy << endl;
    }
}


void Foam::davofNormalModels::ringLinear::modelNormals
(
    const scalar,
    vectorField& n
)
{
    addProfilingInFunction(geometricVoF);
    const vectorField& d = scheme_.ringSum().primitiveField();
    const DynamicField<label>& labels = scheme_.interfaceLabels();

    // Stage 1: the ring sums, the planes positioned along them
    n.resize(labels.size());
    forAll(labels, i)
    {
        n[i] = d[labels[i]];
    }
    positionPlanes(n);

    // Stage 2 (and 3 between the passes)
    for (label pass = 1; pass <= passes_; ++pass)
    {
        fit(n, pass);
        if (pass < passes_)
        {
            positionPlanes(n);
        }
    }
}


// ************************************************************************* //
