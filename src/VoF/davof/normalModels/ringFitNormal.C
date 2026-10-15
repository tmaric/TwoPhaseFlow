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
#include "profiling.H"

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
    addProfilingInFunction(geometricVoF);
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

        // The frame of the raw normal at the cell's own plane centroid;
        // the fit is the base class's mongeFit (the same arithmetic in
        // the same order as the former in-place assembly and solve)
        const vector e3 = m[celli]/mMag;
        vector e1, e2;
        const point o = (mag(normal[celli]) > 0) ? centre[celli] : C[celli];
        const scalar h = Foam::cbrt(V[celli]);

        scalar b[6];
        if (!mongeFit(e3, o, h, xs, ns, ws, wSum, alignTol_, minCells_, b, e1, e2))
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
