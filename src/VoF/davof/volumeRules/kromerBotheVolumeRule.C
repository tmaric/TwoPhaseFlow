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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofVolumeRules::kromerBothe::kromerBothe
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    davofVolumeRule(mesh, dict),
    debug_(dict.getOrDefault<bool>("volumeRuleDebug", false)),
    nDebugPrinted_(0)
{}


namespace
{

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
    if (davof::kbPolyhedronLiquidVolume(faces, S, tolRel*mesh_.V()[celli], v, nLeaves_))
    {
        return v;
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
    const scalar prismVol = mag(d & mesh_.faceAreas()[facei]);
    if (prismVol <= VSMALL) return 0;
    if (S.sphere)
    {
        ++nFallback_;
        return splitPrismLiquidVolume(facei, d, S, tolRel);
    }
    DynamicList<davof::kbFace> faces(12);
    davof::kbPrismFaces(mesh_, facei, d, faces);
    scalar v = 0;
    if (davof::kbPolyhedronLiquidVolume(faces, S, tolRel*prismVol, v, nLeaves_))
    {
        if (debug_ && nDebugPrinted_ < 12)
        {
            const scalar vs = splitPrismLiquidVolume(facei, d, S, tolRel);
            if (mag(v - vs) > 1e-8*prismVol)
            {
                ++nDebugPrinted_;
                Pout<< "kromerBothe prism mismatch: face " << facei << " d " << d
                    << " prismVol " << prismVol << " kb " << v << " split " << vs
                    << " diff/prismVol " << (v - vs)/prismVol << nl;
                printFaces(faces, S, tolRel*prismVol);
            }
        }
        return v;
    }
    ++nFallback_;
    return splitPrismLiquidVolume(facei, d, S, tolRel);
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
