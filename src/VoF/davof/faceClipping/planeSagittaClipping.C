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

#include "planeSagittaClipping.H"
#include "davof.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofFaceClippingModels
{
    defineTypeNameAndDebug(planeSagitta, 0);
    addToRunTimeSelectionTable(davofFaceClipping, planeSagitta, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofFaceClippingModels::planeSagitta::planeSagitta
(
    const dictionary& dict,
    const reconstruction::davof& scheme
)
:
    davofFaceClipping(typeName, dict, scheme),
    minCos_(dict.getOrDefault<scalar>("sagittaMinCos", 0.2)),
    minChord_(dict.getOrDefault<scalar>("sagittaMinChord", 0.05)),
    maxFraction_(dict.getOrDefault<scalar>("sagittaMaxFraction", 0.25)),
    form_(dict.getOrDefault<word>("sagittaForm", "full")),
    tensorSource_(dict.getOrDefault<word>("tensorSource", "ring"))
{
    if (form_ != "full" && form_ != "segment")
    {
        FatalIOErrorInFunction(dict)
            << "faceClipping planeSagitta: unknown sagittaForm '" << form_
            << "'. Valid: full, segment." << exit(FatalIOError);
    }
    if (tensorSource_ != "ring")
    {
        FatalIOErrorInFunction(dict)
            << "faceClipping planeSagitta: unknown tensorSource '"
            << tensorSource_ << "'. Valid: ring." << exit(FatalIOError);
    }
    Info<< "isoAdvection: faceClipping planeSagitta (sagittaForm " << form_
        << ", tensorSource " << tensorSource_ << ", sagittaMinCos " << minCos_
        << ", sagittaMinChord " << minChord_ << ", sagittaMaxFraction "
        << maxFraction_ << ")" << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::davofFaceClippingModels::planeSagitta::correct
(
    const label facei,
    const label planeCell,
    const vector& nLiq,
    const point& xEnd,
    const UList<point>& chord,
    const scalar areaPlane
)
{
    // No trace of the plane on the face at the end of the step: the
    // surface's trace can differ from the chord's absence by a corner
    // triangle of O(h^4) only
    if (chord.size() < 2)
    {
        return areaPlane;
    }

    // The cell's tensor (zero: the fit fell back, counted by the scheme)
    const symmTensor& B = scheme_.curvatureTensor().primitiveField()[planeCell];
    const scalar z0 = scheme_.curvatureOffset().primitiveField()[planeCell];
    if (B == symmTensor(Zero))
    {
        return areaPlane;
    }

    // The chord: its two points, or the two farthest apart of a
    // degenerate set (a vertex on the plane repeats a crossing)
    label i1 = 0, i2 = 1;
    if (chord.size() > 2)
    {
        scalar dMax = -1;
        forAll(chord, i)
        {
            for (label j = i + 1; j < chord.size(); ++j)
            {
                const scalar d = magSqr(chord[j] - chord[i]);
                if (d > dMax)
                {
                    dMax = d;
                    i1 = i;
                    i2 = j;
                }
            }
        }
    }
    const point& p1 = chord[i1];
    const point& p2 = chord[i2];

    const vector& Sf = mesh_.faceAreas()[facei];
    const scalar magSf = mag(Sf);
    const vector nf = Sf/magSf;

    const scalar L = mag(p2 - p1);
    if (L < minChord_*Foam::sqrt(magSf))
    {
        ++nSkippedChord_;
        return areaPlane;
    }
    const vector t = (p2 - p1)/L;

    // The angle between the surface normal and the face plane: the
    // in-plane offset of the trace is the normal offset over cos(phi)
    const vector nu = nLiq - (nLiq & nf)*nf;
    const scalar cosPhi = mag(nu);
    if (cosPhi < minCos_)
    {
        ++nSkippedCos_;
        return areaPlane;
    }

    const vector Bt = B & t;
    const scalar tBt = t & Bt;

    scalar dA = 0;
    if (form_ == "full")
    {
        // The integral of the surface's offset w(l) = z_0 + (1/2) x.B.x
        // from the plane along the chord, x = a + l t, a from the
        // centroid at the end of the step
        const vector a = p1 - xEnd;
        const scalar aBt = a & Bt;
        const scalar aBa = a & (B & a);
        const scalar w = z0*L + 0.5*(aBa*L + aBt*L*L + tBt*L*L*L/3.0);
        dA = -w/cosPhi;
    }
    else
    {
        // The parabolic segment between the chord and the arc through its
        // end points, (2/3) L s with the sagitta s = kappa L^2/8
        dA = tBt*L*L*L/(12.0*cosPhi);
    }

    const scalar dAMax = maxFraction_*magSf;
    if (mag(dA) > dAMax)
    {
        dA = sign(dA)*dAMax;
        ++nCapped_;
    }
    ++nCorrected_;

    scalar area = areaPlane + dA;
    if (area < 0 || area > magSf)
    {
        if (area < -1e-12*magSf || area > (1 + 1e-12)*magSf)
        {
            ++nClamped_;
        }
        area = min(max(area, scalar(0)), magSf);
    }
    return area;
}


// ************************************************************************* //
