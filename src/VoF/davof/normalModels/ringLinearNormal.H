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

Class
    Foam::davofNormalModels::ringLinear

Description
    The area-weighted linear fit of the ring's area vectors (Tomislav,
    2026-11-24). m_j is the area vector of cell j (the scheme's areaNormal,
    as the model ringSum uses it), w_j = |m_j|, x_j the centroid of the PLIC
    polygon of cell j.

    Stage 1: exactly ringSum: the ring sum in every interface cell, the
    plane of every interface cell positioned by alpha along it
    (davof::positionPlane). This gives x_j for every cell holding a plane.

    Stage 2, for every interface cell k, over the ring cells j (the
    point-neighbour ring of ringSum, k included) that hold a plane after
    stage 1: W = sum w_j, abar = sum m_j/W, nbar = abar/|abar|. In 3D
    (t1, t2) is an orthonormal pair orthogonal to nbar and q_j = ((x_j -
    x_k).t1, (x_j - x_k).t2); in 2D (e the empty direction) the e component
    is removed from nbar, t1 = nbar x e normalised and q_j = (x_j - x_k).t1.
    qbar = sum w_j q_j/W, C = sum w_j (q_j - qbar)(q_j - qbar)^T, C y =
    qbar (Cramer's rule in 3D, a division in 2D), a = sum m_j [1/W - (q_j -
    qbar).y], the normal of cell k is a/|a|: the weighted least-squares fit
    of a linear unit-normal field with the weights w_j, evaluated at x_k.
    With y = 0 it is the ring sum over the same cells. Solvability, not a
    tolerance: fewer than 3 ring cells holding a plane in 3D (fewer than 2
    in 2D) or an exactly zero determinant of C gives y = 0. x_k is the
    centroid of the plane of cell k when it holds one, otherwise the cell
    centre.

    Stage 3: the planes repositioned by alpha along the new normals; with
    ringLinearPasses 2 stages 2 and 3 are repeated once with the new x_j
    (m_j unchanged). The scheme positions the final plane along the
    returned normal.

    No fallback, no tolerance. x_j, the planes and m_j of the ring cells on
    other processors come through zoneDistribute's exchange, as in ringFit.

        normalModel            ringLinear;
        ringLinearPasses       1;      // 1 | 2
        ringLinearCorrection   true;   // false: y = 0 everywhere (ringSum)
        ringLinearDiag         false;  // print only, per pass

    The diagnostics (ringLinearDiag): the number of interface cells with
    y = 0 by the solvability condition, the number with a negative weight
    1/W - (q_j - qbar).y for some j, the maximum over the cells of max_j |W
    (q_j - qbar).y| (reduced over the processors).

SourceFiles
    ringLinearNormal.C

\*---------------------------------------------------------------------------*/

#ifndef davofNormalModels_ringLinear_H
#define davofNormalModels_ringLinear_H

#include "davofNormalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofNormalModels
{

/*---------------------------------------------------------------------------*\
                          Class ringLinear Declaration
\*---------------------------------------------------------------------------*/

class ringLinear
:
    public davofNormalModel
{
    // Private Data

        //- 1 | 2
        label passes_;

        //- false: y = 0 everywhere
        bool correction_;

        //- Print only
        bool diag_;


    // Private Member Functions

        //- Position the plane of every interface cell by alpha along n
        void positionPlanes(const vectorField& n);

        //- Stage 2 for every interface cell: the normals a/|a| into n
        void fit(vectorField& n, const label pass);


public:

    //- Runtime type information
    TypeName("ringLinear");


    //- Construct from the davof dictionary and the scheme
    ringLinear
    (
        const dictionary& dict,
        const reconstruction::davof& scheme
    );


    //- Destructor
    virtual ~ringLinear() = default;


    // Member Functions

        virtual void modelNormals(const scalar tState, vectorField& n);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace davofNormalModels
} // End inline namespace twoPhaseFlow
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
