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
    Foam::davofNormalModels::ringFit

Description
    A weighted least-squares quadric through the ring's recovered planes:
    in the frame (e1, e2, e3 = m_c/|m_c|) with the origin at the cell's own
    plane centroid (the cell centre when the plane does not exist), scaled
    by h_c = V_c^(1/3), the Monge patch

        z = a + b x + c y + d x^2 + e x y + f y^2

    fitted to the plane centroids of the ring's interface cells (one
    position equation each) and their unit normals (two gradient equations
    each, dz/dx = -(n.e1)/(n.e3), dz/dy = -(n.e2)/(n.e3), as in leia's
    quadricFit model with normalsOnly), the cell included, with the weights
    w_k = |m_k|/sum_j |m_j| (weight area) or alpha_k (1 - alpha_k)
    normalised (weight alpha). The model normal is that of the fitted
    surface at the origin, (-b, -c, 1) rotated back. Ring cells whose
    normal makes n.e3 <= alignTol with the frame are dropped. Fewer than
    minCells usable ring cells (default 6) or a singular system: the
    ringArea normal of the cell, counted as a fit fallback. Default applyTo
    all.

        normalModel ringFit;
        weight      area;     // area | alpha
        alignTol    0.2;
        minCells    6;

SourceFiles
    ringFitNormal.C

\*---------------------------------------------------------------------------*/

#ifndef davofNormalModels_ringFit_H
#define davofNormalModels_ringFit_H

#include "ringAreaNormal.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofNormalModels
{

/*---------------------------------------------------------------------------*\
                            Class ringFit Declaration
\*---------------------------------------------------------------------------*/

class ringFit
:
    public ringArea
{
    // Private Data

        //- area | alpha
        word weight_;

        //- Ring cells with n_k.e3 <= alignTol_ are dropped
        scalar alignTol_;

        //- Fewer usable ring cells: fallback
        label minCells_;


    // Private Member Functions

        //- Orthonormal frame with e3 given
        static void frameOf(const vector& e3, vector& e1, vector& e2);

        //- In-place Cholesky solve of the SPD system A x = b (row-major
        //  n x n), b overwritten with x; false when not positive-definite
        static bool choleskySolve(scalar* A, scalar* b, const label n);


public:

    //- Runtime type information
    TypeName("ringFit");


    //- Construct from the davof dictionary and the scheme
    ringFit
    (
        const dictionary& dict,
        const reconstruction::davof& scheme
    );


    //- Destructor
    virtual ~ringFit() = default;


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
