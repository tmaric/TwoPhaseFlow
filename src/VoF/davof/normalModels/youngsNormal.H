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
    Foam::davofNormalModels::youngs

Description
    The Youngs normal fvc::grad(alpha) (into alpha = 1 like m_c; gradAlpha's
    convention). Default applyTo check: the slivers and the cells deviating
    by more than maxAngle (Basilisk's mixed Youngs-centred scheme, Aulisa et
    al. 2007, replaces the slivers). A wanted cell where grad(alpha)
    vanishes has no model normal and gets the base class's noNormalFallback
    (keep: m_c; ringArea: the ring normal).

SourceFiles
    youngsNormal.C

\*---------------------------------------------------------------------------*/

#ifndef davofNormalModels_youngs_H
#define davofNormalModels_youngs_H

#include "davofNormalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofNormalModels
{

/*---------------------------------------------------------------------------*\
                            Class youngs Declaration
\*---------------------------------------------------------------------------*/

class youngs
:
    public davofNormalModel
{
public:

    //- Runtime type information
    TypeName("youngs");


    //- Construct from the davof dictionary and the scheme
    youngs
    (
        const dictionary& dict,
        const reconstruction::davof& scheme
    );


    //- Destructor
    virtual ~youngs() = default;


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
