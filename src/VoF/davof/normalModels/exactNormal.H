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
    Foam::davofNormalModels::exact

Description
    The control: the exact normal of the implicit function of the
    exactNormal sub-dictionary (functionType, origin at t = 0, U0, the
    function's keys, as in areaFractionError) translated by U0 t, i.e. grad
    value at the centroid of the plane positioned by alpha along m_c (the
    cell centre when that plane does not exist), at the time tState the
    state belongs to (davof::stateTime). Default applyTo all.

        exactNormal
        {
            functionType sphere;
            origin       (0.5 0.5 0.5);   // at t = 0
            radius       0.15;            // and the function's other keys
            U0           (1 0 0);
        }

SourceFiles
    exactNormal.C

\*---------------------------------------------------------------------------*/

#ifndef davofNormalModels_exact_H
#define davofNormalModels_exact_H

#include "davofNormalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofNormalModels
{

/*---------------------------------------------------------------------------*\
                            Class exact Declaration
\*---------------------------------------------------------------------------*/

class exact
:
    public davofNormalModel
{
    // Private Data

        //- The implicit function's dictionary and its keys
        dictionary exactNormalDict_;
        word functionType_;
        vector origin0_;
        vector U0_;


public:

    //- Runtime type information
    TypeName("exact");


    //- Construct from the davof dictionary and the scheme
    exact
    (
        const dictionary& dict,
        const reconstruction::davof& scheme
    );


    //- Destructor
    virtual ~exact() = default;


    // Member Functions

        virtual void modelNormals(const scalar tState, vectorField& n);

        //- The exact normal is taken at the centroid of the plane
        //  positioned by alpha along m_c
        virtual bool needsPositionedPlanes() const
        {
            return true;
        }

        const word& functionType() const noexcept
        {
            return functionType_;
        }

        const vector& origin0() const noexcept
        {
            return origin0_;
        }

        const vector& U0() const noexcept
        {
            return U0_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace davofNormalModels
} // End inline namespace twoPhaseFlow
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
