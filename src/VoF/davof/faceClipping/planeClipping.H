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
    Foam::davofFaceClippingModels::plane

Description
    The plane clipping: the stored face fraction is the submerged area of
    the face under the moving PLIC plane at the end of the step, as the
    face update stored it before the models existed. corrects() is false,
    so isoAdvection never calls correct() and keeps no chords: the sweep
    and the stored values are bit-identical with the face update without
    the model (the default of the keyword faceClipping).

SourceFiles
    planeClipping.C

\*---------------------------------------------------------------------------*/

#ifndef davofFaceClippingModels_plane_H
#define davofFaceClippingModels_plane_H

#include "davofFaceClipping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofFaceClippingModels
{

/*---------------------------------------------------------------------------*\
                            Class plane Declaration
\*---------------------------------------------------------------------------*/

class plane
:
    public davofFaceClipping
{
public:

    //- Runtime type information
    TypeName("plane");


    //- Construct from the alpha solver dictionary and the scheme
    plane
    (
        const dictionary& dict,
        const reconstruction::davof& scheme
    );


    //- Destructor
    virtual ~plane() = default;


    // Member Functions

        //- The plane-clipped area, unchanged
        virtual scalar correct
        (
            const label facei,
            const label planeCell,
            const vector& nLiq,
            const point& xEnd,
            const UList<point>& chord,
            const scalar areaPlane
        )
        {
            return areaPlane;
        }

        //- Nothing is reported for the plane clipping
        virtual void report(Ostream& os) const
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace davofFaceClippingModels
} // End inline namespace twoPhaseFlow
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
