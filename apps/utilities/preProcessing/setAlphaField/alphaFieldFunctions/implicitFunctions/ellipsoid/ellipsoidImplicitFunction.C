/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2019-2020 DLR
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

#include "ellipsoidImplicitFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace implicitFunctions
    {
        defineTypeNameAndDebug(ellipsoidImplicitFunction, 0);
        addToRunTimeSelectionTable
        (
            implicitFunction,
            ellipsoidImplicitFunction,
            dict
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::implicitFunctions::ellipsoidImplicitFunction::ellipsoidImplicitFunction
(
    const vector& semiAxis,
    const vector& origin
)
:
    semiAxis_(semiAxis),
    origin_(origin)
{}


Foam::implicitFunctions::ellipsoidImplicitFunction::ellipsoidImplicitFunction
(
    const dictionary& dict
)
:
    ellipsoidImplicitFunction
    (
        [&dict]()
        {
            if (dict.found("semiAxis"))
            {
                return dict.get<vector>("semiAxis");
            }

            const scalar equivalentRadius =
                dict.get<scalar>
                (
                    dict.found("effectiveRadius")
                  ? "effectiveRadius"
                  : "equivRadius"
                );

            const scalar horizontalLongAxis =
                dict.get<scalar>("horizontalLongAxis");

            const scalar horizontalSemiAxis = horizontalLongAxis;
            const scalar verticalSemiAxis =
                pow3(equivalentRadius)/sqr(horizontalSemiAxis);

            return vector
            (
                horizontalSemiAxis,
                verticalSemiAxis,
                horizontalSemiAxis
            );
        }(),
        dict.getOrDefault<vector>("origin", Zero)
    )
{}


// ************************************************************************* //
