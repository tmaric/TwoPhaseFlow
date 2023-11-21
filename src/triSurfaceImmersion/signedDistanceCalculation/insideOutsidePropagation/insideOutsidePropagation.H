/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 Tobias Tolle, TU Darmstadt
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
    Foam::TriSurfaceImmersion::insideOutsidePropagation

Author
    Idea:
    Tomislav Maric
    maric@mma.tu-darmstadt.de
    Mathematical Modelling and Analysis Group
    Thermo-Fluids and Interfaces
    TU Darmstadt
    Germany

    Implementation:
    Tobias Tolle
    tolle@mma.tu-darmstadt.de
    Mathematical Modelling and Analysis Group
    Thermo-Fluids and Interfaces
    TU Darmstadt
    Germany

Description
    Propagate inside / outside information from a narrow band around an
    interface to the entire domain by solving a Laplace equation.
    This requires a narrow band thickness of 3 cells on each side of the
    interface. Furthermore, the signed distance for cells not in the narrow
    band must be set to zero.

SourceFiles
    insideOutsidePropagation.C

\*---------------------------------------------------------------------------*/

#ifndef insideOutsidePropagation_H
#define insideOutsidePropagation_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion
{
/*---------------------------------------------------------------------------*\
                         Class insideOutsidePropagation Declaration
\*---------------------------------------------------------------------------*/

class insideOutsidePropagation
{
public:
    // Member Functions
    // - Propagate inside/outside information
    //   Propagate inside/outside information of the given field to the
    //   entire domain. The value for cells without a valid signed distance
    //   needs to be 0.0.
    static tmp<volScalarField> propagateInsideOutside(
        const volScalarField& signedDistance);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
