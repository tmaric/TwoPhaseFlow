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

Class
    Foam::davofVolumeRules::exactSphere

Description
    The exact overlap volume of a sphere with tetrahedra,

        S. Strobl, A. Formella, T. Poeschel, Exact calculation of the
        overlap volume of spheres and mesh elements, Journal of
        Computational Physics 311 (2016) 158-172,
        doi:10.1016/j.jcp.2016.02.003,

    through the authors' reference implementation (thirdParty/overlap/
    overlap.hpp, MIT licence, https://github.com/severinstrobl/overlap,
    header-only on Eigen 3): the cell as the tetrahedra (cell centre, face
    centre, edge) of every face's fan, the prism between a face's pre-image
    and the face as three tetrahedra per fan triangle (the triangular prism
    split along a diagonal of each side). The face fraction is the quadratic
    face cut of davofCurvedGeometry.H, exact for the circle trace up to the
    sliver's 12-point Gauss rule. The rule for the sphere translation case
    (clippingSurface exactSphere); a paraboloid surface falls back to the
    split quadrature (counted).

SourceFiles
    exactSphereVolumeRule.C

\*---------------------------------------------------------------------------*/

#ifndef exactSphereVolumeRule_H
#define exactSphereVolumeRule_H

#include "davofVolumeRule.H"

namespace Foam
{
namespace davofVolumeRules
{

/*---------------------------------------------------------------------------*\
                         Class exactSphere Declaration
\*---------------------------------------------------------------------------*/

class exactSphere
:
    public davofVolumeRule
{
    // Private Member Functions

        //- The overlap of the sphere S with the tetrahedron (a, b, c, d)
        scalar tetOverlap
        (
            const davof::cellSurface& S,
            const point& a, const point& b, const point& c, const point& d
        ) const;


public:

    //- Runtime type information
    TypeName("exactSphere");


    // Constructors

        exactSphere(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~exactSphere() = default;


    // Member Functions

        virtual scalar cellLiquidVolume
        (
            const label celli,
            const davof::cellSurface& S,
            const scalar tolRel
        ) const;

        virtual scalar prismLiquidVolume
        (
            const label facei,
            const vector& d,
            const davof::cellSurface& S,
            const scalar tolRel
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace davofVolumeRules
} // End namespace Foam

#endif

// ************************************************************************* //
