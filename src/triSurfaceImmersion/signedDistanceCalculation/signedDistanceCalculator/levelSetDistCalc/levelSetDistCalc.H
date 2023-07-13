/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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
    Foam::TriSurfaceImmersion::levelSetDistCalc

Description
    Compute the signed distance to a surface implicitly given by a level
    set.

SourceFiles
    levelSetDistCalcI.H
    levelSetDistCalc.cpp

\*---------------------------------------------------------------------------*/

#ifndef levelSetDistCalc_H
#define levelSetDistCalc_H

#include "signedDistanceCalculator.H"

#include "implicitSurfaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion
{
/*---------------------------------------------------------------------------*\
                         Class levelSetDistCalc Declaration
\*---------------------------------------------------------------------------*/

class levelSetDistCalc
:
    public signedDistanceCalculator
{
    // Private Data
    //- Pointer to the implicit surface
    autoPtr<implicitSurface> surfacePtr_;

    //- Maximum number of iterations for gradient descent and closest point
    label maxIt_;

    //- Error threshold: below this value algorithms are considered converged
    scalar epsilon_;

    //- Field of level set values at cell centres
    volScalarField cellLevelSetValues_;

    //- Field of level set values at cell corners
    pointScalarField pointLevelSetValues_;

    //- Labels of all cells in the narrow band
    std::vector<label> narrowBandCells_;


    // Private Member Functions
    //- Return a point on the surface using gradient descend starting at p
    point surfacePoint(const point& p) const;

    //- Find the/a closest point on the surface to p
    point closestPoint(const point& p) const;

    //- Compute level set values at cell centres and cell corners
    void computeMeshLevelSetValues();

    //- Find all cells constituting the narrow band
    void identifyNarrowBandCells();

    //- Compute signed distance at cell centres and corners in the narrow band
    void computeSignedDistances();

    //- Set inside/outside information for points outside of the narrow band
    void setInsideOutside();

public:

    // Static Data Member
    TypeName("levelSet");


    // Constructors
    levelSetDistCalc(const dictionary& configDict, const fvMesh& mesh);


    // Member Functions
    //- Return signed distance of point x to surface
    scalar signedDistance(const point& x) const override;

    //- Return characteristic length of surface
    scalar referenceLength() const override;

    //- Returns 1. Underlyiung surface is not composed of discrete elements.
    label nSurfaceElements() const override;

    //- Returns enclosed volume if surface is closed.
    //  Returns negative value if surface is not closed.
    scalar surfaceEnclosedVolume() const override;

    //- Reference to implicit surface
    inline const implicitSurface& surface() const;

    //- Return the maximum number of iterations for algorithms.
    inline label maxIter() const;

    //- Return the convergence error threshold.
    inline scalar epsilon() const;

    //- Return a reference to the cell centred level set values.
    inline const volScalarField& cellLevelSetValues() const;

    //- Return a reference to the cell corner level set values.
    inline const pointScalarField& pointLevelSetValues() const;

    //- Write the level set fields
    void writeFields() const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "levelSetDistCalcI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
