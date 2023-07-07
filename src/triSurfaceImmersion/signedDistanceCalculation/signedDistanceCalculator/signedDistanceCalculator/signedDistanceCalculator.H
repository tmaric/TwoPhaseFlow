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
    Foam::TriSurfaceImmersion::signedDistanceCalculator

Description
    Unified interface class for signed distance calculation based on
    triangulated surfaces and implicit surfaces defined by a level set.

SourceFiles
    signedDistanceCalculatorI.H
    signedDistanceCalculator.cpp

\*---------------------------------------------------------------------------*/

#ifndef signedDistanceCalculator_H
#define signedDistanceCalculator_H

#include "DynamicList.H"
#include "pointFields.H"
#include "pointIndexHit.H"
#include "pointMesh.H"
#include "runTimeSelectionTables.H"
#include "typeInfo.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion
{

/*---------------------------------------------------------------------------*\
                    Class signedDistanceCalculator Declaration
\*---------------------------------------------------------------------------*/

class signedDistanceCalculator
{
private:
    // Private Data
    //- Dictionary used to construct the object
    const dictionary& dict_;

    //- Reference to the underlying mesh
    const fvMesh& mesh_;

    //- Point mesh formed by the mesh corner points
    pointMesh pMesh_;

    //- Width of the narrow band om each side of interface in number of cells
    scalar narrowBandWidth_;

    //- Value used for points outside of the narrow band
    scalar outOfNarrowBandValue_;


protected:
    //- Signed distance at cell centers.
    volScalarField cellSignedDist_;

    //- Initial signed distance field given by the octree, used to correct the
    //  signed distance propagated by the solution of the Laplace equation.
    volScalarField cellSignedDist0_;

    //- Nearest surface triangle to a cell centre in the narrow band
    DynamicList<pointIndexHit> cellNearestTriangle_;

    //- Signed distance at cell corner points.
    pointScalarField pointSignedDist_;

    //- Nearest surface triangle to a cell centre in the narrow band
    DynamicList<pointIndexHit> pointNearestTriangle_;


public:
    // Static Data Members
    TypeName("signedDistanceCalculatorInterface");

    declareRunTimeSelectionTable(autoPtr,
        signedDistanceCalculator,
        Dictionary,
        (const dictionary& configDict, const fvMesh& mesh),
        (configDict, mesh));


    // Constructors
    explicit signedDistanceCalculator(
        const dictionary& configDict, const fvMesh& mesh);


    // Selectors
    static autoPtr<signedDistanceCalculator> New(
        const dictionary& configDict, const fvMesh& mesh);


    // Member Functions
    //- Acces to configuration dictionary
    inline const dictionary& configDict() const;

    //- Access to underlying mesh
    inline const fvMesh& mesh() const;

    //- Access to point mesh formed by mesh corner points
    inline const pointMesh& pMesh() const;

    //- Width of narrow band in number of cells
    inline scalar narrowBandWidth() const;

    //- Value used for points outside of the narrow band
    inline scalar outOfNarrowBandValue() const;

    //- Map giving the closest surface point for each cell centre
    inline const DynamicList<pointIndexHit>& cellClosestPoint() const;

    //- Map giving the closest surface point for each cell corner point
    inline const DynamicList<pointIndexHit>& pointClosestPoint() const;

    //- Signed distances at cell centres with bulk inside/outside information
    inline const volScalarField& cellSignedDist() const;

    //- Narrow band signed distances at cell centres
    inline const volScalarField& cellSignedDist0() const;

    //- Narrow band signed distances at cell corner points
    inline const pointScalarField& pointSignedDist() const;

    //- Return the signed distance of point x to the interface
    virtual scalar signedDistance(const point& x) const = 0;

    //- Return a characteristic length of the surface
    virtual scalar referenceLength() const = 0;

    //- Return the number of surface elements.
    //  For a triangulated surface, this is the number of triangles.
    //  For an implicit surface, it returns 1.
    virtual label nSurfaceElements() const = 0;

    //- Returns enclosed volume if surface is closed.
    //  If the surface is not closed, a negative value is returned.
    virtual scalar surfaceEnclosedVolume() const = 0;

    //- Set the value to use for points outside the narrow band
    void outOfNarrowBandValue(scalar value);

    //- Set the width of the narrow band as number of cells
    void narrowBandWidth(scalar width);

    //- Write the signed distance fields
    virtual void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "signedDistanceCalculatorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
