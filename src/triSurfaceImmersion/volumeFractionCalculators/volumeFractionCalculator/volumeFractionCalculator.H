/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
 
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
    Foam::TriSurfaceImmersion::volumeFractionCalculator
 
Description
    Interface class enabling runtime type selection for volume fraction
    calculators, namely surface-mesh-cell-intersection (SMCI) and
    surface-mesh-cell-approximation (SMCA).
 
SourceFiles
    volumeFractionCalculator.cpp
    volumeFractionCalculatorI.H
 
Authors (in alphabetical order):
 
    Dirk Gründing   (gruending@mma.tu-darmstadt.de, main developer, Affiliation 2,3)  
    Tomislav Maric  (tomislav.maric@gmx.com, main developer, Affiliation 1,3)
    Tobias Tolle    (bt@lefou-familie.org, main developer, Affiliation 1,3)
 
Additional Contact persons:
 
    Dieter Bothe, bothe@mma.tu-darmstadt.de, head of institute, Affiliation 3
 
Affiliations:
 
    1) Research Group:
       Lagrangian / Eulerian numerical methods for multiphase flows
       Research Group Leader: Tomislav Maric
       Department of Mathematics
       Technische Universität Darmstadt, Germany
 
    2) Research Group:
       Arbitrary Lagrangian Eulerian methods for multiphase flows and
       wetting phenomena
       Research Group Leader: Dirk Gründing
       Department of Mathematics
       Technische Universität Darmstadt, Germany
 
    3) Institute:
       Mathematical modeling and analysis
       Head of institute: Dieter Bothe
       Department of Mathematics
       Technische Universität Darmstadt, Germany

Acknowledgement:
 
    Funded by the German Research Foundation (DFG):
    Project-ID 265191195 – SFB 1194, sub-projects B01, B02 and Z-INF
     
\*---------------------------------------------------------------------------*/

#ifndef volumeFractionCalculator_H
#define volumeFractionCalculator_H

#include "Time.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "pointFields.H"
#include "pointIndexHit.H"
#include "runTimeSelectionTables.H"
#include "surfaceFields.H"
#include "triSurface.H"
#include "typeInfo.H"
#include "volFields.H"

#include "signedDistanceCalculator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion
{
/*---------------------------------------------------------------------------*\
                 Class volumeFractionCalculator Declaration
\*---------------------------------------------------------------------------*/

class volumeFractionCalculator
{
private:

    // Private Data
    //- A reference to the mesh
    const fvMesh& mesh_;

    //- A reference to time
    const Time& runTime_;

    //- Point mesh required to construct point fields
    pointMesh pMesh_;

    //- Write additional geometric objects
    //  Write additional geometric objects used during the computation of
    //  volume fractions in interface cells.
    const bool writeGeometry_;


public:

    // Static Data Members
    TypeName("volumeFractionCalculatorInterface");

    declareRunTimeSelectionTable(autoPtr,
        volumeFractionCalculator,
        Dictionary,
        (const dictionary& configDict, const fvMesh& mesh),
        (configDict, mesh));


    // Constructors
    explicit volumeFractionCalculator(
        const dictionary& configDict, const fvMesh& mesh);


    // Selectors
    static autoPtr<volumeFractionCalculator> New(
        const dictionary& configDict, const fvMesh& mesh);


    // Member Functions
    //- Return reference to time.
    inline const Time& time() const;

    //- Return reference to the underlying volume mesh.
    inline const fvMesh& mesh() const;

    //- Return whether additional geometric objects are written.
    inline bool writeGeometry() const;

    //- Return the average number of triangles per interface cell.
    virtual double nTrianglesPerCell() const = 0;

    //- Return the number of intersected cells (interface cells).
    virtual label nIntersectedCells() const = 0;

    //- Return the maximum refinement level used for tetrahedral refinement.
    virtual label maxRefinementLevel() const = 0;

    //- Return reference to the signed distance calculator.
    virtual const signedDistanceCalculator& sigDistCalc() const = 0;

    //- Identify cells intersected by the interface.
    virtual void findIntersectedCells() = 0;

    //- Compute volume fractions according to signed distance at cell centre.
    //  Sets volume fractions according to the sign of the signed distance
    //  at cell centres for all cells including interface cells.
    void bulkVolumeFraction(volScalarField& alpha);

    //- Compute area fractions according to interpolated signed distance from cell centres.
    //  Sets area fractions according to the sign of the signed distance
    //  at face centres for all faces including those of interface cells.
    void bulkAreaFraction(surfaceScalarField& alpha);

    //- Compute the volume fraction field alpha.
    virtual void calcVolumeFraction(volScalarField& alpha) = 0;

    //- Compute the area fraction field alpha.
    virtual void calcAreaFraction(surfaceScalarField& alpha) = 0;

    //- Write out additional fields used by the volume fraction calculator.
    virtual void writeFields() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volumeFractionCalculatorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
