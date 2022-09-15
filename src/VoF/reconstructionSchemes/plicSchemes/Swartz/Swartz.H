/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 TU Darmstadt 
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
    Foam::reconstruction::Swartz

Description
    Reconstructs an interface (centre and normal vector) consisting of planes
    to match the internal fluid distribution in cells. 

    Reference:
    \verbatim
        Marić, T., Marschall, H., & Bothe, D. (2018). 
        An enhanced un-split face-vertex flux-based VoF method. 
        Journal of Computational Physics, 371, 967–993.
        https://doi.org/10.1016/j.jcp.2018.03.048
    \endverbatim

Authors
    Tomislav Maric, TU Darmstadt (2022)

SourceFiles
    Swartz.H
    Swartz.C

\*---------------------------------------------------------------------------*/

#ifndef Swartz_H
#define Swartz_H

#include "typeInfo.H"
#include "reconstructionSchemes.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "surfaceIteratorPLIC.H"
#include "reconstructedDistanceFunction.H"
#include "zoneDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{

/*---------------------------------------------------------------------------*\
                           Class Swartz Declaration
\*---------------------------------------------------------------------------*/

class Swartz
:
    public reconstructionSchemes
{

    // Private Data
    
    // Number of iterations
    label iteration_;  

    // Reconstruction tolerance.
    scalar isoFaceTol_; // FIXME: shared with plicRDF, move upstream. TM.

    // Critical angle between two normal vectors.
    scalar critAngle_; 

    // Youngs initial normal field. 
    volVectorField youngsNormals_;
    
    // Improved normal field
    volVectorField swartzNormals_; 

    // List of MPI process patches.
    DynamicList<label> procPatchIds_;

    // Interface positioning iterator
    // TODO: Profile + discuss CCS positioning. TM. 
    surfaceIteratorPLIC sIterPLIC_;

public:

    //- Runtime type information
    TypeName("Swartz");

    //- Construct from components
    Swartz
    (
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U,
        const dictionary& dict
    );

    //- Destructor
    virtual ~Swartz() = default;

    // Member Functions

        //- Improve normal orientation 
        virtual void improveNormals(); 

        //- Position the interface.
        virtual void positionInterface(); 

        //- Reconstruct interface
        virtual void reconstruct(bool forceUpdate = true);

        //- Map VoF Field in case of refinement
        virtual void mapAlphaField() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reconstruction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
