/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 DLR
    Modified: 2026 Tomislav Maric, MMA TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is derived from OpenFOAM (plicRDF, v2506).

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
    Foam::reconstruction::efficientPlicRDF

Description
    Reconstructs an interface (centre and normal vector) consisting of planes
    to match the internal fluid distribution in cells. The interface normals
    are estimated by least square gradient scheme on the RDF function (height).
    Uses the normal from the previous times step as intial guess.

    Reference:
    \verbatim
        Henning Scheufler, Johan Roenby,
        Accurate and efficient surface reconstruction from volume
        fraction data on general meshes,
        Journal of Computational Physics, 2019,
        doi 10.1016/j.jcp.2019.01.009

    \endverbatim

    Original code supplied by Henning Scheufler, DLR (2019)

    efficientPlicRDF is OpenFOAM's plicRDF with the same algorithm and the
    same settings, and parallel communication that does not scale with the
    number of processors:
    - the cell-point-cell stencil of the interface band is updated without
      communication (efficientZoneCPCStencil: the neighbours of boundary
      points are synchronised once per mesh),
    - stencil data are exchanged with the stencil neighbour processors only,
      without a global reduction per exchange (efficientZoneDistribute),
    - the band marking synchronises coupled points through globalMeshData
      instead of a reduction on the master
      (efficientReconstructedDistanceFunction).
    Results agree with plicRDF up to the order of the stencil entries
    (round-off).

    Usage, system/fvSolution and system/controlDict:
    \verbatim
    "alpha.water.*"
    {
        advectionScheme      isoAdvection;
        reconstructionScheme efficientPlicRDF;
        ...
    }

    libs (efficientPlicRDF);
    \endverbatim
    Coefficients (isoFaceTol, surfCellTol, tol, relTol, iterations,
    interpolateNormal) as for plicRDF, optionally in efficientPlicRDFCoeffs.

SourceFiles
    efficientPlicRDF.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_efficientPlicRDF_H
#define Foam_efficientPlicRDF_H

#include "typeInfo.H"
#include "reconstructionSchemes.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "surfaceIteratorPLIC.H"
#include "efficientReconstructedDistanceFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{

/*---------------------------------------------------------------------------*\
                           Class efficientPlicRDF Declaration
\*---------------------------------------------------------------------------*/

class efficientPlicRDF
:
    public reconstructionSchemes
{
    // Private Data

        //-  residuals storage
        struct normalRes
        {
            label celli = {};
            scalar normalResidual = {};
            scalar avgAngle = {};
        };

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Interpolation object from cell centres to points
        DynamicField<vector> interfaceNormal_;


    // Switches and tolerances. Tolerances need to go into toleranceSwitches

        //- Tolerance for search of isoFace giving specified VOF value
        scalar isoFaceTol_;

        //- Tolerance for marking of surface cells:
        //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
        scalar surfCellTol_;

        //- Tolerance
        scalar tol_;

        //- Relative tolerance
        scalar relTol_;

        //- Number of times that the interface is reconstructed
        //- has to be bigger than 2
        label iteration_;

        //- Interpolated normal from previous time step
        bool interpolateNormal_;

        //- Calculates the RDF function
        efficientReconstructedDistanceFunction RDF_;

        //- surfaceIterator finds the plane centre for specified VOF value
        surfaceIteratorPLIC sIterPLIC_;


    // Private Member Functions

        //- Set initial normals by interpolation from the previous
        //- timestep or with the Young method
        void setInitNormals(bool interpolate);

        //- compute gradient at the surfaces
        void gradSurf(const volScalarField& phi);

        //- compute the normal residuals
        void calcResidual
        (
            List<normalRes>& normalResidual
        );

        //- interpolation of the normals from the previous time step
        void interpolateNormal();

        //- No copy construct
        efficientPlicRDF(const efficientPlicRDF&) = delete;

        //- No copy assignment
        void operator=(const efficientPlicRDF&) = delete;


public:

    //- Runtime type information
    TypeName("efficientPlicRDF");


    //- Construct from components
    efficientPlicRDF
    (
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U,
        const dictionary& dict
    );


    //- Destructor
    virtual ~efficientPlicRDF() = default;


    // Member Functions

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
