/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 DLR
    Modified: 2026 Tomislav Maric, MMA TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is derived from OpenFOAM (zoneCPCStencil, v2506).

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
    Foam::efficientZoneCPCStencil

Description
    Cell-point-cell stencil of the cells of a zone (zoneCPCStencil of OpenFOAM)
    without parallel communication per stencil update.

    The neighbours of a boundary point (the cells and valid boundary faces
    around it, on all processors that share the point) depend only on the mesh,
    not on the zone. zoneCPCStencil synchronises them for the zone points at
    every update, with variable-size lists whose exchange includes an
    all-to-all of message sizes over all processors. Here they are computed and
    synchronised once per mesh, and a stencil update is purely local.

    The same lists give the processors whose cells or boundary faces can appear
    in a stencil of this processor (neighbourProcs). The relation is symmetric,
    so it can be used for neighbour-only exchanges (efficientZoneDistribute).

    The stencil is a TopologicalMeshObject: it is deleted on a topology change
    (including a redistribution) and rebuilt on the next use.

SourceFiles
    efficientZoneCPCStencil.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_efficientZoneCPCStencil_H
#define Foam_efficientZoneCPCStencil_H

#include "zoneCellStencils.H"
#include "boolList.H"
#include "Map.H"
#include "fvMesh.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class efficientZoneCPCStencil Declaration
\*---------------------------------------------------------------------------*/

class efficientZoneCPCStencil
:
    public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        efficientZoneCPCStencil
    >,
    public zoneCellStencils
{
    // Private Typedefs

        typedef MeshObject
        <
            fvMesh,
            TopologicalMeshObject,
            efficientZoneCPCStencil
        > MeshObject_type;


    // Private Data

        //- Is valid boundary face? (not on empty or coupled patches)
        boolList isValidBFace_;

        //- Cell stencil is already up to date
        boolList uptodate_;

        //- Neighbours (global cells and valid boundary faces) of every
        //- non-empty boundary point, from all processors sharing the point
        Map<labelList> boundaryPointGlobals_;

        //- Processors that can contribute to a stencil of this processor
        labelList neighbourProcs_;


    // Private Member Functions

        //- Compute and synchronise boundaryPointGlobals_, set neighbourProcs_
        void calcBoundaryPointGlobals();

        //- Calculate the stencil of the zone cells that are not up to date
        virtual void calculateStencil
        (
            const boolList& zone,
            labelListList& globalCellCells
        );

        //- No copy construct
        efficientZoneCPCStencil(const efficientZoneCPCStencil&) = delete;

        //- No copy assignment
        void operator=(const efficientZoneCPCStencil&) = delete;


public:

    //- Runtime information
    TypeName("efficientZoneCPCStencil");


    // Constructors

        //- Construct from mesh
        explicit efficientZoneCPCStencil(const fvMesh& mesh);


    // Selectors

        static efficientZoneCPCStencil& New(const fvMesh& mesh);


    // Member Functions

        //- Processors whose cells or boundary faces can appear in a stencil
        //- of this processor (sorted, symmetric relation)
        const labelList& neighbourProcs() const noexcept
        {
            return neighbourProcs_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
