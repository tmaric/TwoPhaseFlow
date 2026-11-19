/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 DLR
    Copyright (C) 2022-2023 OpenCFD Ltd.
    Modified: 2026 Tomislav Maric, MMA TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is derived from OpenFOAM (zoneDistribute, v2506).

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
    Foam::efficientZoneDistribute

Description
    Parallel communication in a narrow band (zoneDistribute of OpenFOAM) with
    neighbour-only exchanges.

    zoneDistribute exchanges through PstreamBuffers::finishedSends with cached
    send/recv connectivity. Every call reduces a "connectivity changed" flag
    over all processors, and the request (setUpCommforZone) and the data
    (getDatafromOtherProc) exchanges flip the send direction, so the cache is
    invalidated and an all-to-all of message sizes follows. Each exchange is
    a global synchronisation point.

    Here both exchanges use PstreamBuffers::finishedNeighbourSends with the
    processors of efficientZoneCPCStencil::neighbourProcs, a symmetric set that
    contains every processor that can own a stencil entry. Sizes and data are
    exchanged with these processors only; there is no global reduction.

SourceFiles
    efficientZoneDistributeI.H
    efficientZoneDistribute.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_efficientZoneDistribute_H
#define Foam_efficientZoneDistribute_H

#include "fvMesh.H"
#include "globalIndex.H"
#include "volFields.H"

#include "efficientZoneCPCStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class efficientZoneDistribute Declaration
\*---------------------------------------------------------------------------*/

class efficientZoneDistribute
:
    public MeshObject<fvMesh, TopologicalMeshObject, efficientZoneDistribute>
{
    // Private Typedefs

        typedef MeshObject
        <
            fvMesh,
            TopologicalMeshObject,
            efficientZoneDistribute
        > MeshObject_type;


    // Private Data

        //- Reference to the zone stencil
        efficientZoneCPCStencil& stencil_;

        //- Global number into index of cells/faces
        const globalIndex& globalNumbering_;

        //- Global cell/face index to send for processor-to-processor comms
        List<labelList> send_;

        //- Requested global indices in the same order sent during setup.
        List<labelList> recv_;

        //- Number of field elements expected from each processor.
        labelList recvSizes_;

        //- Persistent set of exchange buffers
        PstreamBuffers pBufs_;

        //- Cell labels for all cells with a face on a cyclic boundary
        bitSet cyclicBoundaryCells_;

        //- Holds for each cell on cyclic patch the centres of the cells in its
        //  point neighbour stencil.
        //  Used in getCyclicPatches to identify neigbhour patch ID of point
        //  neighbours on other side of a processorPolyPatch.
        autoPtr<Map<vectorField>> cyclicCentres_;


    // Private Member Functions

        //- Exchange the streamed send buffers with the neighbour processors.
        //  Stops if data were streamed to a processor outside the
        //  neighbourhood (PstreamBuffers would drop them silently).
        void finishedNeighbourSends();

        //- Return local volField value at (cell or face) index
        template<typename Type>
        Type getLocalValue
        (
            const VolumeField<Type>& phi,
            const label localIdx
        ) const;

        //- Gives patchNumber and patchFaceNumber for a given
        //- Geometric volume field
        template<typename Type>
        Type faceValue
        (
            const VolumeField<Type>& phi,
            const label localIdx
        ) const;


public:

    //- Runtime information
    TypeName("efficientZoneDistribute");


    // Constructors

        //- Construct from fvMesh
        explicit efficientZoneDistribute(const fvMesh&);

        //- Selector
        static efficientZoneDistribute& New(const fvMesh&);


    //- Destructor
    virtual ~efficientZoneDistribute() = default;


    // Member Functions

        //- Update stencil with boolList the size has to match mesh nCells
        void setUpCommforZone(const boolList& zone, bool updateStencil=true);

        //- Updates stencil with boolList the size has to match mesh nCells
        void updateStencil(const boolList& zone);

        //- Stencil reference
        const labelListList& getStencil() noexcept
        {
            return stencil_;
        }

        //- Addressing reference
        const globalIndex& globalNumbering() const noexcept
        {
            return globalNumbering_;
        }

        //- Finds and returns list of all cyclic patch labels to which celli's
        //  point neighbour cell, globalIdx, belongs. celli and globalIdx touch
        //  in at least one point on these patches. globalIdx typically belongs
        //  to stencil_[celli]. The returned label list is used to transform
        //  positions across cyclic boundaries e.g. to be able to calculate
        //  distances between cell centres and interface centres in plicRDF
        //  across such boundaries.
        List<label> getCyclicPatches
        (
            const label celli,
            const label globalIdx,
            const vector globalIdxCellCentre
        ) const;

        //- Gives patchNumber and patchFaceNumber for a given
        //- Geometric volume field
        template<typename Type>
        Type getValue
        (
            const VolumeField<Type>& phi,
            const Map<Type>& valuesFromOtherProc,
            const label gblIdx
        ) const;

        vector getPosition
        (
            const VolumeField<vector>& positions,
            const Map<vector>& valuesFromOtherProc,
            const label gblIdx,
            const List<label> cyclicPatchID = List<label>()
        ) const;

        //- Returns stencil and provides a Map with globalNumbering
        template<typename Type>
        Map<Field<Type>> getFields
        (
            const boolList& zone,
            const VolumeField<Type>& phi
        );

        //- Returns stencil and provides a Map with globalNumbering
        template<typename Type>
        Map<Field<Type>> getPositionFields
        (
            const boolList& zone,
            const VolumeField<Type>& phi,
            const bool& checkTransformation = false
        );

        //- Returns stencil and provides a Map with globalNumbering
        template<typename Type>
        Map<Type> getDatafromOtherProc
        (
            const boolList& zone,
            const VolumeField<Type>& phi
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "efficientZoneDistributeI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
