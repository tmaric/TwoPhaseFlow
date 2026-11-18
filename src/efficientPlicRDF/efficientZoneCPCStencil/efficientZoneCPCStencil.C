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

\*---------------------------------------------------------------------------*/

#include "efficientZoneCPCStencil.H"
#include "syncTools.H"
#include "ListOps.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(efficientZoneCPCStencil, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::efficientZoneCPCStencil::calcBoundaryPointGlobals()
{
    const fvMesh& mesh = meshRef_;

    // Local neighbours of every non-empty boundary point
    const labelList boundaryPoints(nonEmptyFacesPatch()().meshPoints());

    boundaryPointGlobals_.reserve(boundaryPoints.size());

    labelHashSet pointGlobals;

    for (const label pointi : boundaryPoints)
    {
        boundaryPointGlobals_.insert
        (
            pointi,
            calcFaceCells
            (
                isValidBFace_,
                mesh.pointFaces()[pointi],
                pointGlobals
            )
        );
    }

    // Union over the processors sharing a coupled point. Once per mesh, so
    // the exchange of variable-size lists (all-to-all of sizes) is paid once.
    const labelList& coupledPoints =
        mesh.globalData().coupledPatch().meshPoints();

    List<labelList> coupledGlobals(coupledPoints.size());

    forAll(coupledPoints, i)
    {
        const auto iter = boundaryPointGlobals_.cfind(coupledPoints[i]);

        if (iter.good())
        {
            coupledGlobals[i] = iter.val();
        }
    }

    syncTools::syncPointList
    (
        mesh,
        coupledPoints,
        coupledGlobals,
        ListOps::unionEqOp(),
        labelList(),
        Foam::dummyTransform()
    );

    forAll(coupledPoints, i)
    {
        boundaryPointGlobals_.set
        (
            coupledPoints[i],
            std::move(coupledGlobals[i])
        );
    }

    // Processors owning an entry. A coupled point carries the same union on
    // every processor that shares it, so the relation is symmetric.
    labelHashSet procs;

    forAllConstIters(boundaryPointGlobals_, iter)
    {
        for (const label gblIdx : iter.val())
        {
            const label proci = globalNumbering().whichProcID(gblIdx);

            if (proci != UPstream::myProcNo())
            {
                procs.insert(proci);
            }
        }
    }

    neighbourProcs_ = procs.sortedToc();
}


void Foam::efficientZoneCPCStencil::calculateStencil
(
    const boolList& zone,
    labelListList& globalCellCells
)
{
    const fvMesh& mesh = meshRef_;
    const labelListList& cPoints = mesh.cellPoints();
    const labelListList& pCells = mesh.pointCells();

    labelHashSet stencilSet;

    forAll(zone, celli)
    {
        if (!zone[celli] || uptodate_[celli])
        {
            continue;
        }

        const label globalI = globalNumbering().toGlobal(celli);

        stencilSet.clear();

        bool needsComm = false;

        for (const label pointi : cPoints[celli])
        {
            // Local point neighbours
            for (const label neiCelli : pCells[pointi])
            {
                stencilSet.insert(globalNumbering().toGlobal(neiCelli));
            }

            // Boundary point: neighbours on other processors and valid
            // boundary faces
            const auto iter = boundaryPointGlobals_.cfind(pointi);

            if (iter.good())
            {
                for (const label gblIdx : iter.val())
                {
                    stencilSet.insert(gblIdx);

                    if (!globalNumbering().isLocal(gblIdx))
                    {
                        needsComm = true;
                    }
                }
            }
        }

        stencilSet.erase(globalI);

        // First element is always the cell itself
        labelList& cCells = globalCellCells[celli];
        cCells.resize_nocopy(stencilSet.size() + 1);

        label n = 0;
        cCells[n++] = globalI;

        for (const label gblIdx : stencilSet)
        {
            cCells[n++] = gblIdx;
        }

        if (needsComm)
        {
            needComm_.insert(celli);
        }

        uptodate_[celli] = true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::efficientZoneCPCStencil::efficientZoneCPCStencil(const fvMesh& mesh)
:
    MeshObject_type(mesh),
    zoneCellStencils(mesh),
    uptodate_(mesh.nCells(), false)
{
    // Mark boundary faces to be included in stencil (i.e. not coupled or empty)
    validBoundaryFaces(isValidBFace_);

    calcBoundaryPointGlobals();
}


Foam::efficientZoneCPCStencil& Foam::efficientZoneCPCStencil::New
(
    const fvMesh& mesh
)
{
    auto* ptr =
        mesh.thisDb().getObjectPtr<efficientZoneCPCStencil>
        (
            efficientZoneCPCStencil::typeName
        );

    if (!ptr)
    {
        ptr = new efficientZoneCPCStencil(mesh);
        regIOobject::store(ptr);
    }

    return *ptr;
}


// ************************************************************************* //
