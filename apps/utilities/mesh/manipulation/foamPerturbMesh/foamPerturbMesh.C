/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) Tomislav Maric and TU Darmstadt
     \\/     M anipulation  |
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

Description
    Perturb every point in a hex-cell mesh.

    1. For each mesh point that is not part of a patch, compute the distance
    to the nearest neighbour vertex. A "neighbour vertex" is another vertex in
    the mesh that is connected to a mesh vertex by an edge.
    2. Move the mesh vertex by alpha (<0.5) of the distance to the
    neighbor into a random direction.

    Note that for alpha >= 0.5 it is theoretically possible to produce self-
    intersecting meshes. This probability increases with increasing alpha.

Author
    Dirk Gründing
    gruending@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis Group
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "argList.H"
#include "fvCFD.H"

#include <fstream>
#include <omp.h>
#include <random>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// produce a random number in [0,1]
scalar normalRand()
{
    // shouldnt there be some randomize seed?
    return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}

autoPtr<pointField> randomUnitVec(label size)
{
    autoPtr<pointField> pVecsPtr(new pointField(size, vector::zero));
    auto& pVecs = pVecsPtr.ref();

    forAll(pVecs, pointI)
    {
        pVecs[pointI].x() = normalRand();
        pVecs[pointI].y() = normalRand();
        pVecs[pointI].z() = normalRand();
    }
    pVecs /= mag(pVecs);
    return pVecsPtr;
}

autoPtr<scalarField> pointBoundingBallRadius(const fvMesh& mesh)
{
    const auto& meshPoints(mesh.points());
    // bbRadius short for boundingBoxRadius
    autoPtr<scalarField> bbRadiusPtr(new scalarField(meshPoints.size(), GREAT));
    auto& bbRadius = bbRadiusPtr.ref();

    forAll(bbRadius, pointI)
    {
        const auto curPoint(meshPoints[pointI]);
        const auto& pointPoints(mesh.pointPoints()[pointI]);
        forAll(pointPoints, ppointI)
        {
            const auto ngbPoint(meshPoints[pointPoints[ppointI]]);
            const auto edgeLength(mag(curPoint - ngbPoint));
            if (edgeLength < mag(bbRadius[pointI]))
            {
                bbRadius[pointI] = edgeLength;
            }
        }
    }
    return bbRadiusPtr;
}

void resetBoundaryPertubations(const fvMesh& mesh, vectorField& pertubations)
{
    forAll(mesh.boundaryMesh(), boundaryI)
    {
        const auto& patch(mesh.boundaryMesh()[boundaryI]); // is a polyPatch
        const auto& patchPointLabels(patch.meshPoints());
        forAll(patchPointLabels, pointI)
        {
            label labelI(patchPointLabels[pointI]);
            pertubations[labelI] = vector::zero;
        }
    }
}

// Main program:
int main(int argc, char* argv[])
{
    argList::addOption("alpha", "scalar < 0.5", "Amount of point pertubation.");

    argList args(argc, argv);

    #include "createTime.H"

    fvMesh mesh(IOobject("region0",
        "constant",
        runTime,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE));

    const auto alpha = args.getOrDefault<scalar>("alpha", 0.1);

    auto pertubationVecPtr = randomUnitVec(mesh.points().size());
    auto maxDistancePtr = pointBoundingBallRadius(mesh);
    tmp<vectorField> tpertubations(
        alpha * maxDistancePtr.ref() * pertubationVecPtr.ref());
    auto pertubations(tpertubations());

    resetBoundaryPertubations(mesh, pertubations);
    tmp<vectorField> new_points(mesh.points() + pertubations);
    mesh.movePoints(new_points());
    mesh.write();

    // Trick the mesh into writing for a time step by writing the same points
    // This allows to compare uncut and cut mesh
    // runTime++;
    auto points(new_points());
    mesh.movePoints(new_points());
    mesh.write();

    Info << "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
