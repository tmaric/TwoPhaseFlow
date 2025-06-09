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

\*---------------------------------------------------------------------------*/

#include "AdaptiveTetCellRefinement.H"
#include "surfaceFieldsFwd.H"
#include "surfaceMeshCellApproximation.H"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <tuple>
#include <volFieldsFwd.H>

#include "addToRunTimeSelectionTable.H"

#include "IntersectionCriteria.H"
#include "signedDistanceCalculator.H"
#include "tetVolumeFractionCalculator.H"
#include "triSurfaceDistCalc.H"

namespace Foam::TriSurfaceImmersion
{

defineTypeNameAndDebug(surfaceMeshCellApproximation, 0);
addToRunTimeSelectionTable(
    volumeFractionCalculator, surfaceMeshCellApproximation, Dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
surfaceMeshCellApproximation::cellDecompositionTuple
surfaceMeshCellApproximation::decomposeCell(const label cellID) const
{
    const auto& mesh = this->mesh();
    const auto& thisCell = mesh.cells()[cellID];
    const auto& cellVertexIDs = mesh.cellPoints()[cellID];
    const auto& vertices = mesh.points();
    const auto& pointSignedDist = sigDistCalcPtr_->pointSignedDist();
    const auto& cellSignedDist = sigDistCalcPtr_->cellSignedDist0();

    std::vector<indexedTet> tets(nTets(cellID));

    // Using a barycentric decomposition, the number of unique points
    // is the sum of n_cell_vertices + n_cell_faces + 1 (the cell centre) (TT).
    std::vector<point> points(cellVertexIDs.size() + thisCell.size() + 1);
    std::vector<scalar> sd(points.size());
    std::map<label, label> globalToLocal{};

    // Add vertices to points and their signed distance
    forAll(cellVertexIDs, I)
    {
        points[I] = vertices[cellVertexIDs[I]];
        sd[I] = pointSignedDist[cellVertexIDs[I]];
        globalToLocal[cellVertexIDs[I]] = I;
    }

    // Add the cell centre
    label centre_id = cellVertexIDs.size();
    points[centre_id] = mesh.C()[cellID];
    sd[centre_id] = cellSignedDist[cellID];

    // Add face centres and build the indexed tets
    const auto& faces = mesh.faces();
    label face_centre_id = centre_id + 1;
    label idx_tet = 0;
    for (const auto face_id : thisCell)
    {
        points[face_centre_id] = mesh.Cf()[face_id];
        sd[face_centre_id] =
            sigDistCalcPtr_->signedDistance(mesh.Cf()[face_id]);

        for (const auto& anEdge : faces[face_id].edges())
        {
            tets[idx_tet] = indexedTet{centre_id,
                face_centre_id,
                globalToLocal[anEdge[0]],
                globalToLocal[anEdge[1]]};
            ++idx_tet;
        }
        ++face_centre_id;
    }

    // Signed distance plausibility check
    for (uint idx = 0; idx != points.size(); ++idx)
    {
        assert(mag(sd[idx]) < mesh.bounds().mag());
    }

    return std::make_tuple(tets, points, sd);
}


surfaceMeshCellApproximation::faceDecompositionTuple
surfaceMeshCellApproximation::decomposeFace(label faceID) const
{
    const auto& mesh = this->mesh();
    const auto& vertices = mesh.points();
    const auto& thisFace = mesh.faces()[faceID];
    const auto& pointSignedDist = sigDistCalcPtr_->pointSignedDist();

    // Using a barycentric decomposition, the number of unique points
    // is the number of face vertices + 1 for the face centroid. (TT).
    std::vector<indexedTri> tris(thisFace.size());
    std::vector<point> points(thisFace.size() + 1);
    std::vector<scalar> sd(points.size());

    // TODO: the vertex indices for a face should already be sorted,
    // so two consecutive vertices form an edge.
    // Yet, this has to be verified (TT).
    forAll(thisFace, I)
    {
        label pointID = thisFace[I];
        points[I] = vertices[pointID];
        sd[I] = pointSignedDist[pointID];
    }

    points[points.size() - 1] = mesh.Cf()[faceID];
    sd[sd.size() - 1] = sigDistCalcPtr_->signedDistance(mesh.Cf()[faceID]);

    // TODO: crete indexed tris
    label face_centre_id = points.size() - 1;
    forAll(tris, I)
    {
        label pid_A = I;
        label pid_B = (I + 1)%tris.size();
        tris[I] = indexedTri{pid_A, pid_B, face_centre_id};
    }
    
    return std::make_tuple(tris, points, sd);
}


label surfaceMeshCellApproximation::nTets(const label cellID) const
{
    label nTet = 0;

    const auto& thisCell = this->mesh().cells()[cellID];
    const auto& faces = this->mesh().faces();

    for (const auto faceID : thisCell)
    {
        nTet += faces[faceID].nEdges();
    }

    return nTet;
}


void surfaceMeshCellApproximation::flagInterfaceCellFaces(surfaceScalarField& alpha) const
{
    const auto& cellFaces = this->mesh().cells();

    for (auto cellI : interfaceCellIDs_)
    {
        const auto& interfaceCell = cellFaces[cellI];

        for (auto faceI : interfaceCell)
        {
            alpha[faceI] = -1.0;
        }
    }
} 


label surfaceMeshCellApproximation::interfaceCellVolumeFraction(
        volScalarField& alpha, bool writeDecomposition)
{
    const auto& V = this->mesh().V();
    label maxRefine = 0;

    // TODO (TT): OpenMP disabled for now. Loop does not execute correct with
    // more than two threads. See issue on GitLab.
    //#pragma omp parallel for reduction(max:maxRefine)
    for (const auto cellID : interfaceCellIDs_)
    {
        auto [tets, points, signed_dist] = decomposeCell(cellID);

        adaptiveTetCellRefinement<signedDistanceCalculator,
            boundingBallCriterion>
            refiner{this->sigDistCalc(),
                points,
                signed_dist,
                tets,
                maxAllowedRefinementLevel_};
        tetVolumeFractionCalculator vofCalc{};
        alpha[cellID] =
            vofCalc.accumulatedOmegaPlusVolume(refiner.resultingTets(),
                refiner.signedDistance(),
                refiner.points()) /
            V[cellID];

        // Bound volume fraction field
        alpha[cellID] = max(min(alpha[cellID], 1.0), 0.0);

        maxRefine = std::max(refiner.refinementLevel(), maxRefine);

        if (writeDecomposition)
        {
            refiner.writeTets(cellID);
        }
    }

    return maxRefine;
}


label surfaceMeshCellApproximation::intersectedFacesAreaFraction(
    surfaceScalarField& alpha,
    bool writeDecomposition
)
{
    label maxRefine = 0;
    
    forAll(alpha, faceID)
    {
        if (alpha[faceID] != -1.0)
        {
            continue;
        }
            
        auto [tris, points, signed_dist] = decomposeFace(faceID);
        
        adaptiveTriFaceRefinement<signedDistanceCalculator,
            boundingBallCriterion>
            refiner{this->sigDistCalc(),
                points,
                signed_dist,
                tris,
                maxAllowedRefinementLevel_};
        // TODO: enable once implemented
        // tetVolumeFractionCalculator vofCalc{};
        // alpha[cellID] =
        //     vofCalc.accumulatedOmegaPlusVolume(refiner.resultingTets(),
        //         refiner.signedDistance(),
        //         refiner.points()) /
        //     V[cellID];

        // Bound area fraction field
        alpha[faceID] = max(min(alpha[faceID], 1.0), 0.0);

        maxRefine = std::max(refiner.refinementLevel(), maxRefine);

        if (writeDecomposition)
        {
            refiner.writeTris(faceID);
        }
    }
    
    return maxRefine;
}


void surfaceMeshCellApproximation::determineRefinementLevel(volScalarField& alpha)
{
    /* Currently, there are three ways how the maximum refinement level is set:
        1) User prescribes maximum refinement level directly.
            -> integer value > 0
            -> Nothing more to be done
        2) "Auto mode": compute the level based on tetrahedra and triangle
            edge length. Limited to triangulated surfaces.
            -> value < 0
            -> The refinement class handles the computation of the refinement
                level
        3) "Accuracy driven": compute the refinement level based on a
            user-prescribed target accuracy of the relative global volume
            error and an initial error without refinement. Limited to closed,
            implicit surfaces.
            -> relVolumeErrorThreshold > 0.0
            -> Compute refinement level in this function
    */ 

    if (relVolumeErrorThreshold_ > 0.0)
    {
        // Compute volume fraction without refinement. Only use tets from
        // tetrahedral cell decomposition
        maxAllowedRefinementLevel_ = 0;
        interfaceCellVolumeFraction(alpha, false);

        // TODO (TT): for parallel cases using domain decomposition: one processor
        // should do the calculation and broadcast the computed refinement level.
        scalar Valpha = gSum((this->mesh().V() * alpha)());
        scalar Vref = this->sigDistCalc().surfaceEnclosedVolume();

        auto relativeVolumeError = mag(Valpha - Vref)/Vref;
        auto errorRatio = relativeVolumeError/relVolumeErrorThreshold_;

        // Ideally, this factor should be 2^n, where n is the order of the
        // volume fraction approximation for a single tetrahedron; for the method
        // employed here n=2.
        // In actual numerical experiments for a spherical interface, the value
        // below was observed for consecutive refinement levels.
        const scalar beta = 3.9;

        maxAllowedRefinementLevel_ =
            std::ceil(std::log(errorRatio)/std::log(beta));

        Info<< "Initial relative volume error: " << relativeVolumeError
            << nl
            << "Target relative volume error: " << relVolumeErrorThreshold_
            << nl
            << "Refinement level required: " << maxAllowedRefinementLevel_
            << nl
            << endl;
        }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
surfaceMeshCellApproximation::surfaceMeshCellApproximation(
    const dictionary& configDict, const fvMesh& mesh)
    : volumeFractionCalculator{configDict, mesh},
      sigDistCalcPtr_{
          signedDistanceCalculator::New(configDict.subDict("distCalc"), mesh)},
      interfaceCellIDs_{}, maxAllowedRefinementLevel_{
                               configDict.get<label>("refinementLevel")},
      relVolumeErrorThreshold_{configDict.get<scalar>("relError")}
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void surfaceMeshCellApproximation::calcVolumeFraction(volScalarField& alpha)
{
    bulkVolumeFraction(alpha);
    findIntersectedCells();
    determineRefinementLevel(alpha);

    Info << "Computing volume fraction for interface cells..." << endl;
    Info << "Number of cells flagged as interface cells: "
         << interfaceCellIDs_.size() << endl;

    maxUsedRefinementLevel_ = interfaceCellVolumeFraction(alpha,
        this->writeGeometry());

    Info << "Finished volume fraction calculation" << nl << endl;
}


void surfaceMeshCellApproximation::calcAreaFraction(surfaceScalarField& alpha)
{
    bulkAreaFraction(alpha);    
    if (interfaceCellIDs_.size() == 0)
    {
        findIntersectedCells();
    }
    flagInterfaceCellFaces(alpha);
    // There is no need to determine the maximum allowed refinement level.
    // Either it has been explicitly prescribed by the user, has been
    // computed for the volume fractions before hand or the auto mode is
    // used.
    Info << "Computing area fractions for faces of interface cells..." << endl;
    Info << "Number of cells flagged as interface cells: "
         << interfaceCellIDs_.size() << endl;

    maxUsedRefinementLevel_ = intersectedFacesAreaFraction(alpha,
        this->writeGeometry());

    Info << "Finished volume fraction calculation" << nl << endl;
    
}


void surfaceMeshCellApproximation::findIntersectedCells()
{
    const auto& cellClosestPoint = sigDistCalcPtr_->cellClosestPoint();
    const auto& cellSignedDist = sigDistCalcPtr_->cellSignedDist0();
    const auto& centres = this->mesh().C();
    const auto& points = this->mesh().points();
    const auto& meshCellPoints = this->mesh().cellPoints();

    forAll(cellClosestPoint, cellI)
    {
        auto distSqr = pow(cellSignedDist[cellI], 2.0);

        if
        (
            cellClosestPoint[cellI].hit()
            &&
            considerIntersected(centres[cellI], distSqr, meshCellPoints[cellI],
                points, std::vector<scalar>{}, boundingBallCriterion{})
        )
        {
            interfaceCellIDs_.push_back(cellI);
        }
    }
}


void surfaceMeshCellApproximation::writeFields() const
{
    sigDistCalcPtr_->writeFields();

    // Write identified interface cells as field
    volScalarField interfaceCells{
        "interfaceCells", sigDistCalcPtr_->cellSignedDist()};
    interfaceCells = dimensionedScalar{"interfaceCells", dimLength, 0};

    for (const auto idx : interfaceCellIDs_)
    {
        interfaceCells[idx] = 1.0;
    }

    interfaceCells.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace Foam::TriSurfaceImmersion

// ************************************************************************* //
