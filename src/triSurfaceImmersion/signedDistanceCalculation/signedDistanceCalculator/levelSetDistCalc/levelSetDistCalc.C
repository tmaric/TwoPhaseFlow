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

#include "levelSetDistCalc.H"

#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvcAverage.H"
#include "surfaceInterpolate.H"
#include "volFieldsFwd.H"

namespace Foam::TriSurfaceImmersion
{

defineTypeNameAndDebug(levelSetDistCalc, 0);
addToRunTimeSelectionTable(
    signedDistanceCalculator, levelSetDistCalc, Dictionary);

// Implementation of 'surfacePoint' and 'closestPoint' taken from
// T. Maric's OpenFOAM machine learning project:
// https://gitlab.com/tmaric/openfoam-ml
point levelSetDistCalc::surfacePoint(const point& p) const
{
    vector x_j = p;
    scalar val_j{0.0};
    vector grad_j{0, 0, 0};

    // Descend onto the surface with gradient descent.
    for (label j = 0; j < maxIt_; ++j)
    {
        // Uncomment for debugging info.
        // Info << "K ITERATION                          : " << j << endl;

        val_j = surfacePtr_->value(x_j);
        grad_j = surfacePtr_->grad(x_j);

        // TODO (TT): divide by zero stablization.
        // This is required as coincidence of point p with the centre
        // cannot be ruled out.
        if (mag(grad_j) < SMALL)
        {
            grad_j = vector{SMALL, SMALL, SMALL};
        }

        x_j = x_j - val_j * grad_j / (grad_j & grad_j);

        // Uncomment for debugging info
        // Info << "x_j = " << x_j << endl;
        // Info << "val_j = " << surf.value(x_j) << endl;

        if (mag(val_j) < epsilon_)
        {
            break;
        }
    }

    return x_j;
}


point levelSetDistCalc::closestPoint(const point& p) const
{
    vector p_i = surfacePoint(p);

    // Modified closest-point algorithm from:
    //
    //   Hartmann, Erich. "Geometry and algorithms for computer aided design."
    //   Darmstadt University of Technology 95 (2003).
    //
    // instead of linear tangential projection that fails with zero
    // gradients on some surfaces, a parabolic projection is always
    // used.

    // Move the initial surface point in the tangential direction with respect
    // to x.
    const std::size_t MAX_I = 1;
    for (std::size_t i = 0; i < MAX_I; ++i)
    {
        vector n_i = surfacePtr_->grad(p_i);
        // TODO (TT): divide by zero stablization using SMALL.
        // This is required as coincidence of point p with the centre
        // cannot be ruled out.
        if (mag(n_i) < SMALL)
        {
            n_i = vector{SMALL, SMALL, SMALL};
        }
        else
        {
            n_i /= Foam::mag(n_i);
        }

        vector q_i = p - ((p - p_i) & n_i) * n_i;

        // First-order tangential projection.
        // vector p_i1 = surfacePoint(q_i, surf, maxIt, epsilon);

        // Second-order parabolic projection
        vector f_i = q_i - p_i;
        vector x_i = p_i + f_i - (f_i & f_i) * n_i;
        vector p_i1 = surfacePoint(x_i);

        scalar error = mag(p_i - p_i1);
        p_i = p_i1;

        if (error < epsilon_)
        {
            break;
        }
    }

    // Uncomment for debugging
    // Info << "ps = " << p_i << endl;
    // Info << "val_s = " << surf.value(p_i) << endl;

    return p_i;
}


void levelSetDistCalc::computeMeshLevelSetValues()
{
    const auto& cellCentres = this->mesh().C();

    forAll(cellCentres, I)
    {
        cellLevelSetValues_[I] = surfacePtr_->value(cellCentres[I]);
    }

    const auto& cellCorners = this->mesh().points();

    forAll(cellCorners, I)
    {
        pointLevelSetValues_[I] = surfacePtr_->value(cellCorners[I]);
    }
}


void levelSetDistCalc::identifyNarrowBandCells()
{
    // Active narrow band
    if (this->narrowBandWidth() > 0.0)
    {
        volScalarField isNarrowBandCell{"isNarrowBandCell", cellSignedDist0_};
        isNarrowBandCell *= 0.0;

        // Step 1: tag all cells intersected by the surface
        const auto& mesh = this->mesh();
        const auto& meshCellPoints = mesh.cellPoints();

        forAll(cellLevelSetValues_, cellI)
        {
            const auto& cellDist = cellLevelSetValues_[cellI];
            const auto& cellPoints = meshCellPoints[cellI];

            forAll(cellPoints, pointI)
            {
                if ((pointLevelSetValues_[cellPoints[pointI]] * cellDist) < 0)
                {
                    isNarrowBandCell[cellI] = 1.e15;
                    break;
                }
            }
        }

        // Step 2: use explicit pseudo diffusion to tag cells adjacent to
        // intersected cells as narrow band cells
        auto nbWidth = int(std::ceil(this->narrowBandWidth()));
        for (int i = 0; i != nbWidth; ++i)
        {
            isNarrowBandCell = fvc::average(fvc::interpolate(isNarrowBandCell));
        }

        // Step 3: collect tagged cell IDs
        forAll(isNarrowBandCell, I)
        {
            if (isNarrowBandCell[I] > 0.0)
            {
                narrowBandCells_.push_back(I);
            }
        }
    }
    else
    {
        // Narrow band disabled
        narrowBandCells_.resize(this->mesh().nCells());
        forAll(this->mesh().C(), I)
        {
            narrowBandCells_[I] = I;
        }
    }
}


void levelSetDistCalc::computeSignedDistances()
{
    const auto& cellCentres = this->mesh().C();
    const auto& cellCorners = this->mesh().points();
    const auto cellToPoint = this->mesh().cellPoints();
    dimensionedScalar defaultValue{
        "default", dimLength, this->outOfNarrowBandValue()};
    cellSignedDist0_ = defaultValue;
    pointSignedDist_ = defaultValue;

    for (const auto cellID : narrowBandCells_)
    {
        auto nearestPoint = closestPoint(cellCentres[cellID]);
        cellNearestTriangle_[cellID] = pointIndexHit{true, nearestPoint, 1};
        cellSignedDist0_[cellID] = sign(cellLevelSetValues_[cellID]) *
            mag(cellCentres[cellID] - nearestPoint);

        for (const auto pID : cellToPoint[cellID])
        {
            if (pointSignedDist_[pID] == this->outOfNarrowBandValue())
            {
                nearestPoint = closestPoint(cellCorners[pID]);
                pointNearestTriangle_[pID] =
                    pointIndexHit{true, nearestPoint, 1};
                pointSignedDist_[pID] = sign(pointLevelSetValues_[pID]) *
                    mag(cellCorners[pID] - nearestPoint);
            }
        }
    }
}


void levelSetDistCalc::setInsideOutside()
{
    // Cell centre inside/outside
    forAll(cellSignedDist_, I)
    {
        if (cellNearestTriangle_[I].hit())
        {
            cellSignedDist_[I] = cellSignedDist0_[I];
        }
        else
        {
            cellSignedDist_[I] = cellLevelSetValues_[I];
        }
    }

    // Cell corner inside/outside
    forAll(pointSignedDist_, I)
    {
        if (!pointNearestTriangle_[I].hit())
        {
            pointSignedDist_[I] = pointLevelSetValues_[I];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
levelSetDistCalc::levelSetDistCalc(
    const dictionary& configDict, const fvMesh& mesh)
    : signedDistanceCalculator{configDict, mesh},
      surfacePtr_{implicitSurface::New(configDict)},
      maxIt_{configDict.getOrDefault<label>("maxIter", 100)},
      epsilon_{configDict.getOrDefault<scalar>("epsilon", 1.e-15)},
      cellLevelSetValues_{IOobject("cellLevelSetValues",
                              mesh.time().timeName(),
                              mesh,
                              IOobject::NO_READ,
                              IOobject::NO_WRITE),
          mesh,
          dimensionedScalar("cellLevelSetValue", dimless, 0),
          "zeroGradient"},
      pointLevelSetValues_{IOobject("pointLevelSetValues",
                               mesh.time().timeName(),
                               mesh,
                               IOobject::NO_READ,
                               IOobject::NO_WRITE),
          this->pMesh(),
          dimensionedScalar("pointLevelSetValues", dimless, 0),
          "zeroGradient"},
      narrowBandCells_{}
{
    cellNearestTriangle_.resize(mesh.nCells());
    pointNearestTriangle_.resize(mesh.nPoints());

    computeMeshLevelSetValues();
    identifyNarrowBandCells();
    computeSignedDistances();
    setInsideOutside();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar levelSetDistCalc::signedDistance(const point& x) const
{
    return sign(surfacePtr_->value(x)) * mag(x - closestPoint(x));
}


scalar levelSetDistCalc::referenceLength() const
{
    // TODO (TT): implement a surface dependent reference length,
    // e.g. based on maximum curvature
    return 1.0;
}


label levelSetDistCalc::nSurfaceElements() const
{
    return 1;
}


scalar levelSetDistCalc::surfaceEnclosedVolume() const
{
    return surfacePtr_->volume();
}


void levelSetDistCalc::writeFields() const
{
    this->signedDistanceCalculator::writeFields();
    cellLevelSetValues_.write();
    pointLevelSetValues_.write();

    volScalarField narrowBand{"narrowBand", cellLevelSetValues_};
    narrowBand *= 0.0;

    for (const auto cellID : narrowBandCells_)
    {
        narrowBand[cellID] = 1.0;
    }

    narrowBand.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
