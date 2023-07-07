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

#include "searchDistanceCalculator.H"

#include "dimensionedScalarFwd.H"
#include "fvc.H"

namespace Foam::TriSurfaceImmersion
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
searchDistanceCalculator::searchDistanceCalculator(
    const fvMesh& mesh, scalar searchDistFactor)
    : mesh_{mesh},
      searchDistFactor_{searchDistFactor},
      cellSqrSearchDist_{IOobject("cellSqrSearchDist",
                             mesh_.time().timeName(),
                             mesh_,
                             IOobject::NO_READ,
                             IOobject::NO_WRITE),
          mesh,
          dimensionedScalar{
              "defaultCellSearchDist", pow(dimLength, 2), mesh.bounds().mag()}},
      pMesh_{mesh},
      cellsToPointsInterp_{mesh},
      pointSqrSearchDist_{IOobject("pointSqrSearchDist",
                              mesh_.time().timeName(),
                              mesh_,
                              IOobject::NO_READ,
                              IOobject::NO_WRITE),
          pMesh_,
          dimensionedScalar("pointSqrSearchDist", dimLength, 0),
          "zeroGradient"}
{
    // Narrow band around the interface based on search distances is only
    // enabled for for positive search distance factors. Negative values enable
    // signed distance calculation in the entrire domain. (TT
    if (searchDistFactor_ > 0.0)
    {
        cellSqrSearchDist_ = fvc::average(pow(mesh.deltaCoeffs(), -2)) *
            searchDistFactor_ * searchDistFactor_;
    }

    cellsToPointsInterp_.interpolate(cellSqrSearchDist_, pointSqrSearchDist_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void searchDistanceCalculator::writeFields() const
{
    cellSqrSearchDist_.write();
    pointSqrSearchDist_.write();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
