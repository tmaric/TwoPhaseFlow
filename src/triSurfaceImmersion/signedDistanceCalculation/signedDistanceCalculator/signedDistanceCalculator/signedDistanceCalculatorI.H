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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
const dictionary& signedDistanceCalculator::configDict() const
{
    return dict_;
}


const fvMesh& signedDistanceCalculator::mesh() const
{
    return mesh_;
}


const pointMesh& signedDistanceCalculator::pMesh() const
{
    return pMesh_;
}


scalar signedDistanceCalculator::narrowBandWidth() const
{
    return narrowBandWidth_;
}


const DynamicList<pointIndexHit>& signedDistanceCalculator::cellClosestPoint()
    const
{
    return cellNearestTriangle_;
}


const DynamicList<pointIndexHit>& signedDistanceCalculator::pointClosestPoint()
    const
{
    return pointNearestTriangle_;
}


const volScalarField& signedDistanceCalculator::cellSignedDist() const
{
    return cellSignedDist_;
}


const volScalarField& signedDistanceCalculator::cellSignedDist0() const
{
    return cellSignedDist0_;
}


const pointScalarField& signedDistanceCalculator::pointSignedDist() const
{
    return pointSignedDist_;
}


scalar signedDistanceCalculator::outOfNarrowBandValue() const
{
    return outOfNarrowBandValue_;
}

// ************************************************************************* //
