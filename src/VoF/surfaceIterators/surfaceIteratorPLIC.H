/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 DLR
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
    Foam::surfaceIteratorPLIC

Description
    Finds the cutValue that matches the volume fraction

    Reference:
    \verbatim
        Henning Scheufler, Johan Roenby,
        Accurate and efficient surface reconstruction from volume
        fraction data on general meshes,
        Journal of Computational Physics, 2019,
        doi 10.1016/j.jcp.2019.01.009
    \endverbatim

Author
    Johan Roenby, DHI, all rights reserved.
    Modified Henning Scheufler, DLR

SourceFiles
    surfaceIteratorPLIC.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceIteratorPLIC_H
#define surfaceIteratorPLIC_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cutCellPLIC.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{

/*---------------------------------------------------------------------------*\
                      Class surfaceIteratorPLIC Declaration
\*---------------------------------------------------------------------------*/

class surfaceIteratorPLIC
{
    // Private Data

        //- Mesh whose cells and faces to cut at their intersection with an
        //- isosurface.
        const fvMesh& mesh_;

        //- Cuts the cell and returns volume centre and normal
        cutCellPLIC cutCell_;

        //- Tolerance for marking of surface cells:
        //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
        scalar surfCellTol_;

        //- The plane positioning (readPositioning; 2026-11-06): vofCutCell
        //  (the default, the code below: the cubic through four cuts and
        //  Newton, then the secant on the residual, to tol in alpha) or
        //  ccs (the Consecutive Cubic Spline positioning of
        //  davofCcsPositioning.H on the same cutCellPLIC cut: the cubic
        //  Hermite interpolant through the bracket ends with the exact
        //  slopes dV/dp = -A, A the cut polygon's area, started from the
        //  linear estimate between the vertex extrema, to positioningTol
        //  of V_c or the round-off floor, at most positioningMaxIter cuts;
        //  T. Maric, J. Comput. Phys.: X 11 (2021) 100093,
        //  doi:10.1016/j.jcpx.2021.100093, arXiv:2012.02227). The full and
        //  empty returns of vofCutCell (alpha1 within tol of 0 or 1) are
        //  the same on both paths.
        word positioning_;
        scalar positioningTol_;
        label positioningMaxIter_;

        //- ccs: the counters since the last reportPositioning (the cells
        //  positioned, the cuts, the largest cut count, the largest and
        //  the summed |V_cut - alpha V_c|/V_c of the returned cuts)
        label nPosCells_;
        label nPosEvals_;
        label maxPosEvals_;
        scalar maxPosRes_;
        scalar sumPosRes_;


    // Private Member Functions

        //- The ccs positioning of celli along the unit normal
        label ccsCutCell
        (
            const label celli,
            const scalar alpha1,
            const vector& normal
        );


public:

    // Constructors

        //- Construct from fvMesh and a scalarField
        //  Length of scalarField should equal number of mesh points
        surfaceIteratorPLIC(const fvMesh& mesh, const scalar tol);


    // Member Functions

        //- Determine if a cell is a surface cell
        bool isASurfaceCell(const scalar alpha1) const
        {
            return
            (
                surfCellTol_ < alpha1
             && alpha1 < 1 - surfCellTol_
            );
        }

        //- Select the positioning from the scheme's dictionary (keys
        //  positioning vofCutCell | ccs, positioningTol, positioningMaxIter;
        //  absent: vofCutCell, 1e-14, 30)
        void readPositioning(const dictionary& dict);

        //- The positioning in use
        const word& positioning() const noexcept
        {
            return positioning_;
        }

        //- With ccs: one Info line of the counters (reduced over the
        //  ranks), then the counters reset; nothing with vofCutCell
        void reportPositioning(Ostream& os);

        //- Finds matching cutValue for the given value fraction
        //  \return the cellStatus
        label vofCutCell
        (
            const label celli,
            const scalar alpha1,
            const scalar tol,
            const label maxIter,
            vector normal
        );

        //- The centre point of cutted volume
        const point& subCellCentre() const
        {
            return cutCell_.subCellCentre();
        }

        //- The volume of cutted volume
        scalar subCellVolume() const
        {
            return cutCell_.subCellVolume();
        }

        //- The centre of cutting isosurface
        const point& surfaceCentre() const
        {
            return cutCell_.faceCentre();
        }

        //- The area vector of cutting isosurface
        const vector& surfaceArea() const
        {
            return cutCell_.faceArea();
        }

        //- Volume of Fluid for cellI (subCellVolume_/mesh_.V()[cellI])
        scalar VolumeOfFluid() const
        {
            return cutCell_.VolumeOfFluid();
        }

        //- The cutValue
        scalar cutValue() const
        {
            return cutCell_.cutValue();
        }

        //- The cellStatus
        label cellStatus()
        {
            return cutCell_.cellStatus();
        }

        //- The points of the cutting isosurface in sorted order
        const DynamicList<point>& facePoints()
        {
            return cutCell_.facePoints();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End inline namespace twoPhaseFlow
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
