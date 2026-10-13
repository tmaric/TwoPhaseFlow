/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
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
    Foam::reconstruction::davof

Description
    DAVOF (dual area/volume-of-fluid) interface reconstruction: the PLIC
    normal of every interface cell is recovered from the face liquid-area
    fractions alpha_f of the cell's own faces by the Gauss identity, without
    iteration and without neighbouring cells,

        m_c = sum_f alpha_f S_f^out,   S_f^out = +S_f (owner), -S_f (neighbour),

    and the plane is positioned along m_c so that its cut volume equals
    alpha_c (surfaceIteratorPLIC::vofCutCell): normal_ is the cut polygon's
    area vector, centre_ its centroid.

    Sign conventions (TwoPhaseFlow's): normal_ and m_c point INTO the liquid
    (alpha = 1), as cutCell/cutFacePLIC assume ((p - base).n > 0 is liquid).
    leia's area normal is -m_c (out of the liquid).

    Interface cells: |m_c| > wispTol V_c^(2/3); m_c is set to zero exactly
    in cells whose face fractions are all equal (all wet or all dry: the sum
    of S_f vanishes only to round-off).

    The normal model (davofNormalModel, keyword normalModel; the models
    cellLocal, youngs, plicRDF, exact, ringArea, ringFit, each with applyTo
    all | check, sliverTol and maxAngle): after every plane has been
    positioned by alpha along m_c, the model returns the normal it proposes
    for every interface cell and the cells in which m_c is replaced (all of
    them, or the slivers |m_c| < sliverTol V_c^(2/3) and the cells deviating
    from the model normal by more than maxAngle degrees); the plane of every
    such cell is re-positioned by alpha along the model normal. The ring
    models need the positioned planes of the neighbours, hence the two
    passes. The deprecated keys sliverNormal and normalSource are mapped
    onto the models (davofNormalModel.H).

    Diagnostics (AUTO_WRITE), per cell:

        areaNormalDavof.<group>   m_c [m^2]
        sliverFlag.<group>        1 in sliver cells (0 otherwise)
        pPlane.<group>            the EXPLICIT plane position: the signed
                                  distance from the cell centre along
                                  n_c = m_c/|m_c| to the plane that the Gauss
                                  theorem for x on the liquid sub-volume
                                  implies (exact on planar faces),

            pPlane_c = ( sum_f alpha_f (x_f - x_c).S_f^out - 3 alpha_c V_c ) / |m_c|,

                                  the plane {x : (x - x_c).n_c = pPlane_c},
                                  liquid where (x - x_c).n_c > pPlane_c. This
                                  is cutCellPLIC's cutValue; it equals MINUS
                                  leia's p_c = (3 alpha V - sum_f ...)/|m|,
                                  whose n points out of the liquid.
        planeVolDefect.<group>    |alpha_plane - alpha_c|, alpha_plane the
                                  fraction cut by the explicit plane
        planeAreaDefect.<group>   |A_plane n_plane - m_c| / V_c^(2/3), the
                                  explicit plane's polygon area vector
                                  against m_c

    Both defects are round-off for a plane and O(h^2) on a curved surface
    whose state no single plane realizes (the realizability defect). They
    are always formed from m_c, whatever normal positions the plane.

    The face fractions alphaf.<group> are looked up from the registry; when
    not registered they are read from the time directory (initDavofFields
    writes them). Nothing updates them here: the face update is the
    transport step (isoAdvection stores the moving plane's submerged face
    area at the end of every step as alphaf, see isoAdvection.H), not the
    reconstruction.

    Dictionary (the reconstruction dictionary):
        wispTol      1e-3;      // leia's value (V^(2/3) scaled)
        surfCellTol  1e-8;
        isoFaceTol   1e-8;
        normalModel  cellLocal; // cellLocal | youngs | plicRDF | exact | ringArea | ringFit
        applyTo      check;     // all | check (the model's default otherwise)
        sliverTol    0.1;       // c: a sliver has |m_c| < c V_c^(2/3)
        maxAngle     15;        // degrees (applyTo check)
        exactNormal             // normalModel exact: as areaFractionError
        {
            functionType sphere;
            origin       (0.5 0.5 0.5);   // at t = 0
            radius       0.15;            // and the function's other keys
            U0           (1 0 0);
        }
        weight       area;      // ringFit: area | alpha

    The state time of the exact model: the first reconstruction of a time
    index is the advection step's and describes the state at t - dt, every
    later one in the same index (the function objects at the end of the
    step) the state at t, the constructor's the state at the start time;
    the time used is printed in the Info line.

SourceFiles
    davof.C

\*---------------------------------------------------------------------------*/

#ifndef davof_H
#define davof_H

#include "typeInfo.H"
#include "reconstructionSchemes.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "surfaceIteratorPLIC.H"
#include "cutCellPLIC.H"
#include "davofNormalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace reconstruction
{

/*---------------------------------------------------------------------------*\
                            Class davof Declaration
\*---------------------------------------------------------------------------*/

class davof
:
    public reconstructionSchemes
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Interface cell: |m_c| > wispTol V_c^(2/3)
        scalar wispTol_;

        //- Tolerance for marking of surface cells by alpha (vofCutCell)
        scalar surfCellTol_;

        //- Tolerance of the plane positioning by alpha
        scalar isoFaceTol_;

        //- Time index of the last reconstruction and the number of
        //  reconstructions in it (the state time of the exact model)
        label lastTimeIndex_;
        label nReconstructInIndex_;

        //- Finds the plane position for a given normal and alpha
        surfaceIteratorPLIC sIterPLIC_;

        //- Cuts the explicit plane for the defects
        cutCellPLIC cutCell_;

        //- The face fractions when this scheme had to read them
        autoPtr<surfaceScalarField> alphafPtr_;

        //- The face fractions alphaf.<group> (registered)
        surfaceScalarField& alphaf_;

        //- m_c = sum_f alpha_f S_f^out [m^2]
        volVectorField areaNormal_;

        //- 1 in sliver cells
        volScalarField sliverFlag_;

        //- The explicit plane position (cutValue convention, see above) [m]
        volScalarField pPlane_;

        //- |alpha_plane - alpha_c|
        volScalarField planeVolDefect_;

        //- |A_plane n_plane - m_c| / V_c^(2/3)
        volScalarField planeAreaDefect_;

        //- Maxima of the two defects over the interface cells (reduced)
        scalar maxPlaneVolDefect_;
        scalar maxPlaneAreaDefect_;

        //- The normal model (constructed after the fields, before the
        //  first reconstruction)
        autoPtr<davofNormalModel> normalModel_;


    // Private Member Functions

        //- Read alphaf.<group> from the time directory unless registered
        static autoPtr<surfaceScalarField> readAlphafIfAbsent
        (
            const volScalarField& alpha1
        );

        //- Position the plane of celli by alpha along n; sets normal_ and
        //  centre_ (zero when the cut fails). Returns the success.
        bool positionPlane(const label celli, const vector& n);

        //- The time the current state belongs to (see Description)
        scalar stateTime();

        //- No copy construct
        davof(const davof&) = delete;

        //- No copy assignment
        void operator=(const davof&) = delete;


public:

    //- Runtime type information
    TypeName("davof");


    //- Construct from components
    davof
    (
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U,
        const dictionary& dict
    );


    //- Destructor
    virtual ~davof() = default;


    // Member Functions

        //- Reconstruct interface
        virtual void reconstruct(bool forceUpdate = true);

        //- The mesh
        const fvMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- The volume fractions
        const volScalarField& alpha1() const noexcept
        {
            return alpha1_;
        }

        //- The volumetric flux
        const surfaceScalarField& phi() const noexcept
        {
            return phi_;
        }

        //- The velocity
        const volVectorField& U() const noexcept
        {
            return U_;
        }

        //- The face fractions
        const surfaceScalarField& alphaf() const
        {
            return alphaf_;
        }

        //- The face fractions
        surfaceScalarField& alphaf()
        {
            return alphaf_;
        }

        //- The area normal m_c
        const volVectorField& areaNormal() const
        {
            return areaNormal_;
        }

        //- The sliver flag
        const volScalarField& sliverFlag() const
        {
            return sliverFlag_;
        }

        //- The explicit plane position
        const volScalarField& pPlane() const
        {
            return pPlane_;
        }

        const volScalarField& planeVolDefect() const
        {
            return planeVolDefect_;
        }

        const volScalarField& planeAreaDefect() const
        {
            return planeAreaDefect_;
        }

        scalar maxPlaneVolDefect() const
        {
            return maxPlaneVolDefect_;
        }

        scalar maxPlaneAreaDefect() const
        {
            return maxPlaneAreaDefect_;
        }

        scalar wispTol() const
        {
            return wispTol_;
        }

        //- The normal model
        const davofNormalModel& normalModel() const
        {
            return *normalModel_;
        }

        scalar sliverTol() const
        {
            return normalModel_->sliverTol();
        }

        //- Sliver cells of the last reconstruction
        label nSliver() const
        {
            return normalModel_->nSliver();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reconstruction
} // End inline namespace twoPhaseFlow
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
