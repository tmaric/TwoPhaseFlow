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
    Foam::davofNormalModel

Description
    Run-time selectable source of the normal along which the davof
    reconstruction scheme positions the PLIC plane of an interface cell
    (reconstruction::davof). The scheme first positions every plane by
    alpha along the raw Gauss-identity normal m_c; the model then returns,
    for every interface cell, the normal it proposes (zero where it has
    none), and the base class decides where the raw normal is replaced:

        applyTo all      every interface cell with a model normal
        applyTo check    the sliver cells, |m_c| < sliverTol V_c^(2/3), and
                         the cells whose raw normal deviates from the model
                         normal by more than maxAngle degrees (the rule of
                         plicRDF's residual check)

    The scheme re-positions the plane along the model normal in every cell
    whose normal changed. A cell wanted by the rule whose model normal is
    zero keeps m_c and is counted ("kept m").

    Models (the keyword normalModel of the davof dictionary):

        cellLocal   the raw m_c/|m_c|, nothing replaced (the default)
        youngs      fvc::grad(alpha) (into alpha = 1 like m_c)
        plicRDF     the normal of a helper plicRDF reconstruction of the
                    same alpha (owned by the model; its fields stay
                    unregistered)
        exact       the exact normal of the implicit function of the
                    exactNormal sub-dictionary translated by U0 t (a control)
        ringArea    sum_k m_k / |sum_k m_k| over the point-neighbour ring of
                    the cell (the cell included; interface cells only): the
                    DAVOF normal of the union of the ring's cells, exact for
                    a plane on any mesh
        ringFit     the normal at the cell's own plane centroid of a
                    weighted least-squares quadric (Monge patch in the frame
                    of m_c) through the ring's plane centroids and unit
                    normals; ringArea where the fit is not possible

    Dictionary (the davof scheme's dictionary):
        normalModel cellLocal;   // cellLocal | youngs | plicRDF | exact | ringArea | ringFit
        applyTo     check;       // all | check (default: check for youngs
                                 // and plicRDF, all for the others)
        sliverTol   0.1;         // c: a sliver has |m_c| < c V_c^(2/3)
        maxAngle    15;          // degrees, the check rule
        exactNormal { ... }      // exact: as areaFractionError
        weight      area;        // ringFit: area | alpha

    Deprecated keys sliverNormal (none | Youngs | plicRDF) and normalSource
    (davof | exact | plicRDF) are mapped onto the models when normalModel
    is absent: sliverNormal Youngs -> youngs/check, sliverNormal plicRDF ->
    plicRDF/check (both with maxAngle 180, the slivers-only rule they had),
    normalSource exact -> exact/all, normalSource plicRDF -> plicRDF/all.

SourceFiles
    davofNormalModel.C

\*---------------------------------------------------------------------------*/

#ifndef davofNormalModel_H
#define davofNormalModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{

namespace reconstruction
{
    class davof;
}

/*---------------------------------------------------------------------------*\
                       Class davofNormalModel Declaration
\*---------------------------------------------------------------------------*/

class davofNormalModel
{
protected:

    // Protected Data

        //- The davof scheme (its fields: areaNormal, interface cells,
        //  alphaf, the positioned planes normal()/centre())
        const reconstruction::davof& scheme_;

        //- Reference to mesh
        const fvMesh& mesh_;

        //- all | check
        word applyTo_;

        //- Sliver: |m_c| < sliverTol_ V_c^(2/3)
        scalar sliverTol_;

        //- The check rule's angle [deg] and its cosine
        scalar maxAngle_;
        scalar cosMaxAngle_;

        //- Counts of the last select() (reduced): interface cells,
        //  slivers, normals replaced, cells wanted whose model normal was
        //  zero (kept m_c), fit fallbacks (ringFit)
        label nInterface_;
        label nSliver_;
        label nReplaced_;
        label nKept_;
        label nFallback_;


    // Protected Member Functions

        //- No copy construct
        davofNormalModel(const davofNormalModel&) = delete;

        //- No copy assignment
        void operator=(const davofNormalModel&) = delete;


public:

    //- Runtime type information
    TypeName("davofNormalModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            davofNormalModel,
            dictionary,
            (
                const dictionary& dict,
                const reconstruction::davof& scheme
            ),
            (dict, scheme)
        );


    // Constructors

        //- Construct from the davof dictionary and the scheme; the model's
        //  default applyTo
        davofNormalModel
        (
            const word& type,
            const dictionary& dict,
            const reconstruction::davof& scheme,
            const word& defaultApplyTo
        );


    // Selectors

        //- Select from the keyword normalModel of the davof dictionary
        //  (default cellLocal), mapping the deprecated sliverNormal and
        //  normalSource keys
        static autoPtr<davofNormalModel> New
        (
            const dictionary& dict,
            const reconstruction::davof& scheme
        );


    //- Destructor
    virtual ~davofNormalModel() = default;


    // Member Functions

        //- The model normal of every interface cell, indexed as the
        //  scheme's interfaceLabels (zero where the model has none);
        //  tState the time the state belongs to
        virtual void modelNormals(const scalar tState, vectorField& n) = 0;

        //- Whether the model replaces any normal (cellLocal: no)
        virtual bool replaces() const
        {
            return true;
        }

        //- Whether a cell with |m_c| = mMag is a sliver
        bool isSliver(const label celli, const scalar mMag) const;

        //- For every interface cell (indexed as interfaceLabels) the normal
        //  the scheme positions its plane along (nUse) and whether it
        //  differs from m_c (changed); m the raw area normal per cell
        void select
        (
            const scalar tState,
            const vectorField& m,
            vectorField& nUse,
            boolList& changed
        );

        //- One line of counts of the last select()
        void report(Ostream& os) const;


    // Access

        const word& applyTo() const noexcept
        {
            return applyTo_;
        }

        scalar sliverTol() const noexcept
        {
            return sliverTol_;
        }

        scalar maxAngle() const noexcept
        {
            return maxAngle_;
        }

        label nInterface() const noexcept
        {
            return nInterface_;
        }

        label nSliver() const noexcept
        {
            return nSliver_;
        }

        label nReplaced() const noexcept
        {
            return nReplaced_;
        }

        label nKept() const noexcept
        {
            return nKept_;
        }

        label nFallback() const noexcept
        {
            return nFallback_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End inline namespace twoPhaseFlow
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
