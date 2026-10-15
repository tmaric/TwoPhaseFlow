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
    Foam::davofFaceClipping

Description
    Run-time selectable model of the face update's stored face fraction
    (isoAdvection with the davof reconstruction scheme): what the face
    update stores on a face a plane has set, given the plane-clipped area
    (the submerged area of the face under the moving PLIC plane at the end
    of the step, cutFaceAdvect::subFaceAreaEnd()).

        plane         the plane-clipped area itself (today's face update:
                      bit-identical, no code path of the sweep changes;
                      the default)
        planeSagitta  the plane-clipped area plus the second-order
                      correction between the plane's trace on the face
                      (a chord) and the interface's trace (an arc), from
                      the curvature tensor the scheme fits to the ring
                      (planeSagittaClipping.H)

    Keyword faceClipping of the alpha solver dictionary (the dictionary
    isoAdvection and the davof scheme read), default plane. The model is
    constructed only when the reconstruction scheme is davof.

    isoAdvection calls correct() for every face a plane sets (the downwind
    faces of a surface cell, the upwind faces its plane sweeps, the
    boundary faces, the second ring; the processor sync exchanges the
    final values) with the plane's cell, the unit normal INTO the liquid,
    the plane's centroid translated to the end of the step, the chord
    (cutFaceAdvect::chordEnd(), the plane's trace on the face at the end
    of the step, two points; fewer when the face lies on one side) and
    the plane-clipped area, and stores the returned area as a fraction of
    the face area. The counters of the step are reported in the face
    update's Info line when the model corrects anything.

SourceFiles
    davofFaceClipping.C

\*---------------------------------------------------------------------------*/

#ifndef davofFaceClipping_H
#define davofFaceClipping_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "point.H"
#include "vector.H"

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
                      Class davofFaceClipping Declaration
\*---------------------------------------------------------------------------*/

class davofFaceClipping
{
protected:

    // Protected Data

        //- The davof scheme (the curvature tensor, the positioned planes)
        const reconstruction::davof& scheme_;

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Counters of the step (local; reduced in report()): faces
        //  corrected, corrections capped at the maximum fraction, stored
        //  values clamped to [0, |S_f|], faces skipped because the
        //  interface is nearly parallel to the face (cos), faces skipped
        //  because the chord is too short
        label nCorrected_;
        label nCapped_;
        label nClamped_;
        label nSkippedCos_;
        label nSkippedChord_;


    // Protected Member Functions

        //- No copy construct
        davofFaceClipping(const davofFaceClipping&) = delete;

        //- No copy assignment
        void operator=(const davofFaceClipping&) = delete;


public:

    //- Runtime type information
    TypeName("davofFaceClipping");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            davofFaceClipping,
            dictionary,
            (
                const dictionary& dict,
                const reconstruction::davof& scheme
            ),
            (dict, scheme)
        );


    // Constructors

        //- Construct from the alpha solver dictionary and the scheme
        davofFaceClipping
        (
            const word& type,
            const dictionary& dict,
            const reconstruction::davof& scheme
        );


    // Selectors

        //- Select from the keyword faceClipping of the alpha solver
        //  dictionary (default plane)
        static autoPtr<davofFaceClipping> New
        (
            const dictionary& dict,
            const reconstruction::davof& scheme
        );


    //- Destructor
    virtual ~davofFaceClipping() = default;


    // Member Functions

        //- Whether the model changes the plane-clipped area at all; false
        //  for plane, and isoAdvection then neither keeps the chords nor
        //  calls correct(): the sweep runs as it did without the model
        virtual bool corrects() const
        {
            return false;
        }

        //- Whether the model reads the scheme's curvature tensor, which
        //  isoAdvection then has the scheme fit after the reconstruction
        //  at the start of every step (the tensor of the start-of-step
        //  planes is used for the whole step)
        virtual bool needsCurvatureTensor() const
        {
            return false;
        }

        //- The area to store on face facei, set by the plane of cell
        //  planeCell with the unit normal nLiq INTO the liquid, whose
        //  centroid translated to the end of the step is xEnd; chord the
        //  end points of the plane's trace on the face at the end of the
        //  step (fewer than two: no trace); areaPlane the plane-clipped
        //  area at the end of the step
        virtual scalar correct
        (
            const label facei,
            const label planeCell,
            const vector& nLiq,
            const point& xEnd,
            const UList<point>& chord,
            const scalar areaPlane
        ) = 0;

        //- Reset the counters of a step
        void resetCounters();

        //- The counters of the step, reduced over the processors, as one
        //  clause of the face update's Info line
        virtual void report(Ostream& os) const;


    // Access

        label nCorrected() const noexcept
        {
            return nCorrected_;
        }

        label nCapped() const noexcept
        {
            return nCapped_;
        }

        label nClamped() const noexcept
        {
            return nClamped_;
        }

        label nSkippedCos() const noexcept
        {
            return nSkippedCos_;
        }

        label nSkippedChord() const noexcept
        {
            return nSkippedChord_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End inline namespace twoPhaseFlow
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
