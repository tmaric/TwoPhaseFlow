/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Modified code Copyright (C) 2016-2019 OpenCFD Ltd.
    Modified code Copyright (C) 2018 Johan Roenby
    Modified code Copyright (C) 2019 DLR
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
    Foam::isoAdvection

Description
    Calculates the new VOF (alpha) field after time step dt given the initial
    VOF field and a velocity field U and face fluxes phi. The fluid transport
    calculation is based on an idea of using isosurfaces to estimate the
    internal distribution of fluid in cells and advecting such isosurfaces
    across the mesh faces with the velocity field interpolated to the
    isosurfaces.

    Reference:
    \verbatim
        Roenby, J., Bredmose, H. and Jasak, H. (2016).
        A computational method for sharp interface advection
        Royal Society Open Science, 3
        doi 10.1098/rsos.160405
    \endverbatim

    Original code supplied by Johan Roenby, DHI (2016)
    Modified Henning Scheufler, DLR

SourceFiles
    isoAdvection.C
    isoAdvectionTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef isoAdvection_H
#define isoAdvection_H

#include "advectionSchemes.H"
#include "fvMesh.H"
#include "volFieldsFwd.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "className.H"
#include "reconstructionSchemes.H"
#include "cutFaceAdvect.H"
#include "davofFaceClipping.H"
#include "davofCurvedGeometry.H"
#include "davofVolumeRule.H"
#include "interpolationCellPoint.H"
#include "pointFields.H"
#include "bitSet.H"
#include "zeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace reconstruction
{
    class davof;
}
}

namespace advection
{


/*---------------------------------------------------------------------------*\
                        Class isoAdvection Declaration
\*---------------------------------------------------------------------------*/

class isoAdvection
    : public advectionSchemes
{
    // Private data types

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;


    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Often used reference to alpha1 internal field
        scalarField& alpha1In_;

        //- Face volumetric water transport
        surfaceScalarField dVf_;

        //- Time spent performing interface advection
        scalar advectionTime_;

        //- store timeIndex
        label timeIndex_;

        // Switches and tolerances. Tolerances need to go into toleranceSwitches

            //- Number of alpha bounding steps
            label nAlphaBounds_;

            //- Tolerance for search of isoFace giving specified VOF value
            scalar isoFaceTol_;

            //- Tolerance for marking of surface cells:
            //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
            scalar surfCellTol_;

            //- Print isofaces in a <case>/isoFaces/isoFaces_#N.vtk files.
            //  Intended for debugging
            bool writeIsoFacesToFile_;

        // Cell and face cutting

            //- List of surface cells
            DynamicLabelList surfCells_;

            //- An isoCutFace object to get access to its face cutting functionality
            cutFaceAdvect advectFace_;

            //- Storage for boundary faces downwind to a surface cell
            DynamicLabelList bsFaces_;

            //- Storage for boundary surface iso face centre
            DynamicVectorList bsx0_;

            //- Storage for boundary surface iso face normal
            DynamicVectorList bsn0_;

            //- Storage for boundary surface iso face speed
            DynamicScalarList bsUn0_;

            //- Storage for the velocity at the boundary surface iso face
            //  centre (the face-clipping models: the centroid at the end
            //  of the step); filled when the face fractions are stored
            DynamicVectorList bsU0_;



        // Additional data for parallel runs

            //- List of processor patch labels
            DynamicLabelList procPatchLabels_;

            //- For each patch if it is a processor patch this is a list of the
            //  face labels on this patch that are downwind to a surface cell.
            //  For non-processor patches the list will be empty.
            List<DynamicLabelList> surfaceCellFacesOnProcPatches_;


        // DAVOF face fractions (davof reconstruction): the moving plane's
        // submerged face area at the end of the step stored as alphaf.<group>

            //- Store the end-of-step face fractions (default: the
            //  reconstruction scheme is davof; dictionary key
            //  storeAreaFractions)
            bool storeAreaFractions_;

            //- The registered face fractions alphaf.<group>; nullptr when
            //  none is registered and storeAreaFractions_ is false
            surfaceScalarField* alphafPtr_;

            //- The face fractions created here when storeAreaFractions_
            //  and nothing was registered (initialised upwind)
            autoPtr<surfaceScalarField> alphafOwnPtr_;

            //- Faces set by a plane sweep in this step
            bitSet alphafSetFace_;

            //- The faces alphafSetFace_ marks, as a list (set by a plane
            //  or received in the sync), of this step and of the previous
            //  one; the faces whose sweep fell back to the default
            DynamicLabelList alphafSetFaces_;
            DynamicLabelList alphafPrevSetFaces_;
            DynamicLabelList alphafFallbackFaces_;

            //- The alpha the default face fractions were last evaluated
            //  with (a copy of alpha at the end of the previous step):
            //  every face no plane set holds alphaDefault_ of its upwind
            //  cell (boundary faces: its owner) as long as defaultsValid_
            bool defaultsValid_;
            scalarField alphaDefault_;

            //- The upwind direction of every internal face at the last
            //  evaluation (set: phi >= 0, the owner) and phi's event
            //  number then (phi written since: the directions are
            //  re-scanned)
            bitSet phiPos_;
            label phiEventNo_;

            //- The faces whose fraction is not 0 (kept at every write of
            //  a face value), and their list in ascending face order at
            //  the end of the step, handed to the davof scheme, which
            //  sums the area normal over them in that order
            bitSet wetFace_;
            labelList wetFaces_;

            //- The davof scheme when the reconstruction is one (nullptr
            //  otherwise)
            reconstruction::davof* davofPtr_;

            //- Counters of the step: faces set by their own cell's plane,
            //  faces set by a neighbouring cell's plane (second ring),
            //  sweeps that could not set the face (default kept), values
            //  clipped to [0, 1] by more than 1e-12
            label nAlphafFromPlane_;
            label nAlphafSecondRing_;
            label nAlphafFallback_;
            label nAlphafClipped_;

            //- Faces beyond the downwind faces of the surface cells whose
            //  flux dVf is the plane's (planeFlux)
            label nFluxSecondRing_;

            //- Extend the conservative bounding (limitFluxes) by a second
            //  layer of cells when the face fractions are stored: the
            //  second-ring sweep sets plane fluxes on faces whose other
            //  cell lies beyond the first layer (dictionary key
            //  boundSecondLayer, default true)
            bool boundSecondLayer_;

            //- Cells beyond the first layer the extended bounding corrected
            //  in the step (summed over the bounding passes)
            label nBoundedBeyondRing_;

            //- The per-cell content bound of the plane-set fluxes
            //  (boundCellContent; dictionary key contentBound, default true)
            bool contentBound_;

            //- Per face: the non-surface upwind cell for which planeFlux
            //  set the plane's flux in this step (-1 otherwise), and the
            //  upwind flux it replaced
            labelList planeFluxUpCell_;
            scalarField dVfUpwindPlaneFace_;

            //- Faces and cells the content bound moved in the step
            label nContentBoundFaces_;
            label nContentBoundCells_;

            //- The cell whose content bound is printed when debug is on
            //  (dictionary key debugCell, default -1)
            label debugCell_;

            //- The plane of every surface cell of the step (x0, n0, Un0)
            //  and the velocity U0 at x0 (Un0 = U0.n0), indexed as
            //  surfCells_; cellToSurf_ maps a cell to that index (-1: not
            //  a surface cell)
            DynamicVectorList surfX0_;
            DynamicVectorList surfN0_;
            DynamicScalarList surfUn0_;
            DynamicVectorList surfU0_;
            labelList cellToSurf_;

            //- The face-clipping model (davofFaceClipping, keyword
            //  faceClipping, default plane): what is stored on a face a
            //  plane sets, given the plane-clipped area; constructed when
            //  the reconstruction scheme is davof. clipCorrects_ is its
            //  corrects(): false for plane, and then nothing of the
            //  sweep changes (no chord is kept, correct() is not called)
            autoPtr<davofFaceClipping> faceClippingPtr_;
            bool clipCorrects_;

            //- The faces planeFlux set in the step (planeFluxUpCell_ is
            //  >= 0 exactly there); planeFluxUpCell_, dVfUpwindPlaneFace_
            //  and cellToSurf_ are sized once and reset through this list
            //  and surfCells_ at the end of the step, not over the mesh
            DynamicLabelList planeFluxFaces_;

            //- The second-ring sweep's per-face closest plane (an index
            //  into surfCells_) and its distance, sized once and reset at
            //  the candidate faces after the sweep; the candidate faces
            DynamicLabelList secondRingFaces_;
            labelList secondRingBest_;
            scalarField secondRingDist_;


        // DAVOF curved clipping surface (isoAdvectionCurved.C, 2026-10-23):
        // the dual advection with ONE surface per cell for the volume flux
        // and the face update, replacing the plane sweep when active

            //- none (the plane sweep) | exactSphere | paraboloid (key
            //  clippingSurface)
            word clippingSurface_;

            //- The face state on warped faces (key faceState, 2026-10-25):
            //  scalar (alphaf alone, the default) or vector (alphaf and
            //  the wetted area vector wettedAreaDavof.<group>, the sum over
            //  the fan triangles of the wetted area times the triangle
            //  normal, which carries the normal identity on non-planar
            //  faces; registered by the davof scheme, set here by the
            //  curved face update, the defaults alpha_up S_f and the
            //  processor sync with the sign flipped). Needs a clipping
            //  surface (the plane sweep produces no vector)
            word faceState_;
            surfaceVectorField* wetVecPtr_;

            //- The centroids the ring Monge fit of the curvature tensor
            //  reads (key tensorCentroid): cut (the alpha-positioned
            //  plane's cut-polygon centroid, the default) or anchored (the
            //  plane in the frame of the used normal re-anchored by the
            //  volume rule's cell volume, o + c e3: on a warped cell the
            //  cut polygon closes the wet part of a face by a chord, so
            //  the cut centroid is off the fan polyhedron's plane by
            //  O(warp) and the fitted tensor of an exact plane is not
            //  zero; counted as anchor-centroid failures where the root
            //  is not bracketed)
            word tensorCentroid_;
            label nCentroidFail_;

            //- The prism quadrature: relative tolerance of the flux volume
            //  and the maximum refinement depth (keys curvedTolRel,
            //  curvedMaxDepth)
            scalar curvedTolRel_;
            label curvedMaxDepth_;

            //- The six-coefficient fit: the prior weight of the ring
            //  tensor, the face-residual tolerance, the iteration cap
            //  (keys fitPriorWeight, fitTol, fitMaxIter)
            scalar fitPriorWeight_;
            scalar fitTol_;
            label fitMaxIter_;

            //- faces (the six coefficients from the face fractions, the
            //  surface interpolates them) or volume (the offset shifted
            //  after the fit so that the volume under the surface is
            //  alpha_k V, the feedback of the volume fraction into the
            //  surface; key fitAnchor)
            word fitAnchor_;

            //- six (all six coefficients from the face fractions, the
            //  Gauss-Newton), one (the offset alone from the total wetted
            //  area, the ring tensor and zero slopes: leia's paraboloidFaces)
            //  or none (the frame of the normal and the ring tensor only;
            //  the offset from the volume anchor); key fitDof
            word fitDof_;

            //- The volume anchor's tolerance (of the cell volume; key
            //  anchorTolRel) and the period of the surface-defect
            //  diagnostic in steps (key defectEvery)
            scalar anchorTolRel_;
            label defectEvery_;

            //- The quadrature diagnostic (key curvedDebug): the cell with
            //  the largest anchor residual before the anchor, printed with
            //  its pyramids against the linear tet cut
            bool curvedDebug_;

            //- The volume anchor's face-fraction diagnostics (key
            //  anchorDiagnostics, default off): the face residual of the
            //  surface before the anchor and the face shift after it, two
            //  face-fraction evaluations of every face of every interface
            //  cell per step (14 % of the advection step at N = 64,
            //  2026-10-25), printed in the Info line when on
            bool anchorDiagnostics_;

            //- The volume anchor's root finder (key anchorMethod): ccs (the
            //  default: Consecutive Cubic Spline positioning on the
            //  paraboloid, the root of the cubic Hermite interpolant through
            //  the bracket ends with the exact values and slopes dV/dc, the
            //  projected wetted patch area of the volume rule, re-bracketed
            //  and bisection-safeguarded; T. Maric, J. Comput. Phys.: X 11
            //  (2021) 100093, doi:10.1016/j.jcpx.2021.100093,
            //  arXiv:2012.02227), newton (safeguarded Newton on that slope),
            //  illinois (a bracket by doubling steps and the Illinois rule,
            //  the method until 2026-10-25) or compare (illinois and ccs,
            //  the run continues with ccs, the cells whose offsets differ
            //  counted and printed). ccs and newton take the edge of the
            //  flat root interval for a cell whose target is 0 or V (see
            //  isoAdvectionCurved.C); illinois keeps its bracket end there.
            word anchorMethod_;

            //- The volume flux of a face (key fluxIntegration): prism (the
            //  default, the liquid volume of the prism between the face's
            //  pre-image and the face), trapezoidal (phi_f dt times the
            //  mean of the face fractions at t^n and t^{n+1}) or simpson
            //  (with the half-step pre-image's fraction, phi_f dt
            //  (a^n + 4 a^{n+1/2} + a^{n+1})/6) or simpsonSplit (Simpson's
            //  rule on the sub-intervals between the passages of the trace
            //  over the face's vertices), simpsonSplitS (the same with every
            //  value under the surface S: the face's own fraction at
            //  sigma = 0 instead of the stored one, one-sided limits at the
            //  breakpoints) or gaussSplit (the three-point Gauss rule on
            //  those sub-intervals, no value at a breakpoint and none from
            //  the stored fraction), or sweptFace (the semi-Lagrangian
            //  swept-face flux: the face's vertices and centre traced
            //  backward, the swept volume's triangles clipped by the
            //  donor's PLANE, the flux and the new face fraction from the
            //  same clips; planes only, one donor cell; isoAdvectionSwept.C,
            //  2026-10-28); see isoAdvectionCurved.C
            word fluxIntegration_;

            //- fluxIntegration sweptFace: the tolerance of the plane test
            //  max|b_ij| h of the donor's surface (key sweptPlaneTol,
            //  default 1e-6; a surface above it is a FatalError) and the
            //  counters of the step: the faces computed, the sign flips of
            //  the triangle list, the faces whose pre-image leaves the
            //  donor, the clipped triangles, the largest closure residual
            scalar sweptPlaneTol_;
            label nSweptFaces_;
            label nSweptFlips_;
            label nSweptOutOfDonor_;
            label nSweptClipTris_;
            scalar sweptClosureMax_;

            //- The pre-image classification diagnostic (key sweptClassify,
            //  default false; isoAdvectionSwept.C, 2026-10-29), for the
            //  prism and the swept flux alike: per face the class of the
            //  cells its pre-image points (vertices and centre) reach
            //  beyond the donor, 0 the donor only, 1 an interface cell (a
            //  cell with a surface), 2 a full cell, 3 an empty cell, 4 a
            //  cell outside the donor's vertex ring or outside the mesh
            //  (the priority 1 > 2 > 3 where points reach several), counted
            //  per step; the first cell that leaves [0, 1] by more than
            //  1e-10 after the conservative bounding is printed once with
            //  its faces' fluxes before and after the bounding and their
            //  classes (serial runs); the mass sum(alpha V) after the
            //  advection, the bounding and the clipping is printed every
            //  step relative to its value at t^n
            bool sweptClassify_;
            labelList sweptClass_;
            FixedList<label, 5> nSweptClass_;
            bool sweptEventReported_;
            scalarField dVfBeforeBound_;

            //- fluxIntegration sweptFaceCells (the piecewise clip,
            //  2026-10-29): the counters of the step: the pieces with area,
            //  the faces whose boundary is not covered by the donor's vertex
            //  ring (a part beyond one layer, skipped) and the largest
            //  uncovered share, the faces that fell back to the one-donor
            //  clip on a rank-deficient cap system, the partial cells
            //  without a surface met, the largest closure residual of the
            //  pieces
            label nSweptPieces_;
            label nSweptBeyond_;
            scalar sweptBeyondMax_;
            label nSweptRankFallback_;
            label nSweptPartialCells_;
            scalar sweptPieceClosureMax_;

            //- The kinematic consistency of the swept fluxes (key
            //  sweptScaleToPhi, default false, 2026-10-29): the liquid volume
            //  of the swept volume scaled by |phi_f| dt / V_swept, so that a
            //  wholly wet swept volume gives the flux phi_f dt that the
            //  upwind flux of a full cell uses. Under a rotation the traced
            //  swept volume (straight edges between traced points) differs
            //  from phi_f dt at the order of the angular displacement per
            //  step squared, and a full cell next to an interface cell
            //  received 3e-7 more than it sent in the first step of the
            //  N = 64 rotation (the prism's volume d.S_f IS phi_f dt). The
            //  identity to round-off under a uniform translation.
            bool sweptScaleToPhi_;
            scalar debugR0Max_;
            label debugCellR0_;
            davof::cellSurface debugSurface_;
            scalar debugAlpha_;
            label nAnchorFail_;
            label nAnchorLeaves_;

            //- The cell-volume evaluations of the anchor in the step
            label nAnchorEvals_;

            //- The cells anchored at the edge of the flat interval (a
            //  target of 0 or V) in the step
            label nAnchorEdge_;

            //- anchorMethod compare: the cells printed in the step, the
            //  cells whose offsets differ, those of them at the edge
            label nCompare_;
            label nCompareDiff_;
            label nCompareDiffEdge_;
            label nCompareDiffResolved_;
            scalar anchorShiftMax_;
            scalar anchorShiftSum_;

            //- The refinement levels of the cell volume under the surface
            //  for the surface defect (key defectLevels; -1 switches the
            //  diagnostic off)
            label defectLevels_;

            //- Write the per-cell fit residual and surface defect fields
            //  (key writeCurvedFields)
            bool writeCurvedFields_;

            //- The surface of every cell that has one (interface cells)
            List<davof::cellSurface> cellSurfaces_;
            bitSet hasSurface_;

            //- The run-time selected rule for the cell volume, the prism
            //  volume and the face fraction under a surface (key
            //  volumeRule: splitQuadrature | kromerBothe | exactSphere)
            autoPtr<davofVolumeRule> volumeRulePtr_;
            label nRuleFallback_;

            //- Counters of the step
            label nCurvedFaces_;
            label nCurvedSecondRing_;
            label nCurvedBoundary_;
            label nCurvedSkipped_;
            label nCurvedFallback_;
            label nCurvedLeaves_;
            label nFitCells_;
            label nFitConverged_;
            label nFitStalled_;
            label nFitSingular_;
            label nFitNoBracket_;
            label nFitNoFrame_;
            scalar fitResMax_;
            scalar fitResSum_;
            scalar defectMax_;
            scalar defectSum_;
            label nDefectCells_;
            autoPtr<volScalarField> fitResidualPtr_;
            autoPtr<volScalarField> paraboloidDefectPtr_;


    // Private Member Functions

        //- No copy construct
        isoAdvection(const isoAdvection&) = delete;

        //- No copy assignment
        void operator=(const isoAdvection&) = delete;


        // Advection functions

            //- Extend markedCell with cell-face-cell.
            void extendMarkedCells(bitSet& markedCell) const;

            void setProcessorPatches();

            //- For each face calculate volumetric face transport during dt
            void timeIntegratedFlux();

            //- For a given cell return labels of faces fluxing out of this cell
            //  (based on sign of phi)
            void setDownwindFaces
            (
                const label celli,
                DynamicLabelList& downwindFaces
            ) const;

            // LimitFluxes
            template < class SpType, class SuType >
            void limitFluxes
            (
                const SpType& Sp,
                const SuType& Su
            );

            // Bound fluxes
            template < class SpType, class SuType >
            void boundFlux
            (
                const bitSet& nextToInterface,
                surfaceScalarField& dVfcorrectionValues,
                DynamicLabelList& correctedFaces,
                const SpType& Sp,
                const SuType& Su
            );

            //- Given the face volume transport dVf calculates the total volume
            //  leaving a given cell. Note: cannot use dVf member because
            //  netFlux is called also for corrected dVf
            scalar netFlux
            (
                const surfaceScalarField& dVf,
                const label celli
            ) const;

            //- Determine if a cell is a surface cell
            bool isASurfaceCell(const label celli) const
            {
                return
                (
                    surfCellTol_ < alpha1In_[celli]
                 && alpha1In_[celli] < 1 - surfCellTol_
                );
            }

            //- Clear out isoFace data
            void clearIsoFaceData()
            {
                surfCells_.clear();
                bsFaces_.clear();
                bsx0_.clear();
                bsn0_.clear();
                bsUn0_.clear();
                bsU0_.clear();
            }

        // Face value functions needed for random face access where the face
        // can be either internal or boundary face

            //- Return face value for a given Geometric surface field
            template<typename Type>
            Type faceValue
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
                const label facei
            ) const;

            //- Set face value for a given Geometric surface field
            template<typename Type>
            void setFaceValue
            (
                GeometricField<Type, fvsPatchField, surfaceMesh>& f,
                const label facei,
                const Type& value
            ) const;


        // Parallel run handling functions

            //- Synchronize dVf across processor boundaries using upwind value;
            //  alphaf (optional) is exchanged on the same faces, the upwind
            //  side's face fraction copied without the sign change
            DynamicList<label> syncProcPatches
            (
                surfaceScalarField& dVf,
                const surfaceScalarField& phi,
                bool returnSyncedFaces=false,
                surfaceScalarField* alphaf = nullptr,
                surfaceVectorField* wetVec = nullptr
            );


        // DAVOF face fractions

            //- Look up alphaf.<group> in the registry; create it (upwind)
            //  when storeAreaFractions_ and none is registered
            void initAreaFractions();

            //- The default of one face from the alpha field given: the
            //  upwind cell's value (phi >= 0: the owner), boundary faces
            //  the owner's; faces of empty patches carry no value
            void setDefaultAreaFraction
            (
                const label facei,
                const scalarField& alpha
            );

            //- After the bounding: every face not set by a plane in this
            //  step holds the upwind cell's new alpha. Evaluated on the
            //  faces whose stored value can differ from it only (the
            //  faces set in the previous step, the fallbacks of this one,
            //  the faces of the cells whose alpha changed bitwise against
            //  alphaDefault_, the faces whose upwind direction changed);
            //  on every face when defaultsValid_ is false
            void defaultAreaFractions();

            //- The list of the wet faces from wetFace_ (bitSet::toc(),
            //  ascending), handed to the davof scheme
            void publishWetFaces();

            //- Store the end area of the last sweep on facei as the face
            //  fraction (clipped to [0, 1]); a negative end area writes
            //  the default (the upwind cell's alpha at the start of the
            //  step) and counts as a fallback. The plane that swept the
            //  face (its cell planeCell, x0, n0, U0 at x0) and dt are
            //  handed to the face-clipping model when it corrects the
            //  plane-clipped area (clipCorrects_); otherwise unread
            void storeAreaFraction
            (
                const label facei,
                const scalar areaEnd,
                const bool secondRing,
                const label planeCell,
                const point& x0,
                const vector& n0,
                const vector& U0,
                const scalar dt
            );

            //- The faces not set by any plane whose cell shares a face with
            //  a surface cell: swept by the neighbouring plane whose
            //  end-of-step position is closest to the face
            void sweepSecondRing(const scalar dt);

            //- After a plane sweep of a face whose upwind cell is not a
            //  surface cell: does the plane speak for the face? Yes for a
            //  partial (sliver) upwind cell, and for a pure one when the
            //  plane crosses the face during the step (faceSwept)
            bool planeSpeaks(const label upCell) const;

            //- The flux of such a face: the plane's where it speaks, else
            //  the upwind flux; records the face as plane-set for upCell
            scalar planeFlux
            (
                const label facei,
                const label upCell,
                const scalar dVfPlane,
                const scalar dVfUpwind
            );

            //- The content bound: for every non-surface cell W with faces
            //  whose flux a plane set (planeFlux), the liquid leaving W
            //  through its faces may not exceed alpha_W V_W plus the liquid
            //  entering it in this step, nor the gas (1 - alpha_W) V_W plus
            //  the entering gas; when violated, the plane-set fluxes of W
            //  are moved towards their upwind values by one factor
            //  lambda in [0, 1] chosen so that the bound holds with
            //  equality (the upwind fluxes satisfy it for cell Courant
            //  numbers below one). Counts the bounded faces and cells.
            void boundCellContent(const scalar dt);

            //- Check if the face is on processor patch and append it to the
            //  list of surface cell faces on processor patches
            void checkIfOnProcPatch(const label facei);

            //- Apply the bounding based on user inputs
            void applyBruteForceBounding();

public:

        // DAVOF curved clipping surface (isoAdvectionCurved.C)

            //- The surface of every interface cell for the state at tState
            //  (exactSphere: the translated sphere; paraboloid: the
            //  six-coefficient fit to the cell's face fractions)
            void buildCellSurfaces(const scalar tState);

            //- The volume flux and the face fraction of every face of the
            //  interface cells and their neighbours from the surface of the
            //  upwind cell (second ring: the closest neighbouring surface)
            void curvedFlux(const scalar dt);

            //- tensorCentroid anchored: for every interface cell the plane
            //  in the frame of the used normal through the cut centroid
            //  re-anchored by the volume rule's cell volume (the fan
            //  polyhedron) to alpha_k V, the point o + c e3 handed to the
            //  scheme as the centroid of the ring Monge fit
            void anchorTensorCentroids();

            //- The Info line of the step
            void reportCurved() const;

            //- fluxIntegration sweptFace (isoAdvectionSwept.C): the liquid
            //  volume vLiq of the volume swept by facei between its traced
            //  pre-image and itself under the PLANE S of the cell donor
            //  (positive; the caller gives it the sign of phiF), the
            //  pre-image's face fraction af and wetted area vector aWet in
            //  the face's orientation. sweeps false (the displacement
            //  parallel to the face): the fraction only, no volume, no
            //  orientation count. UInterp interpolates U in a cell, Up is
            //  the volPointInterpolation of U at the vertices.
            void sweptFaceFlux
            (
                const label facei,
                const label donor,
                const davof::cellSurface& S,
                const scalar dt,
                const scalar phiF,
                const bool sweeps,
                const interpolationCellPoint<vector>& UInterp,
                const pointVectorField& Up,
                scalar& af,
                vector& aWet,
                scalar& vLiq,
                pointField* qOut = nullptr,
                point* xqOut = nullptr
            );

            //- fluxIntegration sweptFaceCells (isoAdvectionSwept.C,
            //  davofSweptPieces.H): the same arguments and results as
            //  sweptFaceFlux, the swept volume cut into the pieces inside
            //  the donor and its vertex-ring cells, each wet in total (a
            //  full cell), dry (an empty cell) or clipped by that cell's
            //  own plane (an interface cell); the face fraction from the
            //  same clip of the pre-image fan
            void sweptFaceCellsFlux
            (
                const label facei,
                const label donor,
                const davof::cellSurface& S,
                const scalar dt,
                const scalar phiF,
                const bool sweeps,
                const interpolationCellPoint<vector>& UInterp,
                const pointVectorField& Up,
                scalar& af,
                vector& aWet,
                scalar& vLiq,
                pointField* qOut = nullptr,
                point* xqOut = nullptr
            );

            //- The traced pre-image (q, xq) and the outward triangle list T
            //  of the swept volume of facei (shared by the two swept
            //  fluxes); preSign is the orientation of the stored pre-image
            //  triangles against the face's own order; counts the flips,
            //  the closure residual and the out-of-donor faces
            void sweptTriangles
            (
                const label facei,
                const label donor,
                const scalar dt,
                const scalar phiF,
                const bool sweeps,
                const interpolationCellPoint<vector>& UInterp,
                const pointVectorField& Up,
                List<FixedList<point, 3>>& T,
                pointField& q,
                point& xq,
                label& preSign,
                scalar& vTotal
            );

            //- The class of the cells the pre-image points q (vertices)
            //  and xq (centre) of a face of the cell donor reach beyond
            //  the donor (sweptClassify_; the codes above)
            label classifyPreImage
            (
                const label donor,
                const pointField& q,
                const point& xq
            ) const;


public:

    //- Runtime type information
    TypeName("isoAdvection");

    //- Constructors

        //- Construct given alpha, phi and velocity field. Note: phi should be
        // divergence free up to a sufficient tolerance
        isoAdvection
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~isoAdvection() = default;


    // Member functions

        //- Advect the free surface. Updates alpha field, taking into account
        //  multiple calls within a single time step.
        template < class SpType, class SuType >
        void advect(const SpType& Sp, const SuType& Su);

        //- advection of the interface
        virtual void advect()
        {
            zeroField Sp;
            zeroField Su;
            advect<zeroField,zeroField>(Sp,Su);
        }

        //- advection of the interface with source terms
        virtual void advect
        (
            const volScalarField::Internal& Sp,
            const volScalarField::Internal& Su
        )
        {
            advect<volScalarField::Internal,volScalarField::Internal>(Sp,Su);
        }


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End advection Foam
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "isoAdvectionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
