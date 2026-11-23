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
    Foam::davofNormalModels::ringSumRDF

Description
    A control model (Tomislav, 2026-11-23): the ring-sum normal followed by
    rdfIterations updates of the reconstructed distance function (RDF) of
    plicRDF.

    Step 1: the ring sum d_c = sum_j m_j of every interface cell (the
    normal of the model ringSum), the plane of every interface cell
    positioned by alpha along it as the davof scheme positions its planes
    (davof::positionPlane).

    Step 2, repeated rdfIterations times: the RDF of these planes on the
    point-neighbour stencil (reconstructedDistanceFunction::constructRDF,
    the call of plicRDF::reconstruct()), the least-squares gradient of the
    RDF at the interface cells (leastSquareGrad polyDegree1 over the
    stencil, the body of plicRDF::gradSurf()), normalised; for all but the
    last repetition the plane of every interface cell is re-positioned by
    alpha along the new normal. The normal of the last repetition is
    returned; the scheme positions the plane along it by alpha.

    rdfIterations 0 returns the ring sum itself (the model ringSum).
    No residual tolerance, no too-coarse rule, no fallback: a cell whose
    gradient vanishes has no model normal (the base class's
    noNormalFallback decides, default keep). The plicRDF scheme is not
    modified: plicRDF::gradSurf() is private and reconstruct() starts from
    its own interpolated normals, so the model calls the public classes
    that plicRDF calls (reconstructedDistanceFunction, zoneDistribute,
    leastSquareGrad) in the same order with the same arguments.

        normalModel   ringSumRDF;
        rdfIterations 1;          // 0 | 1 | 2 ...

    Reference:
    \verbatim
        Henning Scheufler, Johan Roenby,
        Accurate and efficient surface reconstruction from volume
        fraction data on general meshes,
        Journal of Computational Physics 383 (2019) 1-23,
        doi:10.1016/j.jcp.2019.01.009, arXiv:1801.05382
    \endverbatim

SourceFiles
    ringSumRDFNormal.C

\*---------------------------------------------------------------------------*/

#ifndef davofNormalModels_ringSumRDF_H
#define davofNormalModels_ringSumRDF_H

#include "davofNormalModel.H"
#include "reconstructedDistanceFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
namespace davofNormalModels
{

/*---------------------------------------------------------------------------*\
                          Class ringSumRDF Declaration
\*---------------------------------------------------------------------------*/

class ringSumRDF
:
    public davofNormalModel
{
    // Private Data

        //- Number of RDF updates (step 2)
        label rdfIterations_;

        //- The RDF object plicRDF uses (mesh singleton)
        reconstructedDistanceFunction& RDF_;


    // Private Member Functions

        //- Position the plane of every interface cell by alpha along n
        void positionPlanes(const vectorField& n);

        //- The least-squares gradient of phi at the interface cells
        //  (plicRDF::gradSurf)
        void gradSurf(const volScalarField& phi, vectorField& g);


public:

    //- Runtime type information
    TypeName("ringSumRDF");


    //- Construct from the davof dictionary and the scheme
    ringSumRDF
    (
        const dictionary& dict,
        const reconstruction::davof& scheme
    );


    //- Destructor
    virtual ~ringSumRDF() = default;


    // Member Functions

        //- The normal of the last RDF update (the ring sum for 0)
        virtual void modelNormals(const scalar tState, vectorField& n);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace davofNormalModels
} // End inline namespace twoPhaseFlow
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
