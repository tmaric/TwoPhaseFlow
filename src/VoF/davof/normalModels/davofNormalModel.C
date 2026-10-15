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

\*---------------------------------------------------------------------------*/

#include "davofNormalModel.H"
#include "davof.H"
#include "zoneDistribute.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
inline namespace twoPhaseFlow
{
    defineTypeNameAndDebug(davofNormalModel, 0);
    defineRunTimeSelectionTable(davofNormalModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofNormalModel::davofNormalModel
(
    const word& type,
    const dictionary& dict,
    const reconstruction::davof& scheme,
    const word& defaultApplyTo
)
:
    scheme_(scheme),
    mesh_(scheme.mesh()),
    applyTo_(dict.getOrDefault<word>("applyTo", defaultApplyTo)),
    sliverTol_(dict.getOrDefault<scalar>("sliverTol", 0.1)),
    maxAngle_(dict.getOrDefault<scalar>("maxAngle", 15)),
    cosMaxAngle_(Foam::cos(degToRad(maxAngle_))),
    noNormalFallback_(dict.getOrDefault<word>("noNormalFallback", "keep")),
    mRingPtr_(nullptr),
    mRingFilled_(),
    nInterface_(0),
    nSliver_(0),
    nReplaced_(0),
    nKept_(0),
    nFallback_(0),
    nRingFallback_(0)
{
    if (applyTo_ != "all" && applyTo_ != "check")
    {
        FatalIOErrorInFunction(dict)
            << "normalModel " << type << ": unknown applyTo '" << applyTo_
            << "'. Valid: all, check." << exit(FatalIOError);
    }
    if (noNormalFallback_ != "keep" && noNormalFallback_ != "ringArea")
    {
        FatalIOErrorInFunction(dict)
            << "normalModel " << type << ": unknown noNormalFallback '"
            << noNormalFallback_ << "'. Valid: keep, ringArea."
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::davofNormalModel> Foam::davofNormalModel::New
(
    const dictionary& dict,
    const reconstruction::davof& scheme
)
{
    word modelType(dict.getOrDefault<word>("normalModel", ""));

    // The dictionary the model reads: the scheme's, or a copy carrying the
    // mapping of the deprecated keys.
    dictionary mapped;
    const dictionary* dictPtr = &dict;

    if (modelType.empty())
    {
        const word sliverNormal
        (
            dict.getOrDefault<word>("sliverNormal", "none")
        );
        const word normalSource
        (
            dict.getOrDefault<word>("normalSource", "davof")
        );

        modelType = "cellLocal";
        word applyTo;
        bool sliversOnly = false;

        if (normalSource == "exact")
        {
            modelType = "exact";
            applyTo = "all";
        }
        else if (normalSource == "plicRDF")
        {
            modelType = "plicRDF";
            applyTo = "all";
        }
        else if (normalSource != "davof")
        {
            FatalIOErrorInFunction(dict)
                << "Unknown normalSource '" << normalSource
                << "'. Valid: davof, exact, plicRDF (deprecated: use"
                << " normalModel)." << exit(FatalIOError);
        }
        else if (sliverNormal == "Youngs")
        {
            modelType = "youngs";
            applyTo = "check";
            sliversOnly = true;
        }
        else if (sliverNormal == "plicRDF")
        {
            modelType = "plicRDF";
            applyTo = "check";
            sliversOnly = true;
        }
        else if (sliverNormal != "none")
        {
            FatalIOErrorInFunction(dict)
                << "Unknown sliverNormal '" << sliverNormal
                << "'. Valid: none, Youngs, plicRDF (deprecated: use"
                << " normalModel)." << exit(FatalIOError);
        }

        if (modelType != "cellLocal")
        {
            mapped = dict;
            mapped.set("normalModel", modelType);
            mapped.set("applyTo", applyTo);
            if (sliversOnly)
            {
                // The old rule replaced the slivers only: no angle check
                mapped.set("maxAngle", scalar(180));
            }
            dictPtr = &mapped;
            Info<< "davof: the keys sliverNormal/normalSource are deprecated;"
                << " sliverNormal " << sliverNormal << ", normalSource "
                << normalSource << " mapped onto normalModel " << modelType
                << ", applyTo " << applyTo
                << (sliversOnly ? " (slivers only, maxAngle 180)" : "")
                << endl;
        }
    }

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "normalModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<davofNormalModel>(ctorPtr(*dictPtr, scheme));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::davofNormalModel::isSliver
(
    const label celli,
    const scalar mMag
) const
{
    return mMag < sliverTol_*Foam::pow(mesh_.V()[celli], 2.0/3.0);
}


void Foam::davofNormalModel::fillInterfaceNormals
(
    const reconstruction::davof& scheme,
    volVectorField& mInterface,
    labelList& filled
)
{
    const vectorField& m = scheme.areaNormal().primitiveField();
    const DynamicField<label>& labels = scheme.interfaceLabels();

    vectorField& mI = mInterface.primitiveFieldRef();
    for (const label celli : filled)
    {
        mI[celli] = Zero;
    }
    forAll(labels, i)
    {
        mI[labels[i]] = m[labels[i]];
    }
    filled = labels;
}


void Foam::davofNormalModel::ringAreaNormals
(
    const reconstruction::davof& scheme,
    volVectorField& mInterface,
    labelList& filled,
    vectorField& n
)
{
    fillInterfaceNormals(scheme, mInterface, filled);

    // The stencil and the exchange pattern are set up at every call:
    // zoneCPCStencil keeps the point-neighbour stencil of every cell it
    // has computed (uptodate_), so updateStencil costs one scan of the
    // zone plus the new cells, and the pattern (one all-to-all in
    // parallel) must follow the zone, which changes every step; it is
    // not cached across calls because zoneDistribute is one object per
    // mesh, shared with the helper plicRDF of the plicRDF model (its own
    // setUpCommforZone, with its own zone, can run in between).
    const fvMesh& mesh = scheme.mesh();
    zoneDistribute& exchange = zoneDistribute::New(mesh);
    const boolList& zone = scheme.interfaceCell();
    exchange.setUpCommforZone(zone, true);

    const Map<vector> mapM
    (
        exchange.getDatafromOtherProc(zone, mInterface)
    );
    const labelListList& stencil = exchange.getStencil();
    const DynamicField<label>& labels = scheme.interfaceLabels();

    n.resize(labels.size());
    forAll(labels, i)
    {
        const label celli = labels[i];
        vector sum(Zero);
        for (const label gblIdx : stencil[celli])
        {
            sum += exchange.getValue(mInterface, mapM, gblIdx);
        }
        n[i] = sum;
    }
}


void Foam::davofNormalModel::select
(
    const scalar tState,
    const vectorField& m,
    vectorField& nUse,
    boolList& changed
)
{
    const DynamicField<label>& labels = scheme_.interfaceLabels();

    nUse.resize(labels.size());
    changed.resize(labels.size());
    nInterface_ = labels.size();
    nSliver_ = 0;
    nReplaced_ = 0;
    nKept_ = 0;
    nFallback_ = 0;
    nRingFallback_ = 0;

    vectorField nModel(labels.size(), Zero);
    if (replaces())
    {
        modelNormals(tState, nModel);
    }

    // The ring normals of the fallback, computed up front (the exchange is
    // collective: every processor takes part whether or not it has a cell
    // without a model normal)
    const bool ringFallback = (noNormalFallback_ == "ringArea");
    vectorField nRing;
    if (replaces() && ringFallback)
    {
        if (!mRingPtr_)
        {
            mRingPtr_.reset
            (
                new volVectorField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "davofRingFallback",
                            scheme_.alpha1().group()
                        ),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        IOobject::NO_REGISTER
                    ),
                    mesh_,
                    dimensionedVector(dimArea, Zero)
                )
            );
        }
        ringAreaNormals(scheme_, *mRingPtr_, mRingFilled_, nRing);
    }

    const bool all = (applyTo_ == "all");

    forAll(labels, i)
    {
        const label celli = labels[i];
        const scalar mMag = mag(m[celli]);
        const bool sliver = isSliver(celli, mMag);
        if (sliver)
        {
            ++nSliver_;
        }

        nUse[i] = m[celli];
        changed[i] = false;

        if (!replaces())
        {
            continue;
        }

        const vector& g = nModel[i];
        const scalar gMag = mag(g);
        const bool usable = (gMag > 0 && gMag < ROOTVGREAT);

        bool wanted = all || sliver;
        if (!wanted && usable && mMag > 0)
        {
            const scalar cosAngle = (m[celli] & g)/(mMag*gMag);
            wanted = (cosAngle < cosMaxAngle_);
        }
        if (!wanted)
        {
            continue;
        }

        if (usable)
        {
            nUse[i] = g;
            changed[i] = true;
            ++nReplaced_;
            continue;
        }

        // No model normal: the ring normal of the cell, or m_c
        if (ringFallback)
        {
            const vector& r = nRing[i];
            const scalar rMag = mag(r);
            if (rMag > 0 && rMag < ROOTVGREAT)
            {
                nUse[i] = r;
                changed[i] = true;
                ++nRingFallback_;
                continue;
            }
        }
        ++nKept_;
    }

    nInterface_ = returnReduce(nInterface_, sumOp<label>());
    nSliver_ = returnReduce(nSliver_, sumOp<label>());
    nReplaced_ = returnReduce(nReplaced_, sumOp<label>());
    nKept_ = returnReduce(nKept_, sumOp<label>());
    nFallback_ = returnReduce(nFallback_, sumOp<label>());
    nRingFallback_ = returnReduce(nRingFallback_, sumOp<label>());
}


void Foam::davofNormalModel::report(Ostream& os) const
{
    os  << "normalModel " << type() << " (applyTo " << applyTo_
        << ", sliverTol " << sliverTol_;
    if (applyTo_ == "check")
    {
        os  << ", maxAngle " << maxAngle_;
    }
    if (noNormalFallback_ != "keep")
    {
        os  << ", noNormalFallback " << noNormalFallback_;
    }
    os  << "): " << nSliver_ << " slivers, " << nReplaced_ << " replaced";
    if (nKept_)
    {
        os  << ", " << nKept_ << " kept m";
    }
    if (nFallback_)
    {
        os  << ", " << nFallback_ << " fit fallbacks";
    }
    if (nRingFallback_)
    {
        os  << ", " << nRingFallback_ << " ring fallbacks";
    }
}


// ************************************************************************* //
