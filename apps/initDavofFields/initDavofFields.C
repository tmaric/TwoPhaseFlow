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

Application
    initDavofFields

Description
    Initialises the DAVOF (dual area/volume-of-fluid) state from an implicit
    function: the cell volume fractions alpha.<phase> and the face liquid-area
    fractions alphaf.<phase> (IOobject::groupName("alphaf", alpha.group())),
    both written to the current time directory.

    Reads system/initAlphaFieldDict as initAlphaField does (field, the
    implicit function `type` and its parameters), plus

        davofMethod  detrixheAslam | quadraticFaces | exactSphere;
        wispTol      1e-6;      // interface cell: |m_c| > wispTol V_c^(2/3)
        signedDistance false;   // Newton-projected signed distance (ellipsoid)
        invert       false;
        faceState    scalar;    // vector: also the wetted area vectors
                                // wettedAreaDavof.<phase> (2026-10-25), the
                                // sum over every face's fan triangles of the
                                // wetted area vectors of the same cuts

    Prints sum(alpha V), the maximum clip of alpha to [0, 1], the number of
    faces that fell back to the Detrixhe-Aslam value (quadraticFaces) and the
    maximum per-cell inconsistency
    |sum_f alphaf S_f^out - mTet_c| / V_c^(2/3) over the interface cells
    (mTet_c the tet zero-set area vector of the Detrixhe-Aslam cut; round-off
    for detrixheAslam, the O(h) chord-to-quadratic / chord-to-exact state
    difference for the other two).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "implicitFunction.H"
#include "davofInit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Initialise the DAVOF state (alpha and the face fractions alphaf)"
        " from the implicit function of system/initAlphaFieldDict"
    );
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    IOdictionary initAlphaFieldDict
    (
        IOobject
        (
            "initAlphaFieldDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word nameField(initAlphaFieldDict.get<word>("field"));
    const word davofMethod(initAlphaFieldDict.get<word>("davofMethod"));
    const scalar wispTol =
        initAlphaFieldDict.getOrDefault<scalar>("wispTol", 1e-6);
    const bool invert = initAlphaFieldDict.getOrDefault<bool>("invert", false);
    const word faceState
    (
        initAlphaFieldDict.getOrDefault<word>("faceState", "scalar")
    );
    if (faceState != "scalar" && faceState != "vector")
    {
        FatalIOErrorInFunction(initAlphaFieldDict)
            << "faceState '" << faceState << "' unknown. Valid: scalar, vector."
            << exit(FatalIOError);
    }

    Info<< "Reading field " << nameField << nl << endl;
    volScalarField alpha1
    (
        IOobject
        (
            nameField,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    const word alphafName(IOobject::groupName("alphaf", alpha1.group()));
    Info<< "Creating face-fraction field " << alphafName << nl << endl;
    surfaceScalarField alphaf
    (
        IOobject
        (
            alphafName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        calculatedFvsPatchScalarField::typeName
    );

    Foam::autoPtr<Foam::implicitFunction> func = implicitFunction::New
    (
        initAlphaFieldDict.get<word>("type"),
        initAlphaFieldDict
    );

    davofInit init(mesh, func(), initAlphaFieldDict);

    Info<< "davofMethod " << davofMethod << endl;
    if (davofMethod == "detrixheAslam" || davofMethod == "linearInterpolant")
    {
        init.detrixheAslam(alpha1, alphaf);
    }
    else if (davofMethod == "quadraticFaces")
    {
        init.quadraticFaces(alpha1, alphaf);
    }
    else if (davofMethod == "exactSphere")
    {
        init.exactSphere(alpha1, alphaf);
    }
    else
    {
        FatalIOErrorInFunction(initAlphaFieldDict)
            << "Unknown davofMethod '" << davofMethod
            << "'. Valid: detrixheAslam, quadraticFaces, exactSphere."
            << exit(FatalIOError);
    }

    // The consistency check before any inversion (mTet is the liquid's).
    const scalar maxConsistency = init.maxConsistency(alphaf, wispTol);

    // faceState vector: the wetted area vectors of the same construction
    autoPtr<surfaceVectorField> wetVecPtr;
    scalar maxConsistencyVec = -1;
    if (faceState == "vector")
    {
        wetVecPtr.reset
        (
            new surfaceVectorField
            (
                IOobject
                (
                    IOobject::groupName("wettedAreaDavof", alpha1.group()),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector(dimArea, Zero),
                calculatedFvsPatchVectorField::typeName
            )
        );
        init.fillWettedArea(*wetVecPtr);
        maxConsistencyVec = init.maxConsistency(alphaf, *wetVecPtr, wispTol);
    }

    if (invert)
    {
        alpha1 = scalar(1) - alpha1;
        alphaf = scalar(1) - alphaf;
        if (wetVecPtr)
        {
            *wetVecPtr = mesh.Sf() - *wetVecPtr;
        }
    }

    // Round trip of the fractions through ASCII at full precision.
    IOstream::defaultPrecision(18);
    alpha1.write();
    alphaf.write();
    if (wetVecPtr)
    {
        wetVecPtr->write();
        Info<< "davof maxConsistency (wetted area vectors) = "
            << maxConsistencyVec << endl;
    }

    const scalarField& alpha = alpha1.internalField();
    Info<< "sum(alpha*V) = " << gSum(mesh.V()*alpha)
        << ", 1-max(alpha1) = " << 1 - gMax(alpha)
        << "\t min(alpha1) = " << gMin(alpha) << nl
        << "davof maxAlphaClip = " << init.maxAlphaClip() << nl
        << "davof nFaceFallback = " << init.nFaceFallback() << nl
        << "davof maxConsistency = " << maxConsistency
        << " (max |sum_f alphaf S_f^out - mTet| / V^(2/3), interface cells, wispTol "
        << wispTol << ")" << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
