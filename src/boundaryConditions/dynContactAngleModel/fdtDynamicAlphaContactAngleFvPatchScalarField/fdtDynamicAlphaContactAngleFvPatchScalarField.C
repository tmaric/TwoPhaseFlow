/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "fdtDynamicAlphaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "scalarField.H"
#include "volMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fdtDynamicAlphaContactAngleFvPatchScalarField::
fdtDynamicAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    theta0_(0.0),
    thetaA_(0.0),
    thetaR_(0.0),
    dxdy_(1.0)
{}


Foam::fdtDynamicAlphaContactAngleFvPatchScalarField::
fdtDynamicAlphaContactAngleFvPatchScalarField
(
    const fdtDynamicAlphaContactAngleFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, p, iF, mapper),
    theta0_(gcpsf.theta0_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_),
    dxdy_(gcpsf.dxdy_)
{}


Foam::fdtDynamicAlphaContactAngleFvPatchScalarField::
fdtDynamicAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF, dict),
    theta0_(dict.get<scalar>("theta0")),
    thetaA_(dict.get<scalar>("thetaA")),
    thetaR_(dict.get<scalar>("thetaR")),
    dxdy_(dict.get<scalar>("dxdy"))
{
    evaluate();
}


Foam::fdtDynamicAlphaContactAngleFvPatchScalarField::
fdtDynamicAlphaContactAngleFvPatchScalarField
(
    const fdtDynamicAlphaContactAngleFvPatchScalarField& gcpsf
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf),
    theta0_(gcpsf.theta0_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_),
    dxdy_(gcpsf.dxdy_)
{}


Foam::fdtDynamicAlphaContactAngleFvPatchScalarField::
fdtDynamicAlphaContactAngleFvPatchScalarField
(
    const fdtDynamicAlphaContactAngleFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, iF),
    theta0_(gcpsf.theta0_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_),
    dxdy_(gcpsf.dxdy_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::fdtDynamicAlphaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{
    // Get the interface normals at the wall.
    const vectorField nf(patch().nf());

    // Find the direction of the interface parallel to the wall
    vectorField nWall(nHat - (nf & nHat)*nf);
    // Normalise nWall
    nWall /= (mag(nWall) + SMALL);

    // Calculate contact line velocity relative to the wall 
    // velocity for moving walls. 
    vectorField Uwall(Up.patchInternalField() - Up);
    // Make contact linen velocity Uwall tangential to the wall.
    Uwall -= (nf & Uwall)*nf;

    // Calculate component of the contact line velocity Uwal in the direction
    // of the interface normal tagential to the wall.
    scalarField uwall(nWall & Uwall);

    // Fetch physical properties
    // TODO(TM): calculate single-field properties?
    // TODO(TM): check if \rho1\nu1 > \rho2\nu2

    const dictionary& transportProperties =
    this->db().objectRegistry::lookupObject<IOdictionary>
    (
        "transportProperties"
    );

    // Choose the larger dynamic viscosity from available two phases.
    word phase1Name (wordList(transportProperties.lookup("phases"))[0]);
    word phase2Name (wordList(transportProperties.lookup("phases"))[1]);

    // Get constant phase-specific densities and kinematic viscosities.
    dimensionedScalar rho1(transportProperties.subDict(phase1Name).get<dimensionedScalar>("rho"));
    dimensionedScalar nu1c(transportProperties.subDict(phase1Name).get<dimensionedScalar>("nu"));

    dimensionedScalar rho2(transportProperties.subDict(phase2Name).get<dimensionedScalar>("rho"));
    dimensionedScalar nu2c(transportProperties.subDict(phase2Name).get<dimensionedScalar>("nu"));

    word nuName; 
    dimensionedScalar rho(rho1); 
    // If the dynamic viscosity of phase1 is larger
    if (rho1*nu1c > rho2*nu2c)
    {
        nuName = "nu1";
        rho = rho1;
    }
    else
    {
        nuName = "nu2";
        rho = rho2;
    }
    const volScalarField& nu =
        this->db().objectRegistry::lookupObject<volScalarField>(nuName);

    // Fetch the wall kinematic viscosity of phase1 
    const label patchi = this->patch().index();
    const fvPatchScalarField& nup = nu.boundaryField()[patchi];
    scalarField muwall (nup*rho.value());

    // Wall Capillary number
    dimensionedScalar sigmap(transportProperties.get<dimensionedScalar>("sigma"));
    scalarField CaWall(muwall*uwall/sigmap.value());

    // Compute the contact angles at the wall.
    tmp<scalarField> thetafTmp = Foam::radToDeg(Foam::acos(nHat & nf));
    scalarField& thetaf = thetafTmp.ref();

    // For all boundary faces
    forAll(thetaf, faceI)
    {
        // If we are in a contact-line cell
        if (mag(nHat[faceI]) > 0)
        {
            if (thetaf[faceI] < thetaR_) // Receding regime
            {
                // TODO(TM): - empirical dynamic contact angle
                //           - see how to make this work with Navier-Slip
                //           - GNBC?
                thetaf[faceI] = thetaR_;
            }
            else if (thetaf[faceI] > thetaA_) // Advancing regime
            {
                // TODO(TM): - empirical dynamic contact angle
                //           - see how to make this work with Navier-Slip
                //           - GNBC?
                thetaf[faceI] = thetaA_;
            }
            else // Hysteresis regime
            {
                // Equation 32 in the manuscript.
                scalar Cstar = dxdy_ * (thetaA_ - thetaR_) * muwall[faceI]  / 
                    (mag(cos(degToRad(thetaA_)) - cos(degToRad(thetaR_)))
                     *sigmap.value());

                // Equation 31 in the manuscript.
                scalar dtheta = Cstar * mag(uwall[faceI]);

                if (uwall[faceI] < 0)
                {
                    thetaf[faceI] -= dtheta;
                }
                if (uwall[faceI] > 0)
                {
                    thetaf[faceI] += dtheta;
                }
            }
        }
    }

    return thetafTmp; 
}


void Foam::fdtDynamicAlphaContactAngleFvPatchScalarField::write(Ostream& os) const
{
    alphaContactAngleTwoPhaseFvPatchScalarField::write(os);
    os.writeEntry("theta0", theta0_);
    os.writeEntry("thetaA", thetaA_);
    os.writeEntry("thetaR", thetaR_);
    os.writeEntry("dxdy", dxdy_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fdtDynamicAlphaContactAngleFvPatchScalarField
    );
}


// ************************************************************************* //
