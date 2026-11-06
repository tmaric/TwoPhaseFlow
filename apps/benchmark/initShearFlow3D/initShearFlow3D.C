/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
    unofficial extension to OpenFOAM.

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
    initShearFlow3D

Description
    Writes the velocity U (cell centres) and the volumetric flux phi (face
    centres, phi_f = U(x_f) . S_f) of the 3D shear flow at t = 0,

        u =  sin(2 pi y) sin^2(pi x),
        v = -sin(2 pi x) sin^2(pi y),
        w =  umax (1 - r/R)^2,   r = sqrt((x - 1/2)^2 + (y - 1/2)^2),

    umax = 1, R = 1/2: the shear test of P. Liovic, M. Rudman, J.-L. Liow,
    D. Lakehal, D. Kothe, A 3D unsplit-advection volume tracking algorithm
    with planarity-preserving interface reconstruction, Comput. Fluids 35
    (2006) 1011-1032, doi:10.1016/j.compfluid.2005.09.003, in the form of
    eq. (4.5) of T. Maric, H. Marschall, D. Bothe, voFoam - A geometrical
    Volume of Fluid algorithm on arbitrary unstructured meshes with local
    dynamic adaptive mesh refinement using OpenFOAM, arXiv:1305.3417 (2013),
    and the field of this repository's testsuite
    advection/test-deformationSphere (generateUDeform.C, J. Roenby). The time
    factor cos(pi t / T) is applied by advectorVoF (spirallingFlow T of the
    isoAdvector dictionary): U = U0 cos(pi (t + dt/2) / T), phi likewise, so
    the flow reverses at T/2 and the exact state at t = T is the initial one.

    Prints the reference speed s = max over the cells of 0.5 sum_f |phi_f| /
    V_c times V_c^(1/3), so that dt = Co h / s gives OpenFOAM's Courant
    number Co at t = 0 (REFERENCE_SPEED of the leia case isoShear3D).

    Reads U from the time directory (its boundary conditions) and writes U and
    phi there.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

void shearField(const vectorField& X, vectorField& Uf)
{
    const scalar pi = Foam::constant::mathematical::pi;
    const scalar umax = 1;
    const scalar R = 0.5;
    Uf.setSize(X.size());
    forAll(X, i)
    {
        const scalar x = X[i].x(), y = X[i].y();
        const scalar r = Foam::sqrt(sqr(x - 0.5) + sqr(y - 0.5));
        Uf[i] = vector
        (
            Foam::sin(2*pi*y)*sqr(Foam::sin(pi*x)),
           -Foam::sin(2*pi*x)*sqr(Foam::sin(pi*y)),
            umax*sqr(1 - r/R)
        );
    }
}

} // End anonymous namespace


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Write U and phi of the 3D shear flow (Liovic et al. 2006) at t = 0"
    );
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading field U (its boundary conditions)\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    vectorField Uc;
    shearField(mesh.C().primitiveField(), Uc);
    U.primitiveFieldRef() = Uc;
    U.correctBoundaryConditions();

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U)
    );

    vectorField Uf;
    shearField(mesh.Cf().primitiveField(), Uf);
    phi.primitiveFieldRef() = Uf & mesh.Sf().primitiveField();

    surfaceScalarField::Boundary& phib = phi.boundaryFieldRef();
    forAll(phib, patchi)
    {
        vectorField Ufp;
        shearField(mesh.Cf().boundaryField()[patchi], Ufp);
        phib[patchi] = Ufp & mesh.Sf().boundaryField()[patchi];
    }

    // The reference speed of the Courant number at t = 0
    scalar s = 0;
    scalar divMax = 0;
    const labelUList& own = mesh.faceOwner();
    const labelUList& nei = mesh.faceNeighbour();
    scalarField sumAbs(mesh.nCells(), Zero), sumPhi(mesh.nCells(), Zero);
    // the internal faces only (faceOwner lists every face, faceNeighbour the
    // internal ones: the loop over every face indexed nei beyond its size)
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const scalar p = phi.primitiveField()[facei];
        sumAbs[own[facei]] += mag(p); sumAbs[nei[facei]] += mag(p);
        sumPhi[own[facei]] += p; sumPhi[nei[facei]] -= p;
    }
    forAll(phib, patchi)
    {
        const label start = mesh.boundaryMesh()[patchi].start();
        forAll(phib[patchi], i)
        {
            const scalar p = phib[patchi][i];
            sumAbs[own[start + i]] += mag(p); sumPhi[own[start + i]] += p;
        }
    }
    forAll(sumAbs, celli)
    {
        const scalar V = mesh.V()[celli];
        s = max(s, 0.5*sumAbs[celli]/V*Foam::cbrt(V));
        divMax = max(divMax, mag(sumPhi[celli])/V);
    }
    reduce(s, maxOp<scalar>());
    reduce(divMax, maxOp<scalar>());

    U.write();
    phi.write();

    Info<< "initShearFlow3D: reference speed (max over the cells of"
        << " 0.5 sum_f |phi_f| / V times V^(1/3)) " << s
        << ", max |sum_f phi_f| / V " << divMax
        << ", max |U| " << gMax(mag(U.primitiveField())) << nl
        << "End" << endl;

    return 0;
}


// ************************************************************************* //
