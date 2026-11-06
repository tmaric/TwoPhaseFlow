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
    initVortex2D

Description
    Writes the velocity U (cell centres) and the volumetric flux phi of the
    2D single vortex at t = 0 in the x-z plane of a one-cell-thick mesh
    (the thickness along y, front and back patches empty),
        u = -sin(2 pi z) sin^2(pi x),   w = sin(2 pi x) sin^2(pi z),
    the stream function psi = (1/pi) sin^2(pi x) sin^2(pi z): the sheared
    disc of W.J. Rider, D.B. Kothe, Reconstructing volume tracking,
    J. Comput. Phys. 141 (1998) 112-152, doi:10.1006/jcph.1998.5906, and
    the field of this repository's testsuite advection/test-vortexShearedDisc
    (generateUVortex2D.C, J. Roenby). The flux of every face is the exact
    line integral of the stream function around the face,
        phi_f = sum over the face's edges of (psi_1 + psi_2)/2 (y_2 - y_1),
    which is divergence-free to round-off on any mesh extruded along y
    (the in-plane edges have y_2 = y_1, the edges along y carry psi at
    their (x, z)). The time factor cos(pi t / T) is applied by advectorVoF
    (spirallingFlow T of the isoAdvector dictionary): U = U0 cos(pi (t +
    dt/2) / T), phi likewise, so the flow reverses at T/2 and the exact
    state at t = T is the initial disc.

    Prints the reference speed s = max over the cells of 0.5 sum_f |phi_f| /
    V_c times V_c^(1/3), so that dt = Co h / s gives OpenFOAM's Courant
    number Co at t = 0 (REFERENCE_SPEED of the leia case isoShear2D), and
    the largest |sum_f phi_f| / V_c.

    Reads U from the time directory (its boundary conditions) and writes U
    and phi there.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

void vortexField(const vectorField& X, vectorField& Uf)
{
    const scalar pi = Foam::constant::mathematical::pi;
    Uf.setSize(X.size());
    forAll(X, i)
    {
        const scalar x = X[i].x(), z = X[i].z();
        Uf[i] = vector
        (
           -Foam::sin(2*pi*z)*sqr(Foam::sin(pi*x)),
            0,
            Foam::sin(2*pi*x)*sqr(Foam::sin(pi*z))
        );
    }
}

scalar streamFunction(const point& p)
{
    const scalar pi = Foam::constant::mathematical::pi;
    return sqr(Foam::sin(pi*p.x()))*sqr(Foam::sin(pi*p.z()))/pi;
}

scalar faceFlux(const face& f, const pointField& P)
{
    scalar phi = 0;
    forAll(f, i)
    {
        const point& p1 = P[f[i]];
        const point& p2 = P[f[(i + 1) % f.size()]];
        phi += 0.5*(streamFunction(p1) + streamFunction(p2))*(p2.y() - p1.y());
    }
    return phi;
}

} // End anonymous namespace


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Write U and phi of the 2D single vortex (Rider & Kothe 1998) at t = 0"
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
    vortexField(mesh.C().primitiveField(), Uc);
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
    const pointField& P = mesh.points();
    const faceList& faces = mesh.faces();
    forAll(phi.primitiveFieldRef(), facei)
    {
        phi.primitiveFieldRef()[facei] = faceFlux(faces[facei], P);
    }
    surfaceScalarField::Boundary& phib = phi.boundaryFieldRef();
    forAll(phib, patchi)
    {
        const label start = mesh.boundaryMesh()[patchi].start();
        forAll(phib[patchi], i)
        {
            phib[patchi][i] = faceFlux(faces[start + i], P);
        }
    }

    // The reference speed of the Courant number at t = 0 and the divergence
    scalar s = 0;
    scalar divMax = 0;
    const labelUList& own = mesh.faceOwner();
    const labelUList& nei = mesh.faceNeighbour();
    scalarField sumAbs(mesh.nCells(), Zero), sumPhi(mesh.nCells(), Zero);
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

    Info<< "initVortex2D: reference speed (max over the cells of"
        << " 0.5 sum_f |phi_f| / V times V^(1/3)) " << s
        << ", max |sum_f phi_f| / V " << divMax
        << ", max |U| " << gMax(mag(U.primitiveField())) << nl
        << "End" << endl;

    return 0;
}


// ************************************************************************* //
