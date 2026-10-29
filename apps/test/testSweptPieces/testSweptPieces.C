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
    testSweptPieces

Description
    The unit check of davofSweptPieces.H (the piecewise clip of the swept
    volume, fluxIntegration sweptFaceCells) against the Python reference
    leia/docs/davof/davof-swept-volume/study_multicell.py
    (data/multicell.csv, 2026-10-29): the face x = 1 of the unit donor cell
    [0,1]^3, planar and warped, traced backward through the uniform velocity
    (1, 0.6, 0.4) over dt = 0.35 (the pre-image is the exact translation, so
    no trajectory integrator is needed), cut into the pieces of the donor,
    the neighbours y- and z- and the edge neighbour y-z-, each clipped by
    its cell's tangent-plane PLIC of the sphere of radius 2.2 about
    (-1, 0.6, 1.4) (liquid inside). The reference values are the signed
    tetrahedra clipped as convex polyhedra with explicit caps, which the
    closure method reproduced to 1e-16 there. Exits 1 above 1e-12 relative
    to the swept volume.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "davofSweptPieces.H"
#include "FixedList.H"
#include "List.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::noCheckProcessorDirectories();
    argList args(argc, argv, false, false, false);

    const vector U(1, 0.6, 0.4);
    const scalar dt = 0.35;
    const vector d = U*dt;
    const point cSph(-1.0, 0.6, 1.4);
    const scalar rSph = 2.2;

    // the cells about the donor's edge y = 0, z = 0 as half-space lists
    // ((x - x0).n <= 0 inside) and their centres
    struct cellDef { word name; List<davof::halfSpace> hs; point centre; scalar ref[2]; };
    auto H = [](const vector& n, const point& x0) { davof::halfSpace h; h.n = n; h.x0 = x0; return h; };
    const davof::halfSpace Y0 = H(vector(0, 1, 0), point(0, 0, 0));
    const davof::halfSpace Y1 = H(vector(0, -1, 0), point(0, 0, 0));
    const davof::halfSpace Z0 = H(vector(0, 0, 1), point(0, 0, 0));
    const davof::halfSpace Z1 = H(vector(0, 0, -1), point(0, 0, 0));
    // the reference wet volumes of data/multicell.csv: planar, warped
    List<cellDef> cellsDef(4);
    cellsDef[0] = {"donor", List<davof::halfSpace>({Y1, Z1}), point(0.5, 0.5, 0.5), {2.305029066657e-01, 2.236824058813e-01}};
    cellsDef[1] = {"neighbour y-", List<davof::halfSpace>({Y0, Z1}), point(0.5, -0.5, 0.5), {2.881191314660e-02, 2.810220144542e-02}};
    cellsDef[2] = {"neighbour z-", List<davof::halfSpace>({Y1, Z0}), point(0.5, 0.5, -0.5), {4.540402083841e-03, 7.658143502765e-03}};
    cellsDef[3] = {"edge y-z-", List<davof::halfSpace>({Y0, Z0}), point(0.5, -0.5, -0.5), {8.693329138745e-06, 0.0}};
    const scalar totalRef[2] = {2.638639152252e-01, 2.594427508295e-01};

    bool ok = true;
    for (label cfg = 0; cfg < 2; ++cfg)
    {
        List<point> face(4);
        face[0] = point(1, 0, 0); face[1] = point(1, 1, 0);
        face[2] = point(1, 1, 1); face[3] = point(1, 0, 1);
        if (cfg == 1)
        {
            face[0].x() += 0.1; face[2].x() += 0.1; face[1].x() -= 0.1; face[3].x() -= 0.1;
        }
        point xf = 0.25*(face[0] + face[1] + face[2] + face[3]);
        List<point> pre(4);
        forAll(face, i) pre[i] = face[i] - d;
        const point xq = xf - d;
        // the outward list for phi_f > 0 (the prototype's construction)
        List<FixedList<point, 3>> T(16);
        for (label i = 0; i < 4; ++i)
        {
            const label j = (i + 1) % 4;
            T[4*i][0] = xf; T[4*i][1] = face[i]; T[4*i][2] = face[j];
            T[4*i + 1][0] = xq; T[4*i + 1][1] = pre[j]; T[4*i + 1][2] = pre[i];
            T[4*i + 2][0] = face[i]; T[4*i + 2][1] = pre[i]; T[4*i + 2][2] = pre[j];
            T[4*i + 3][0] = face[i]; T[4*i + 3][1] = pre[j]; T[4*i + 3][2] = face[j];
        }
        scalar sixV = 0;
        forAll(T, k) sixV += T[k][0] & (T[k][1] ^ T[k][2]);
        const scalar Vtot = sixV/6;
        Info<< (cfg == 0 ? "oblique, planar face" : "oblique, warped face")
            << ": swept volume " << Vtot << nl;
        scalar total = 0;
        forAll(cellsDef, c)
        {
            const cellDef& cd = cellsDef[c];
            // the tangent-plane PLIC: the normal out of the liquid at the
            // sphere point nearest to the cell centre
            vector n = cd.centre - cSph;
            n /= mag(n);
            davof::halfSpace plic;
            plic.n = n;
            plic.x0 = cSph + rSph*n;
            const davof::pieceResult r = davof::sweptPieceVolume(T, cd.hs, true, plic);
            const scalar rel = mag(r.volume - cd.ref[cfg])/Vtot;
            Info<< "    " << cd.name << ": wet volume " << r.volume << " reference " << cd.ref[cfg]
                << " rel " << rel << " caps " << r.nCaps << " closure " << r.closureResidual
                << (r.ok ? "" : " RANK DEFICIENT") << nl;
            if (rel > 1e-12 || !r.ok) ok = false;
            total += r.volume;
        }
        const scalar relT = mag(total - totalRef[cfg])/Vtot;
        Info<< "    all pieces: " << total << " reference " << totalRef[cfg] << " rel " << relT << nl;
        if (relT > 1e-12) ok = false;
    }
    Info<< (ok ? "testSweptPieces: PASS" : "testSweptPieces: FAIL") << endl;
    return ok ? 0 : 1;
}


// ************************************************************************* //
