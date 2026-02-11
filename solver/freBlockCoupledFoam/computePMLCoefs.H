{
    const auto& Ccs = mesh.C();

    forAll(Ccs, i)
    {
        const scalar Omega = 2*Foam::constant::mathematical::pi*f.value();
        const scalar Dist = mag(Ccs[i] - origin.value());

        scalar sigmaLocal = 0.0;
        if ((Rmax.value() - Dist) < L.value())
        {
            scalar s = (L.value() - (Rmax.value() - Dist))/L.value();
            s = min(max(s, scalar(0)), scalar(1));
            sigmaLocal = sigmaMax.value()*Foam::pow(s, p.value());
        }

        sigma[i] = sigmaLocal;

        const scalar omega2 = sqr(Omega) + sqr(alphaPml.value());
        const scalar denom = sqr(sigmaLocal) + omega2;
        const scalar sigmaCubeOverOmega = (omega2 > VSMALL)
            ? (sigmaLocal*sigmaLocal*sigmaLocal*Omega/omega2)
            : 0.0;
        const scalar tc1 = (sigmaLocal > 0 || Omega > 0)
            ? (sigmaCubeOverOmega + Omega*sigmaLocal)/denom
            : 0.0;

        TC1[i] = tc1;

        SC0[i] = (3*sigmaLocal*sigmaLocal)/sqr(cg).value();
        SC1[i] =
            (Omega*(3*sigmaLocal) - sigmaCubeOverOmega)/sqr(cg).value();
    }
}

//-----> compute coefficients in boundary field
const polyBoundaryMesh& patches = mesh.boundaryMesh();
const surfaceVectorField& Cf = mesh.Cf();

forAll(patches, pi)
{
    const polyPatch& patchi = patches[pi];
    Info << "patch name: " << patchi.name() << endl;
    if
    (
        !isA<emptyPolyPatch>(patchi) && !isA<wedgePolyPatch>(patchi)
    )
    {
        forAll(patchi, fi)
        {
            const scalar Omega = 2*Foam::constant::mathematical::pi*f.value();
            const scalar Dist = mag(Cf.boundaryField()[pi][fi] - origin.value());

            scalar sigmaLocal = 0.0;
            if ((Rmax.value() - Dist) < L.value())
            {
                scalar s = (L.value() - (Rmax.value() - Dist))/L.value();
                s = min(max(s, scalar(0)), scalar(1));
                sigmaLocal = sigmaMax.value()*Foam::pow(s, p.value());
            }

            sigma.boundaryFieldRef()[pi][fi] = sigmaLocal;

            const scalar omega2 = sqr(Omega) + sqr(alphaPml.value());
            const scalar denom = sqr(sigmaLocal) + omega2;
            const scalar sigmaCubeOverOmega = (omega2 > VSMALL)
                ? (sigmaLocal*sigmaLocal*sigmaLocal*Omega/omega2)
                : 0.0;
            const scalar tc1 = (sigmaLocal > 0 || Omega > 0)
                ? (sigmaCubeOverOmega + Omega*sigmaLocal)/denom
                : 0.0;

            TC1.boundaryFieldRef()[pi][fi] = tc1;

            SC0.boundaryFieldRef()[pi][fi] = (3*sigmaLocal*sigmaLocal)/sqr(cg).value();
            SC1.boundaryFieldRef()[pi][fi] =
                (Omega*(3*sigmaLocal) - sigmaCubeOverOmega)/sqr(cg).value();
        }
    }
}
