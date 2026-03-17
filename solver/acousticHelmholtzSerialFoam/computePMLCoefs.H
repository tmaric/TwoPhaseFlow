{
    const auto& Ccs = mesh.C();

    forAll(Ccs, i)
    {
        const scalar Omega = 2*Foam::constant::mathematical::pi*f.value();

        switch (selectedModel)
        {
            case RECTANGLE:
            {
                vector DpFunc(Zero);

                for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
                {
                    if
                    (
                        (PMLMax.component(cmpt) - Ccs[i].component(cmpt)) < L
                     || (Ccs[i].component(cmpt) - PMLMin.component(cmpt)) < L
                    )
                    {
                        const scalar distCmpt =
                        (
                            (
                                (PMLMax.component(cmpt) - Ccs[i].component(cmpt)) < L
                            )
                          ? (PMLMax.component(cmpt) - Ccs[i].component(cmpt))
                          : (Ccs[i].component(cmpt) - PMLMin.component(cmpt))
                        ).value();

                        scalar s = (L.value() - distCmpt)/L.value();
                        s = min(max(s, scalar(0)), scalar(1));
                        DpFunc.component(cmpt) =
                            sigmaMax.value()*Foam::pow(s, po.value());
                    }
                }

                const scalar dx = DpFunc.x();
                const scalar dy = DpFunc.y();
                const scalar sigmaLocal = max(dx, dy);

                TC1[i] = 0;
                SC0[i] = sqr(sigmaLocal)/sqr(cg).value();
                SC1[i] = (2*Omega*sigmaLocal)/sqr(cg).value();
                sigma[i] = sigmaLocal;
                
                break;
            }

            case SPHERICAL:
            {
                const scalar Dist = mag(Ccs[i] - origin.value());
                scalar sigmaLocal = 0.0;

                if ((Rmax.value() - Dist) < L.value())
                {
                    scalar s = (L.value() - (Rmax.value() - Dist))/L.value();
                    s = min(max(s, scalar(0)), scalar(1));
                    sigmaLocal = sigmaMax.value()*Foam::pow(s, po.value());
                }

                sigma[i] = sigmaLocal;

                const scalar denom = sqr(sigmaLocal) + sqr(Omega);
                const scalar sigmaCubeOverOmega = (mag(Omega) > VSMALL)
                    ? (sigmaLocal*sigmaLocal*sigmaLocal/Omega)
                    : 0.0;
                const scalar tc1 = (sigmaLocal > 0 || Omega > 0)
                    ? (sigmaCubeOverOmega + Omega*sigmaLocal)/(denom + VSMALL)
                    : 0.0;

                TC1[i] = tc1;
                SC0[i] = (3*sigmaLocal*sigmaLocal)/sqr(cg).value();
                SC1[i] = (Omega*(3*sigmaLocal) - sigmaCubeOverOmega)/sqr(cg).value();
                break;
            }

            case UNDEFINED:
            {
                FatalErrorInFunction
                    << "Internal error: PML type is undefined."
                    << exit(FatalError);
                break;
            }
        }
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

            switch (selectedModel)
            {
                case RECTANGLE:
                {
                    vector DpFunc(Zero);

                    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
                    {
                        if
                        (
                            (PMLMax.component(cmpt) - Cf.boundaryField()[pi][fi].component(cmpt)) < L
                         || (Cf.boundaryField()[pi][fi].component(cmpt) - PMLMin.component(cmpt)) < L
                        )
                        {
                            const scalar distCmpt =
                            (
                                (
                                    (PMLMax.component(cmpt) - Cf.boundaryField()[pi][fi].component(cmpt)) < L
                                )
                              ? (PMLMax.component(cmpt) - Cf.boundaryField()[pi][fi].component(cmpt))
                              : (Cf.boundaryField()[pi][fi].component(cmpt) - PMLMin.component(cmpt))
                            ).value();

                            scalar s = (L.value() - distCmpt)/L.value();
                            s = min(max(s, scalar(0)), scalar(1));
                            DpFunc.component(cmpt) =
                                sigmaMax.value()*Foam::pow(s, po.value());
                        }
                    }

                    const scalar dx = DpFunc.x();
                    const scalar dy = DpFunc.y();
                    const scalar sigmaLocal = max(dx, dy);

                    TC1.boundaryFieldRef()[pi][fi] = 0;
                    SC0.boundaryFieldRef()[pi][fi] =
                        sqr(sigmaLocal)/sqr(cg).value();
                    SC1.boundaryFieldRef()[pi][fi] =
                        (2*Omega*sigmaLocal)/sqr(cg).value();
                    sigma.boundaryFieldRef()[pi][fi] = sigmaLocal;
                    break;
                }

                case SPHERICAL:
                {
                    const scalar Dist = mag(Cf.boundaryField()[pi][fi] - origin.value());
                    scalar sigmaLocal = 0.0;

                    if ((Rmax.value() - Dist) < L.value())
                    {
                        scalar s = (L.value() - (Rmax.value() - Dist))/L.value();
                        s = min(max(s, scalar(0)), scalar(1));
                        sigmaLocal = sigmaMax.value()*Foam::pow(s, po.value());
                    }

                    sigma.boundaryFieldRef()[pi][fi] = sigmaLocal;

                    const scalar denom = sqr(sigmaLocal) + sqr(Omega);
                    const scalar sigmaCubeOverOmega = (mag(Omega) > VSMALL)
                        ? (sigmaLocal*sigmaLocal*sigmaLocal/Omega)
                        : 0.0;
                    const scalar tc1 = (sigmaLocal > 0 || Omega > 0)
                        ? (sigmaCubeOverOmega + Omega*sigmaLocal)/(denom + VSMALL)
                        : 0.0;

                    TC1.boundaryFieldRef()[pi][fi] = tc1;
                    SC0.boundaryFieldRef()[pi][fi] =
                        (3*sigmaLocal*sigmaLocal)/sqr(cg).value();
                    SC1.boundaryFieldRef()[pi][fi] =
                        (Omega*(3*sigmaLocal) - sigmaCubeOverOmega)/sqr(cg).value();
                    break;
                }

                case UNDEFINED:
                {
                    FatalErrorInFunction
                        << "Internal error: PML type is undefined."
                        << exit(FatalError);
                    break;
                }
            }
        }
    }
}
