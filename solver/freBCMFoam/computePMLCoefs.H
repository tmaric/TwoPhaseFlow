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
                const scalar dz = DpFunc.z();
                const scalar dxdyDzOverOmega =
                    (mag(Omega) > VSMALL) ? (dx*dy*dz/Omega) : 0.0;

                const scalar tc1x = (dxdyDzOverOmega + Omega*(dy + dz - dx))
                    /(sqr(dx) + sqr(Omega) + VSMALL);
                const scalar tc1y = (dxdyDzOverOmega + Omega*(dx + dz - dy))
                    /(sqr(dy) + sqr(Omega) + VSMALL);
                const scalar tc1z = (dxdyDzOverOmega + Omega*(dy + dx - dz))
                    /(sqr(dz) + sqr(Omega) + VSMALL);

                TC1[i] = (tc1x + tc1y + tc1z)/3.0;
                SC0[i] = (dx*dy + dx*dz + dy*dz)/sqr(cg).value();
                SC1[i] = (Omega*(dx + dy + dz) - dxdyDzOverOmega)/sqr(cg).value();
                sigma[i] = max(dx, max(dy, dz));
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
        !isA<emptyPolyPatch>(patchi)
     && !isA<wedgePolyPatch>(patchi)
     && !isA<processorPolyPatch>(patchi)
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
                    const scalar dz = DpFunc.z();
                    const scalar dxdyDzOverOmega =
                        (mag(Omega) > VSMALL) ? (dx*dy*dz/Omega) : 0.0;

                    const scalar tc1x = (dxdyDzOverOmega + Omega*(dy + dz - dx))
                        /(sqr(dx) + sqr(Omega) + VSMALL);
                    const scalar tc1y = (dxdyDzOverOmega + Omega*(dx + dz - dy))
                        /(sqr(dy) + sqr(Omega) + VSMALL);
                    const scalar tc1z = (dxdyDzOverOmega + Omega*(dy + dx - dz))
                        /(sqr(dz) + sqr(Omega) + VSMALL);

                    TC1.boundaryFieldRef()[pi][fi] = (tc1x + tc1y + tc1z)/3.0;
                    SC0.boundaryFieldRef()[pi][fi] =
                        (dx*dy + dx*dz + dy*dz)/sqr(cg).value();
                    SC1.boundaryFieldRef()[pi][fi] =
                        (Omega*(dx + dy + dz) - dxdyDzOverOmega)/sqr(cg).value();
                    sigma.boundaryFieldRef()[pi][fi] = max(dx, max(dy, dz));
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

// Synchronize coupled boundaries (notably processor patches) from internal
// values after coefficient recomputation.
TC1.correctBoundaryConditions();
SC0.correctBoundaryConditions();
SC1.correctBoundaryConditions();
sigma.correctBoundaryConditions();
