{
    const auto updateCoefficients =
    [&]
    (
        const tensor& sigmaValue,
        scalar& sc0,
        scalar& sc1,
        scalar& sc2,
        tensor& tc0,
        tensor& tc1,
        tensor& tc2
    )
    {
        const scalar dx = sigmaValue.xx();
        const scalar dy = sigmaValue.yy();
        const scalar dz = sigmaValue.zz();

        tc0 = tensor(dx, 0, 0, 0, dy, 0, 0, 0, dz);
        tc1 = tensor
        (
            dx - dy - dz, 0, 0,
            0, dy - dx - dz, 0,
            0, 0, dz - dx - dy
        );
        tc2 = tensor
        (
            dy*dz, 0, 0,
            0, dx*dz, 0,
            0, 0, dx*dy
        );

        const scalar invCgSqr = 1/sqr(cg.value());
        sc0 = (dx + dy + dz)*invCgSqr;
        sc1 = (dx*dy + dx*dz + dy*dz)*invCgSqr;
        sc2 = dx*dy*dz*invCgSqr;
    };

    forAll(sigma, celli)
    {
        updateCoefficients
        (
            sigma[celli],
            SC0[celli],
            SC1[celli],
            SC2[celli],
            TC0[celli],
            TC1[celli],
            TC2[celli]
        );
    }

    forAll(sigma.boundaryField(), patchi)
    {
        forAll(sigma.boundaryField()[patchi], facei)
        {
            updateCoefficients
            (
                sigma.boundaryField()[patchi][facei],
                SC0.boundaryFieldRef()[patchi][facei],
                SC1.boundaryFieldRef()[patchi][facei],
                SC2.boundaryFieldRef()[patchi][facei],
                TC0.boundaryFieldRef()[patchi][facei],
                TC1.boundaryFieldRef()[patchi][facei],
                TC2.boundaryFieldRef()[patchi][facei]
            );
        }
    }
}
