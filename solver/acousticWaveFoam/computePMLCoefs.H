{
    const auto updateCoefficients =
    [&]
    (
        const tensor& sigmaValue,
        scalar& c0,
        scalar& c1,
        scalar& c2,
        tensor& tc0,
        tensor& t1,
        tensor& t2
    )
    {
        // The shared frequency-domain field setup may store sigma with a
        // convention-dependent sign. Time-domain damping must be nonnegative.
        const scalar dx = mag(sigmaValue.xx());
        const scalar dy = mag(sigmaValue.yy());
        const scalar dz = mag(sigmaValue.zz());

        tc0 = tensor(dx, 0, 0, 0, dy, 0, 0, 0, dz);
        t1 = tensor
        (
            dx - dy - dz, 0, 0,
            0, dy - dx - dz, 0,
            0, 0, dz - dx - dy
        );
        t2 = tensor
        (
            dy*dz, 0, 0,
            0, dx*dz, 0,
            0, 0, dx*dy
        );

        const scalar invCgSqr = 1/sqr(cg.value());
        c0 = (dx + dy + dz)*invCgSqr;
        c1 = (dx*dy + dx*dz + dy*dz)*invCgSqr;
        c2 = dx*dy*dz*invCgSqr;
    };

    forAll(sigma, celli)
    {
        updateCoefficients
        (
            sigma[celli],
            C0[celli],
            C1[celli],
            C2[celli],
            TC0[celli],
            T1[celli],
            T2[celli]
        );
    }

    forAll(sigma.boundaryField(), patchi)
    {
        forAll(sigma.boundaryField()[patchi], facei)
        {
            updateCoefficients
            (
                sigma.boundaryField()[patchi][facei],
                C0.boundaryFieldRef()[patchi][facei],
                C1.boundaryFieldRef()[patchi][facei],
                C2.boundaryFieldRef()[patchi][facei],
                TC0.boundaryFieldRef()[patchi][facei],
                T1.boundaryFieldRef()[patchi][facei],
                T2.boundaryFieldRef()[patchi][facei]
            );
        }
    }
}
