{
    const scalar omega = constant::mathematical::twoPi*f.value();

    const auto updateCoefficients =
    [&]
    (
        const tensor& sigmaValue,
        const scalar alpha,
        scalar& c0,
        scalar& c1,
        tensor& t0,
        tensor& t1
    )
    {
        const scalar dx = sigmaValue.xx();
        const scalar dy = sigmaValue.yy();
        const scalar dz = sigmaValue.zz();
        const scalar invSoundSpeedSqr =
            alpha/sqr(cl.value()) + (1 - alpha)/sqr(cg.value());
        const scalar dxyzOverOmega =
            mag(omega) > VSMALL ? dx*dy*dz/omega : 0;

        t0 = tensor
        (
            (dx*(dy + dz - dx) - dy*dz)/(sqr(omega) + sqr(dx) + VSMALL),
            0, 0, 0,
            (dy*(dx + dz - dy) - dx*dz)/(sqr(omega) + sqr(dy) + VSMALL),
            0, 0, 0,
            (dz*(dx + dy - dz) - dx*dy)/(sqr(omega) + sqr(dz) + VSMALL)
        );
        t1 = tensor
        (
            (dxyzOverOmega + omega*(dy + dz - dx))
           /(sqr(omega) + sqr(dx) + VSMALL),
            0, 0, 0,
            (dxyzOverOmega + omega*(dx + dz - dy))
           /(sqr(omega) + sqr(dy) + VSMALL),
            0, 0, 0,
            (dxyzOverOmega + omega*(dx + dy - dz))
           /(sqr(omega) + sqr(dz) + VSMALL)
        );

        c0 = (dx*dy + dx*dz + dy*dz)*invSoundSpeedSqr;
        c1 =
            (omega*(dx + dy + dz) - dxyzOverOmega)*invSoundSpeedSqr;
    };

    forAll(sigma, celli)
    {
        updateCoefficients
        (
            sigma[celli],
            min(max(alpha1[celli], scalar(0)), scalar(1)),
            C0[celli],
            C1[celli],
            T0[celli],
            T1[celli]
        );
    }

    forAll(sigma.boundaryField(), patchi)
    {
        forAll(sigma.boundaryField()[patchi], facei)
        {
            updateCoefficients
            (
                sigma.boundaryField()[patchi][facei],
                min
                (
                    max(alpha1.boundaryField()[patchi][facei], scalar(0)),
                    scalar(1)
                ),
                C0.boundaryFieldRef()[patchi][facei],
                C1.boundaryFieldRef()[patchi][facei],
                T0.boundaryFieldRef()[patchi][facei],
                T1.boundaryFieldRef()[patchi][facei]
            );
        }
    }
}
