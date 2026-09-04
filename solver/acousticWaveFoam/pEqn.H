while (pimple.correctNonOrthogonal())
{
    if (timeIntegrationMethod == "thetaNewmark")
    {
        const dimensionedScalar deltaT = runTime.deltaT();
        const scalar theta = timeIntegrationTheta;

        // Retain the non-dissipative update in the physical domain and damp
        // unresolved modes with the bounded implicit update only in the PML.
        fvScalarMatrix vEqn
        (
            fvm::ddt(v) - p
        );
        vEqn.solve();

        fvVectorMatrix uEqn
        (
            fvm::ddt(u)
          + (TC0 & u)
          + (T1 & fvc::grad(p))
          - (T2 & fvc::grad(v))
        );
        uEqn.solve();

        fvScalarMatrix pEqn
        (
            (1.0/(theta*deltaT))
           *fvm::ddt(physicalTimeWeight*rho*compressibility, p)
          - (1.0/(theta*deltaT))
           *physicalTimeWeight*rho*compressibility*pDot.oldTime()
          + fvm::d2dt2
            (
                pmlImplicitWeight*rho*compressibility,
                p
            )
          - rho*fvm::laplacian
            (
                invRhof
               *(theta*physicalTimeWeightf + pmlImplicitWeightf),
                p,
                "laplacian(invRhof,p)"
            )
          + (1.0 - theta)*previousPressureLaplacian()
          + fvm::ddt(C0, p)
          + fvm::Sp(C1, p)
          + C2*v
          - fvc::div(u, "div(u)")
        );
        pEqn.solve();
    }
    else
    {
        // Fully implicit reference update.
        fvScalarMatrix vEqn
        (
            fvm::ddt(v) - p
        );
        vEqn.solve();

        fvVectorMatrix uEqn
        (
            fvm::ddt(u)
          + (TC0 & u)
          + (T1 & fvc::grad(p))
          - (T2 & fvc::grad(v))
        );
        uEqn.solve();

        fvScalarMatrix pEqn
        (
            fvm::d2dt2(rho*compressibility, p)
          - rho*fvm::laplacian(invRhof, p)
          + fvm::ddt(C0, p)
          + fvm::Sp(C1, p)
          + C2*v
          - fvc::div(u)
        );
        pEqn.solve();
    }
}
