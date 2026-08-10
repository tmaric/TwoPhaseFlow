while (pimple.correctNonOrthogonal())
{
    if (laplacianCorrection > 0.5)
    {
        const dimensionedScalar deltaT = runTime.deltaT();
        const dimensionedTensor identity
        (
            "identity",
            dimless,
            tensor::I
        );
        const scalar theta = timeIntegrationTheta;

        tmp<volScalarField> pTheta =
            theta*p + (1.0 - theta)*p.oldTime();

        // Theta updates of the PML memory variables.
        v == v.oldTime() + deltaT*pTheta();
        tmp<volScalarField> vTheta =
            theta*v + (1.0 - theta)*v.oldTime();

        tmp<volTensorField> pmlInv = inv(identity + theta*deltaT*TC0);
        tmp<volVectorField> uTheta = pmlInv()
          &
            (
                u.oldTime()
              - theta*deltaT*(T1 & fvc::grad(pTheta()))
              + theta*deltaT*(T2 & fvc::grad(vTheta()))
            );
        u ==
            (uTheta() - (1.0 - theta)*u.oldTime())/theta;

        // Theta=0.5 is the non-dissipative midpoint/Newmark update. A small
        // off-centering above 0.5 damps unresolved modes in the PML.
        fvScalarMatrix pEqn
        (
            (1.0/(theta*deltaT))*fvm::ddt(invc2/rho, p)
          - (1.0/(theta*deltaT))*(invc2/rho)*pDot.oldTime()
          - theta*fvm::laplacian(invRhof, p)
          - (1.0 - theta)*fvc::laplacian(invRhof, p.oldTime())
          + fvm::ddt(C0/rho, p)
          + theta*fvm::Sp(C1/rho, p)
          + (1.0 - theta)*(C1/rho)*p.oldTime()
          + (C2/rho)*vTheta()
          - fvc::div(uTheta(), "div(u)")/rho
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
            fvm::d2dt2(invc2/rho, p)
          - fvm::laplacian(invRhof, p)
          + fvm::ddt(C0/rho, p)
          + fvm::Sp(C1/rho, p)
          + (C2/rho)*v

          - fvc::div(u)/rho
        );
        pEqn.solve();
    }
}
