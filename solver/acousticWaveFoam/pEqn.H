while (pimple.correctNonOrthogonal())
{
    // eq.14 in Kaltenbacher2013
    fvScalarMatrix vEqn 
    (
        fvm::ddt(v) - p
    );
    vEqn.solve();
    Info << "Solving vEqn done \n" << endl;
    // eq.13 in Kaltenbacher2013
    fvVectorMatrix uEqn
    (
        fvm::ddt(u) 
      + (TC0 & u) 
      + (TC1 & fvc::grad(p)) 
      - (TC2 & fvc::grad(v))
    );
    uEqn.solve();

    // Trapezoidal/Newmark averaging of the Laplacian avoids the artificial
    // modal damping caused by using the Laplacian fully at t^{n+1}.
    fvScalarMatrix pEqn
    (
        fvm::d2dt2(invc2/rho, p)
      - 0.25*fvm::laplacian(invRhof,p)
      - 0.5*fvc::laplacian(invRhof,p.oldTime())
      - 0.25*fvc::laplacian(invRhof,p.oldTime().oldTime())
      + fvm::ddt(SC0/rho, p)
      + fvm::Sp(SC1/rho, p)
      + (SC2/rho)*v

      - fvc::div(u)/rho
      //+ sqr(c0)*(rhol - rho)/rhol*(fvc::grad(alpha1) & fvc::grad(p)));
    );
    pEqn.relax();
     pEqn.solve();

     if (pimple.finalNonOrthogonalIter())
    {
        p.relax();
        solve (fvm::ddt(rho, U) == - fvc::grad(p));
    }
}
