while (pimple.correctNonOrthogonal())
{
    // eq.14 in Kaltenbacher2013
    fvScalarMatrix vEqn //solve
    (
        fvm::ddt(v) - p
    );
    //vEqn.relax();
    vEqn.solve();
    //Info << "Solving vEqn done \n" << endl;
    // eq.13 in Kaltenbacher2013
    fvVectorMatrix uEqn
    (
        fvm::ddt(u) 
      + (TC0 & u) 
      + (TC1 & fvc::grad(p)) 
      - (TC2 & fvc::grad(v))
    );
    uEqn.solve();

    fvScalarMatrix pEqn
    (
        fvm::d2dt2( (1 + ((kl - kg)/kg) * alpha1)/sqr(cg), p)
      + fvm::ddt( SC0, p )
      + fvm::Sp( SC1, p )
      + SC2*v 
      + fvm::laplacian((rhol - rhog)/rhol*alphaf - 1, p) 
      - fvc::div(u) 
      //+ sqr(c0)*(rhol - rho)/rhol*(fvc::grad(alpha1) & fvc::grad(p))
     );
     pEqn.relax();
     pEqn.solve();

     if (pimple.finalNonOrthogonalIter())
    {
        p.relax();
        solve (fvm::ddt(rho, U) == - fvc::grad(p));
    }
}


