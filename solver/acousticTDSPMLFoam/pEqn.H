
/*fvScalarMatrix paEqn
(
    fvm::d2dt2(pa) - sqr(c0)*fvc::laplacian(pa)
);

if (solvePrimaryRegion)
{
    paEqn.relax();
    paEqn.solve();
    solve (fvm::ddt(rho, U) == - fvc::grad(pa));
}
else
{
    pa.correctBoundaryConditions();
    U.correctBoundaryConditions();
}*/

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
    // implicitly add the TC0 damping term to the fvmatrix
    // Info << "Adding TC0 to the diagonal of uEqn: Start \n" << endl;
    //Info << "##TC0.component(tensor::XX) type" << TC0.component(tensor::XX).type() << endl;
    //Info << "##uEqn.diag()[0]: " << uEqn.diag()[0] << endl;
    // uEqn.diag().component(0).cref()+=TC0.component(tensor::XX).ref().primitiveFieldRef();
    // uEqn.diag().component(1).ref()+=TC0.component(tensor::YY).ref().primitiveFieldRef();
    // uEqn.diag().component(2).ref()+=TC0.component(tensor::ZZ).ref().primitiveFieldRef();
    // Info << "Adding TC0 to the diagonal of uEqn: End \n" << endl;
    //uEqn.relax();
    uEqn.solve();
    //solve(uEqn);
    //Info << "Solving uEqn done \n" << endl;
    //uEqn.solve();

    // eq.12 in Kaltenbacher2013
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


