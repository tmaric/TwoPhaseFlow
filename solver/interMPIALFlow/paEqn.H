
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix paEqn
    (
        fvm::d2dt2( 1 + ((kl - kg)/kg) * alpha1, pa) 
      + sqr(cg)*fvm::laplacian((rho1 - rho2)/rho1*alphaf - 1, pa) 
     );
     paEqn.relax();
     paEqn.solve();

     if (pimple.finalNonOrthogonalIter())
    {
        pa.relax();
        solve
        (
            fvm::ddt(rho, Ua)
         ==
            fvc::reconstruct((-fvc::snGrad(pa))*mesh.magSf())
        );

    }
}

