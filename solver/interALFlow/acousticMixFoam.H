{
    #include "computeAlphaf.H"

        while (pimple.correct())
        {
            #include "paEqn.H"
        }

	pr == 1.0/(2.0*rho*sqr(cg))*pa*pa - rho/2.0*(Ua&Ua);
	momFlux == rho*(Ua*Ua);
        //iter=iter+1.0;
        //AcousticForcePotential == (iter - 1)/iter * AcousticForcePotential + 1/iter * (2*Foam::constant::mathematical::pi*pow(R,3)*( 1/(3.0*rho*sqr(cg))*p*p - rho/2*(U&U)));
        //pr == (iter - 1)/iter * pr + 1/iter * ( 1.0/(2.0*rho*sqr(cg))*pa*pa - rho/2.0*(Ua&Ua)); // NOTE: use new form to balance <p_2> and pa.
        //momFlux == (iter - 1)/iter * momFlux + 1/iter * ( rho*(Ua*Ua));
        //runTime.write();
}


// ************************************************************************* //
