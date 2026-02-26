{
    #include "computeAlphaf.H"

        while (pimple.correct())
        {
            #include "paEqn.H"
        }

	pr == 1.0/(2.0*rho*sqr(cg))*pa*pa - rho/2.0*(Ua&Ua);
	momFlux == rho*(Ua*Ua);
}


// ************************************************************************* //
