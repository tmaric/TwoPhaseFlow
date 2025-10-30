{
	const auto& Ccs = mesh.C();

	forAll(Ccs, i)
	{
		//vector Dist=vector(0, 0, 0);
		vector DpFunc=vector(0, 0, 0);
		//dimensionedVector Ccsi(dimLength, Ccs[i]);
switch (selectedModel)
{
    case RECTANGLE:
    {
		vector Dist=vector(0, 0, 0);
                for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
                {
                        if( (PMLMax.component(cmpt) - Ccs[i].component(cmpt)) < L || (Ccs[i].component(cmpt) - PMLMin.component(cmpt)) < L)
                        {
                                Dist.component(cmpt) = ((PMLMax.component(cmpt) - Ccs[i].component(cmpt) < L) ? (PMLMax.component(cmpt) - Ccs[i].component(cmpt) ).value(): (Ccs[i].component(cmpt) - PMLMin.component(cmpt)).value());
                                DpFunc.component(cmpt) = (cg*(1/(Dist.component(cmpt)+VSMALL) - 1/L)).value(); // inverse distance damping function, see Kaltenbacher2013, p410
                                // DpFunc.component(cmpt) = (4*cg/(cg/(f*10))*sqr((L-Dist.component(cmpt))/L)).value();// square distance ratio from Wang2025
                        }
                }

        break; // <-- Don't forget the 'break'!
    }

    case CYLINDRICAL:
    {
		dimensionedScalar Dist;
		Dist = mag(Ccs[i] - origin.value());

                if( (Dist - Rmin) < L && (Rmax - Dist) < L)
                        {
                                DpFunc.x() = (cg*(1/(Rmax - Dist +VSMALL) - 1/L)).value(); 
				DpFunc.y() = (cg*(1/(Rmax - Dist +VSMALL) - 1/L)).value(); 
				DpFunc.z() = (cg*(1/(Rmax - Dist +VSMALL) - 1/L)).value(); // inverse distance damping function, see Kaltenbacher2013, p410
                                // DpFunc.component(cmpt) = (4*cg/(cg/(f*10))*sqr((L-Dist.component(cmpt))/L)).value();// square distance ratio from Wang2025
                        }

        break;
    }
    case UNDEFINED:
    {
        // This case should not be reachable because of the 
        // FatalError, but it's good practice for completeness.
        FatalErrorInFunction
            << "Internal error: PML type is undefined."
            << exit(FatalError);
        break;
    }
}		

		TC0[i].xx()=DpFunc.x();
		TC0[i].yy()=DpFunc.y();
		TC0[i].zz()=DpFunc.z();

		TC1[i].xx()=DpFunc.x() - DpFunc.y() - DpFunc.z();
		TC1[i].yy()=DpFunc.y() - DpFunc.x() - DpFunc.z();
		TC1[i].zz()=DpFunc.z() - DpFunc.x() - DpFunc.y();

		TC2[i].xx()=DpFunc.y() * DpFunc.z();
		TC2[i].yy()=DpFunc.x() * DpFunc.z();
		TC2[i].zz()=DpFunc.x() * DpFunc.y();

		SC0[i] =(DpFunc.x() + DpFunc.y() + DpFunc.z())/(sqr(cg)).value(); 
		SC1[i] =(DpFunc.x()*DpFunc.y() + DpFunc.x()*DpFunc.z() + DpFunc.y()*DpFunc.z())/(sqr(cg)).value(); 
		SC2[i] =(DpFunc.x()*DpFunc.y()*DpFunc.z())/(sqr(cg)).value();
	}
	
	//update the boundary fields
	const polyBoundaryMesh& patches = mesh.boundaryMesh();
	const surfaceVectorField& Cf = mesh.Cf();

	forAll(patches, pi)
	{
		const polyPatch& patchi= patches[pi];
		Info << "patch name: " << patchi.name() << endl;
		if
		(
			!isA<emptyPolyPatch>(patchi) && !isA<wedgePolyPatch>(patchi)
		)
		{
			forAll(patchi, fi)
			{
//				vector Dist=vector(0, 0, 0);
				vector DpFunc=vector(0, 0, 0);

switch (selectedModel)
{
    case RECTANGLE:
    {

                vector Dist=vector(0, 0, 0);

                for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
                {
                        if( (PMLMax.component(cmpt) - Cf.boundaryField()[pi][fi].component(cmpt)) < L || (Cf.boundaryField()[pi][fi].component(cmpt) - PMLMin.component(cmpt)) < L)
                        {
                                Dist.component(cmpt) = ((PMLMax.component(cmpt) - Cf.boundaryField()[pi][fi].component(cmpt) < L) ? (PMLMax.component(cmpt) - Cf.boundaryField()[pi][fi].component(cmpt)).value() : (Cf.boundaryField()[pi][fi].component(cmpt) - PMLMin.component(cmpt))).value(); 
                                DpFunc.component(cmpt) = (cg*(1/(Dist.component(cmpt) + VSMALL) - 1/L)).value(); //- cg/(L*Omega); // inverse distance damping function, see Kaltenbache>
        // DpFunc.component(cmpt) = (4*cg/(cg/(f*10))*sqr((L-Dist.component(cmpt))/L)).value();
                        }
                }

        break; // <-- Don't forget the 'break'!
    }

    case CYLINDRICAL:
    {
                dimensionedScalar Dist;
                Dist = mag(Cf.boundaryField()[pi][fi] - origin.value());

                if( (Dist - Rmin) < L && (Rmax - Dist) < L)
                        {
                                DpFunc.x() = (cg*(1/(Rmax - Dist +VSMALL) - 1/L)).value(); 
                                DpFunc.y() = (cg*(1/(Rmax - Dist +VSMALL) - 1/L)).value(); 
                                DpFunc.z() = (cg*(1/(Rmax - Dist +VSMALL) - 1/L)).value(); // inverse distance damping function, see Kaltenbacher2013, p410
                                // DpFunc.component(cmpt) = (4*cg/(cg/(f*10))*sqr((L-Dist.component(cmpt))/L)).value();// square distance ratio from Wang2025
                        }

        break;
    }
    case UNDEFINED:
    {
        // This case should not be reachable because of the 
        // FatalError, but it's good practice for completeness.
        FatalErrorInFunction
            << "Internal error: PML type is undefined."
            << exit(FatalError);
        break;
    }
}
				
				TC0.boundaryFieldRef()[pi][fi].xx()=DpFunc.x();
				TC0.boundaryFieldRef()[pi][fi].yy()=DpFunc.y();
				TC0.boundaryFieldRef()[pi][fi].zz()=DpFunc.z();

				TC1.boundaryFieldRef()[pi][fi].xx()=DpFunc.x() - DpFunc.y() - DpFunc.z();
				TC1.boundaryFieldRef()[pi][fi].yy()=DpFunc.y() - DpFunc.x() - DpFunc.z();
				TC1.boundaryFieldRef()[pi][fi].zz()=DpFunc.z() - DpFunc.x() - DpFunc.y();

				TC2.boundaryFieldRef()[pi][fi].xx()=DpFunc.y() * DpFunc.z();
				TC2.boundaryFieldRef()[pi][fi].yy()=DpFunc.x() * DpFunc.z();
				TC2.boundaryFieldRef()[pi][fi].zz()=DpFunc.x() * DpFunc.y();

				SC0.boundaryFieldRef()[pi][fi] =(DpFunc.x() + DpFunc.y() + DpFunc.z())/(sqr(cg)).value(); 
				SC1.boundaryFieldRef()[pi][fi] =(DpFunc.x()*DpFunc.y() + DpFunc.x()*DpFunc.z() + DpFunc.y()*DpFunc.z())/(sqr(cg)).value(); 
				SC2.boundaryFieldRef()[pi][fi] =(DpFunc.x()*DpFunc.y()*DpFunc.z())/(sqr(cg)).value();
			}
		}
	}
}
