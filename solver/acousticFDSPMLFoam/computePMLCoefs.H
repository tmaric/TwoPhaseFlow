{
	const auto& Ccs = mesh.C();

	forAll(Ccs, i)
	{
		auto Omega = 2*Foam::constant::mathematical::pi*(f.value());
		vector DpFunc = vector(0, 0, 0);

		switch (selectedModel)
		{
			case RECTANGLE:
			{
				vector Dist = vector(0, 0, 0);
				for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
				{
					if( (PMLMax.component(cmpt) - Ccs[i].component(cmpt)) < L || (Ccs[i].component(cmpt) - PMLMin.component(cmpt)) < L)
					{
						Dist.component(cmpt) = (((PMLMax.component(cmpt) - Ccs[i].component(cmpt)) < L) ? PMLMax.component(cmpt) - Ccs[i].component(cmpt) : Ccs[i].component(cmpt) - PMLMin.component(cmpt)).value();
						DpFunc.component(cmpt) = (cg/(Dist.component(cmpt) + VSMALL) - cg/L).value(); // inverse distance damping function, see Kaltenbacher2013, p410
						// DpFunc.component(cmpt) = (4*cg/(cg/(f*10))*sqr((L-Dist.component(cmpt))/L)).value();// square distance ratio from Wang2025
					}
				}
				break;
			}
			case CYLINDRICAL:
			{
				dimensionedScalar Dist;
				Dist = mag(Ccs[i] - origin.value());
				if( (Dist - Rmin) < L && (Rmax - Dist) < L)
				{
					DpFunc.x() = (cg*(1/(Rmax - Dist + VSMALL) - 1/L)).value();
					DpFunc.y() = (cg*(1/(Rmax - Dist + VSMALL) - 1/L)).value();
					DpFunc.z() = (cg*(1/(Rmax - Dist + VSMALL) - 1/L)).value();
				}
				break;
			}
			case UNDEFINED:
			{
				FatalErrorInFunction
					<< "Internal error: PML type is undefined."
					<< exit(FatalError);
				break;
			}
		}
		TC0[i].xx()= ((DpFunc.x()*(DpFunc.y() + DpFunc.z() - DpFunc.x()) - DpFunc.y()*DpFunc.z())/(sqr(DpFunc.x())+sqr(Omega)));//.value();
		TC1[i].xx()= ((DpFunc.x()*DpFunc.y()*DpFunc.z()/Omega + Omega*(DpFunc.y() + DpFunc.z() - DpFunc.x()))/(sqr(DpFunc.x())+sqr(Omega)));//.value();

	        TC0[i].yy()= ((DpFunc.y()*(DpFunc.x() + DpFunc.z() - DpFunc.y()) - DpFunc.x()*DpFunc.z())/(sqr(DpFunc.y())+sqr(Omega)));//.value();
	        TC1[i].yy()= ((DpFunc.x()*DpFunc.y()*DpFunc.z()/Omega + Omega*(DpFunc.x() + DpFunc.z() - DpFunc.y()))/(sqr(DpFunc.y())+sqr(Omega)));//.value();

	        TC0[i].zz()= ((DpFunc.z()*(DpFunc.y() + DpFunc.x() - DpFunc.z()) - DpFunc.y()*DpFunc.x())/(sqr(DpFunc.z())+sqr(Omega)));//.value();
	        TC1[i].zz()= ((DpFunc.x()*DpFunc.y()*DpFunc.z()/Omega + Omega*(DpFunc.y() + DpFunc.x() - DpFunc.z()))/(sqr(DpFunc.z())+sqr(Omega)));//.value();


		SC0[i] =((DpFunc.x()*DpFunc.y() + DpFunc.x()*DpFunc.z() + DpFunc.y()*DpFunc.z())/sqr(cg)).value();
		SC1[i] =((Omega*(DpFunc.x() + DpFunc.y() + DpFunc.z()) - (DpFunc.x()*DpFunc.y()*DpFunc.z())/Omega)/sqr(cg)).value();
}

//-----> compute Rhof in boundary field
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
	        auto Omega = 2*Foam::constant::mathematical::pi*f.value();
	        vector DpFunc = vector(0, 0, 0);

		switch (selectedModel)
		{
			case RECTANGLE:
			{
				vector Dist = vector(0, 0, 0);
				for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
				{
					if( (PMLMax.component(cmpt) - Cf.boundaryField()[pi][fi].component(cmpt)) < L || (Cf.boundaryField()[pi][fi].component(cmpt) - PMLMin.component(cmpt)) < L)
					{
						Dist.component(cmpt) = (((PMLMax.component(cmpt) - Cf.boundaryField()[pi][fi].component(cmpt)) < L) ? PMLMax.component(cmpt) - Cf.boundaryField()[pi][fi].component(cmpt) : Cf.boundaryField()[pi][fi].component(cmpt) - PMLMin.component(cmpt)).value();
						DpFunc.component(cmpt) = (cg/(Dist.component(cmpt) + VSMALL) - cg/L).value(); // inverse distance damping function, see Kaltenbacher2013, p410
						// DpFunc.component(cmpt) = (4*cg/(cg/(f*10))*sqr((L-Dist.component(cmpt))/L)).value();
					}
				}
				break;
			}
			case CYLINDRICAL:
			{
				dimensionedScalar Dist;
				Dist = mag(Cf.boundaryField()[pi][fi] - origin.value());
				if( (Dist - Rmin) < L && (Rmax - Dist) < L)
				{
					DpFunc.x() = (cg*(1/(Rmax - Dist + VSMALL) - 1/L)).value();
					DpFunc.y() = (cg*(1/(Rmax - Dist + VSMALL) - 1/L)).value();
					DpFunc.z() = (cg*(1/(Rmax - Dist + VSMALL) - 1/L)).value();
				}
				break;
			}
			case UNDEFINED:
			{
				FatalErrorInFunction
					<< "Internal error: PML type is undefined."
					<< exit(FatalError);
				break;
			}
		}
	                TC0.boundaryFieldRef()[pi][fi].xx()= ((DpFunc.x()*(DpFunc.y() + DpFunc.z() - DpFunc.x()) - DpFunc.y()*DpFunc.z())/(sqr(DpFunc.x())+sqr(Omega)));//.value();
	                TC1.boundaryFieldRef()[pi][fi].xx()= ((DpFunc.x()*DpFunc.y()*DpFunc.z()/Omega + Omega*(DpFunc.y() + DpFunc.z() - DpFunc.x()))/(sqr(DpFunc.x())+sqr(Omega)));//.value();

	                TC0.boundaryFieldRef()[pi][fi].yy()= ((DpFunc.y()*(DpFunc.x() + DpFunc.z() - DpFunc.y()) - DpFunc.x()*DpFunc.z())/(sqr(DpFunc.y())+sqr(Omega)));//.value();
	                TC1.boundaryFieldRef()[pi][fi].yy()= ((DpFunc.x()*DpFunc.y()*DpFunc.z()/Omega + Omega*(DpFunc.x() + DpFunc.z() - DpFunc.y()))/(sqr(DpFunc.y())+sqr(Omega)));//.value();

	                TC0.boundaryFieldRef()[pi][fi].zz()= ((DpFunc.z()*(DpFunc.y() + DpFunc.x() - DpFunc.z()) - DpFunc.y()*DpFunc.x())/(sqr(DpFunc.z())+sqr(Omega)));//.value();
	                TC1.boundaryFieldRef()[pi][fi].zz()= ((DpFunc.x()*DpFunc.y()*DpFunc.z()/Omega + Omega*(DpFunc.y() + DpFunc.x() - DpFunc.z()))/(sqr(DpFunc.z())+sqr(Omega)));//.value();

	                SC0.boundaryFieldRef()[pi][fi] =((DpFunc.x()*DpFunc.y() + DpFunc.x()*DpFunc.z() + DpFunc.y()*DpFunc.z())/sqr(cg)).value();
	                SC1.boundaryFieldRef()[pi][fi] =((Omega*(DpFunc.x() + DpFunc.y() + DpFunc.z()) - (DpFunc.x()*DpFunc.y()*DpFunc.z())/Omega)/sqr(cg)).value();	
}
}
}
}
