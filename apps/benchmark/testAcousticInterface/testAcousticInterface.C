#include "fvCFD.H"
#include "transmissionFit.H"
#include <cmath>
using namespace Foam::acoustic;
int main(int argc,char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    volScalarField alpha(IOobject("testAlpha",runTime.timeName(),mesh),mesh,
        dimensionedScalar(dimless,Zero));
    cutFaceAdvect cutter(mesh,alpha);
    label checks=0;
    scalar worst=0;
    auto check=[&](bool ok,const char* message)
    {
        if(!ok) FatalErrorInFunction << message << " check=" << checks << exit(FatalError);
        ++checks;
    };
    pointField square({point(0,-1,-1),point(0,1,-1),point(0,1,1),point(0,-1,1)});
    face poly({0,1,2,3});
    for(scalar offset:{-1.0,-0.9,0.0,0.3,0.999,1.0})
    {
        auto a=cutArea(poly,square,vector(4,0,0),point(0,offset,0),vector(0,1,0),cutter);
        check(mag(a.fraction-(1-offset)/2)<1e-12,"square cut");
        auto b=cutArea(poly.reverseFace(),square,vector(-4,0,0),
            point(0,offset,0),vector(0,-1,0),cutter);
        check(mag(a.fraction+b.fraction-1)<1e-12,"phase complement");
        check(mag(a.liquid+a.gas-vector(4,0,0))<1e-12,"area partition");
    }
    pointField warped(square); warped[2].x()=0.3;
    vector ws=poly.areaNormal(warped);
    auto wa=cutArea(poly,warped,ws,point(0,0.1,0),vector(0,1,0),cutter);
    auto wb=cutArea(poly.reverseFace(),warped,-ws,point(0,0.1,0),vector(0,1,0),cutter);
    check(mag(wa.fraction-wb.fraction)<1e-12 && mag(wa.liquid+wb.liquid)<1e-12,"warped reversal");

    for(label dimension:{1,2,3}) for(scalar ratio:{1.0,10.0,1000.0,1e6})
    for(label angle=0;angle<13;++angle) for(scalar offset:{-0.4,-0.1,0.0,0.1,0.4})
    {
        vectorField basis(dimension,vector::zero);
        forAll(basis,j) basis[j][j]=1;
        scalar theta=dimension>1 ? constant::mathematical::pi*angle/12:0;
        vector n(Foam::cos(theta),Foam::sin(theta),0);
        if(dimension==3)
        {
            n.z()=0.31*Foam::sin(2*theta+0.3); n/=mag(n);
        }
        DynamicList<point> cloud;
        cloud.append(point(-0.5,0,0)); cloud.append(point(0.5,0,0));
        cloud.append(point(-1.5,0,0)); cloud.append(point(1.5,0,0));
        if(dimension>=2) for(scalar x:{-0.5,0.5}) for(scalar y:{-1.,1.})
            cloud.append(point(x,y,0));
        if(dimension==3) for(scalar x:{-0.5,0.5}) for(scalar z:{-1.,1.})
            cloud.append(point(x,0,z));
        pointField pts(cloud);
        if(dimension>1 && angle%3)
        {
            pts[0].y()=-0.13; pts[1].y()=0.21;
            pts[2].x()-=0.15; pts[3].x()+=0.22;
            if(dimension==3) { pts[0].z()=0.09; pts[1].z()=-0.11; }
        }
        point xc=offset*n;
        auto area=cutArea(poly,square,vector(4,0,0),xc,n,cutter);
        auto fit=transmissionFit(pts,point::zero,xc,n,area.liquid,area.gas,
            basis,ratio,1,1e-12,1e10);
        check(fit.valid,"fit rank/conditioning");
        vector t(0.31,-0.47,dimension==3 ? 0.27:0);
        for(direction d=dimension;d<3;++d) t[d]=0;
        t-=(t&n)*n;
        scalar q=0.37/ratio;
        scalar calculated=0,constant=0;
        forAll(pts,k)
        {
            vector r=pts[k]-xc;
            scalar pressure=0.23+(t&r)+(r&n)*((r&n)>=0 ? ratio:1)*q;
            calculated+=fit.weights[k]*pressure;
            constant+=fit.weights[k];
        }
        scalar exact=((t/ratio+q*n)&area.liquid)+((t+q*n)&area.gas);
        scalar error=mag(calculated-exact);
        worst=max(worst,error);
        check(error<1e-10,"piecewise affine transmission");
        check(mag(constant)<1e-10,"constant preservation");
        if(ratio==1 && (dimension==1 || angle%3==0))
            check(mag(fit.weights[0]+4/ratio)<1e-10
               && mag(fit.weights[1]-4/ratio)<1e-10,"two point limit");
        std::swap(pts[0],pts[1]);
        auto reverse=transmissionFit(pts,point::zero,xc,n,-area.liquid,-area.gas,
            basis,ratio,1,1e-12,1e10);
        check(reverse.valid,"reverse fit");
        std::swap(reverse.weights[0],reverse.weights[1]);
        check(gMax(mag(reverse.weights+fit.weights))<1e-10,"owner reversal");
    }
    // Orthogonal coincident interface, unequal distances and densities.
    pointField pts({point(-0.25,0,0),point(0.75,0,0)});
    vectorField basis({vector(1,0,0)});
    auto f=transmissionFit(pts,point::zero,point::zero,vector(1,0,0),
        vector(1,0,0),vector::zero,basis,1000,1.2,1e-12,1e10);
    check(f.valid && mag(f.weights[1]-1/(0.25*1.2+0.75*1000))<1e-12,
        "unequal distance harmonic limit");
    // Degenerate stencil must fail, never silently use an arbitrary inverse.
    auto fail=transmissionFit(pts,point::zero,point::zero,vector(1,0,0),
        vector(1,0,0),vector::zero,vectorField({vector(1,0,0),vector(0,1,0)}),
        1000,1.2,1e-12,1e10);
    check(!fail.valid,"degenerate fit rejection");
    Info<< "ACOUSTIC_INTERFACE_TESTS_PASSED checks=" << checks
        << " maxAffineError=" << worst << nl;
    return 0;
}
