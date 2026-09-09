#include "transmissionFit.H"
#include "SVD.H"
#include <algorithm>
#include <cmath>
#include <vector>

namespace Foam { namespace acoustic
{
namespace
{
bool lessPoint(const point& a, const point& b)
{
    for (direction d=0; d<3; ++d)
    {
        if (a[d] < b[d]) return true;
        if (a[d] > b[d]) return false;
    }
    return false;
}
bool goodSvd(const SVD& s, label rank, scalar tol, scalar limit, scalar& cond)
{
    if (!s.converged()) return false;
    std::vector<scalar> values;
    scalar largest = 0;
    forAll(s.S(), i) largest = max(largest, mag(s.S()[i]));
    if (!(largest > VSMALL) || !std::isfinite(largest)) return false;
    forAll(s.S(), i)
        if (mag(s.S()[i]) > tol*largest) values.push_back(mag(s.S()[i]));
    if (label(values.size()) != rank) return false;
    cond = largest / *std::min_element(values.begin(), values.end());
    return std::isfinite(cond) && cond <= limit;
}
}

CutArea cutArea
(
    const face& f, const pointField& allPoints, const vector& Sf,
    const point& centre, const vector& normal, cutFaceAdvect& cutter
)
{
    if (f.size() < 3 || mag(Sf) <= VSMALL)
        FatalErrorInFunction << "Degenerate acoustic face" << exit(FatalError);
    // Canonical starting vertex and traversal remove owner/proc ordering.
    label first = 0;
    forAll(f, i) if (lessPoint(allPoints[f[i]], allPoints[f[first]])) first=i;
    const label step = lessPoint
    (
        allPoints[f[(first+1)%f.size()]],
        allPoints[f[(first+f.size()-1)%f.size()]]
    ) ? 1 : -1;
    pointField pts(f.size());
    forAll(pts, i) pts[i]=allPoints[f[(first+step*i+f.size())%f.size()]];
    const scalar h = sqrt(mag(Sf));
    scalarField distances(pts.size());
    forAll(pts, i)
    {
        distances[i] = ((pts[i]-centre)&normal)/h;
        if (mag(distances[i]) < 1e-12) distances[i]=0;
    }

    // Choose the first valid fan in canonical order. Concave faces that have
    // no valid vertex fan are rejected, rather than integrating overlaps.
    label base=-1;
    forAll(pts, b)
    {
        bool ok=true;
        scalar orientation=0;
        for (label i=1; i<pts.size()-1; ++i)
        {
            vector a=0.5*((pts[(b+i)%pts.size()]-pts[b])
                         ^(pts[(b+i+1)%pts.size()]-pts[b]));
            scalar dot=a&Sf;
            if (mag(a) <= 1e-14*mag(Sf)) continue;
            if (mag(dot) <= 1e-14*mag(a)*mag(Sf)) { ok=false; break; }
            if (!orientation) orientation=dot;
            if (dot*orientation < 0) { ok=false; break; }
        }
        if (ok && orientation) { base=b; break; }
    }
    if (base < 0)
        FatalErrorInFunction << "No non-overlapping acoustic face fan" << exit(FatalError);
    scalar total=0, liquid=0;
    vector areaSum=vector::zero, liquidVector=vector::zero;
    for (label i=1; i<pts.size()-1; ++i)
    {
        face tri(3);
        tri[0]=base; tri[1]=(base+i)%pts.size(); tri[2]=(base+i+1)%pts.size();
        vector a=tri.areaNormal(pts);
        if (mag(a) <= 1e-14*mag(Sf)) continue;
        if ((a&Sf)<0) { std::swap(tri[1],tri[2]); a=-a; }
        cutter.calcSubFace(tri,pts,distances,0);
        vector l=cutter.subFaceArea();
        if ((l&a)<0) l=-l;
        total+=mag(a); liquid+=mag(l);
        areaSum+=a; liquidVector+=l;
    }
    if (total<=VSMALL || mag(areaSum-Sf)>1e-10*mag(Sf))
        FatalErrorInFunction << "Acoustic triangulation area mismatch" << exit(FatalError);
    CutArea out;
    out.fraction=liquid/total;
    if (!std::isfinite(out.fraction) || out.fraction < -1e-12 || out.fraction > 1+1e-12)
        FatalErrorInFunction << "Invalid acoustic liquid area fraction" << exit(FatalError);
    out.fraction=min(max(out.fraction,scalar(0)),scalar(1));
    out.liquid=liquidVector;
    out.gas=Sf-liquidVector;
    return out;
}

FitResult transmissionFit
(
    const pointField& pts, const point& cf, const point& xc, const vector& nInput,
    const vector& sl, const vector& sg, const vectorField& basis,
    scalar rhoL, scalar rhoG, scalar tol, scalar limit
)
{
    FitResult out;
    const label dim=basis.size(), m=dim+1, N=pts.size();
    if (dim<1 || dim>3 || N<m || rhoL<=0 || rhoG<=0) return out;
    vector n=vector::zero;
    forAll(basis,i) n+=(basis[i]&nInput)*basis[i];
    if (mag(n)<1e-12) return out;
    n/=mag(n);
    vectorField tangent(dim-1);
    label nt=0;
    forAll(basis,i)
    {
        vector t=basis[i]-(basis[i]&n)*n;
        for (label j=0;j<nt;++j) t-=(t&tangent[j])*tangent[j];
        if (mag(t)>1e-10 && nt<dim-1) tangent[nt++]=t/mag(t);
    }
    if (nt!=dim-1) return out;
    scalar h=mag(pts[1]-pts[0]), ref=max(rhoL,rhoG);
    if (h<=VSMALL) return out;
    scalarRectangularMatrix B(N,m,Zero);
    forAll(pts,i)
    {
        vector r=(pts[i]-xc)/h;
        scalar s=r&n;
        scalar rho=s>=0 ? rhoL:rhoG;
        B(i,0)=1;
        forAll(tangent,j) B(i,j+1)=tangent[j]&r;
        B(i,m-1)=(rho/ref)*s;
    }
    // Padding gives a square constraint matrix, including its full nullspace,
    // without requiring an SVD implementation supporting wide matrices.
    scalarRectangularMatrix C(m,m,Zero);
    for (label i=0;i<2;++i) for(label j=0;j<m;++j) C(i,j)=B(i,j);
    SVD cs(C,tol);
    scalar cc;
    if (!goodSvd(cs,2,tol,limit,cc)) return out;
    const scalarRectangularMatrix cp=cs.VSinvUt();
    scalarRectangularMatrix X(m,N,Zero);
    for(label j=0;j<m;++j) { X(j,0)=cp(j,0); X(j,1)=cp(j,1); }
    scalarRectangularMatrix Z(m,m-2,Zero);
    label nz=0;
    scalar smax=0;
    forAll(cs.S(),i) smax=max(smax,mag(cs.S()[i]));
    forAll(cs.S(),k) if (mag(cs.S()[k])<=tol*smax)
    {
        for(label j=0;j<m;++j) Z(j,nz)=cs.V()(j,k);
        ++nz;
    }
    scalar dc=1;
    if (m>2)
    {
        scalarRectangularMatrix D(N-2,m-2,Zero), R(N-2,N,Zero);
        for(label i=2;i<N;++i)
        {
            scalar wi=1/max(mag(pts[i]-cf)/h,scalar(0.25));
            for(label j=0;j<m-2;++j)
                for(label k=0;k<m;++k) D(i-2,j)+=wi*B(i,k)*Z(k,j);
            R(i-2,i)=wi;
            for(label j=0;j<N;++j)
                for(label k=0;k<m;++k) R(i-2,j)-=wi*B(i,k)*X(k,j);
        }
        SVD ds(D,tol);
        if (!goodSvd(ds,m-2,tol,limit,dc)) return out;
        scalarRectangularMatrix pinv=ds.VSinvUt();
        for(label j=0;j<m;++j) for(label k=0;k<N;++k)
            for(label a=0;a<m-2;++a) for(label b=0;b<N-2;++b)
                X(j,k)+=Z(j,a)*pinv(a,b)*R(b,k);
    }
    scalarField f(m,Zero);
    forAll(tangent,j) f[j+1]=((sl/rhoL+sg/rhoG)&tangent[j])/h;
    f[m-1]=((sl+sg)&n)/(ref*h);
    out.weights.setSize(N,Zero);
    for(label k=0;k<N;++k)
        for(label j=0;j<m;++j) out.weights[k]+=f[j]*X(j,k);
    scalar defect=0, scale=0;
    forAll(out.weights,k)
    {
        if (!std::isfinite(out.weights[k])) return out;
        defect+=out.weights[k]; scale+=mag(out.weights[k]);
    }
    if (mag(defect)>1e-10*max(scale,VSMALL)) return out;
    out.condition=max(cc,dc);
    out.valid=true;
    return out;
}
}}
