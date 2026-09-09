#include "acousticInterface.H"
#include "reconstructionSchemes.H"
#include "upwind.H"
#include "syncTools.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "PstreamBuffers.H"
#include "UOPstream.H"
#include "UIPstream.H"
#include "OFstream.H"
#include <set>
#include <cmath>
#include <algorithm>

namespace Foam { namespace acoustic
{
namespace
{
#include "legacy/processorBC.H"
void legacyAreas
(
    const fvMesh& mesh,volScalarField& alpha1,surfaceScalarField& alphaf,
    const surfaceScalarField& phi,cutFaceAdvect& cutFace
)
{
#include "legacy/computeAlphaf.H"
}
bool symmetryPatch(const polyPatch& p)
{
    return isA<emptyPolyPatch>(p) || isA<wedgePolyPatch>(p)
        || isA<symmetryPolyPatch>(p) || isA<symmetryPlanePolyPatch>(p);
}
template<class T> List<T> asList(const std::set<T>& input)
{
    List<T> out(input.size());
    label i=0; for(const T& x:input) out[i++]=x;
    return out;
}
bool finiteVector(const vector& v)
{
    return std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z());
}
}

acousticInterface::acousticInterface(const fvMesh& mesh,const globalIndex& global)
:
    mesh_(mesh),global_(global),areaName_("legacy"),fluxName_("legacy"),
    diagnostics_(false),maxRings_(3),svdTol_(1e-12),maxCondition_(1e10)
{
    const dictionary& schemes=mesh.schemesDict();
    if (schemes.found("acousticInterface"))
    {
        const dictionary& d=schemes.subDict("acousticInterface");
        areaName_=d.lookupOrDefault<word>("areaFraction","legacy");
        fluxName_=d.lookupOrDefault<word>("flux","legacy");
        diagnostics_=d.lookupOrDefault<Switch>("writeDiagnostics",false);
        if (d.found("plicTransmissionCoeffs"))
        {
            const dictionary& c=d.subDict("plicTransmissionCoeffs");
            maxRings_=c.lookupOrDefault<label>("maxStencilRings",3);
            svdTol_=c.lookupOrDefault<scalar>("svdRelativeTolerance",1e-12);
            maxCondition_=c.lookupOrDefault<scalar>("maxConditionNumber",1e10);
        }
    }
    if ((areaName_!="legacy" && areaName_!="plicAverage")
     || (fluxName_!="legacy" && fluxName_!="plicTransmission"))
        FatalErrorInFunction << "Unknown acousticInterface method. areaFraction: "
            << "legacy or plicAverage; flux: legacy or plicTransmission." << exit(FatalError);
    if (transmission() && areaName_!="plicAverage")
        FatalErrorInFunction << "plicTransmission requires areaFraction plicAverage."
            << " Supported combinations: legacy/legacy, plicAverage/legacy,"
            << " plicAverage/plicTransmission." << exit(FatalError);
    if (maxRings_<1 || !(svdTol_>0 && svdTol_<1)
        || !(maxCondition_>=1 && std::isfinite(maxCondition_)))
        FatalErrorInFunction << "Invalid plicTransmissionCoeffs" << exit(FatalError);
    Info<< "acousticInterface: areaFraction=" << areaName_ << " flux=" << fluxName_
        << " writeDiagnostics=" << diagnostics_ << " maxStencilRings=" << maxRings_
        << " svdRelativeTolerance=" << svdTol_
        << " maxConditionNumber=" << maxCondition_ << nl;
    if (areaName_!="legacy")
    {
        forAll(mesh.boundaryMesh(),p)
            // processorCyclic also derives from processorPolyPatch. Only
            // ordinary processor patches have the face pairing used below.
            if ((mesh.boundaryMesh()[p].coupled()
              && mesh.boundaryMesh()[p].type()!=processorPolyPatch::typeName)
             || mesh.boundaryMesh()[p].type()=="overset")
                FatalErrorInFunction << "New acoustic interface methods do not support patch "
                    << mesh.boundaryMesh()[p].name() << " of type "
                    << mesh.boundaryMesh()[p].type() << exit(FatalError);
    }

    // Empty dimensions are excluded. Wedges additionally constrain gradients
    // to their centre plane, even though geometricD includes all 3 dimensions.
    vector normal=vector::zero;
    forAll(mesh.boundaryMesh(),p)
        if (isA<wedgePolyPatch>(mesh.boundaryMesh()[p]))
        {
            normal=refCast<const wedgePolyPatch>(mesh.boundaryMesh()[p]).centreNormal();
            break;
        }
    DynamicList<vector> active;
    for(direction d=0;d<3;++d) if (mesh.geometricD()[d]>0)
    {
        vector e=vector::zero; e[d]=1;
        if (mag(normal)>SMALL) e-=(e&normal)*normal;
        forAll(active,j) e-=(e&active[j])*active[j];
        if(mag(e)>1e-10) active.append(e/mag(e));
    }
    basis_=vectorField(active);
}

void acousticInterface::computeAreas
(
    volScalarField& alpha,surfaceScalarField& af,
    const surfaceScalarField& phi,cutFaceAdvect& cutter
)
{
    if(areaName_=="legacy") legacyAreas(mesh_,alpha,af,phi,cutter);
    else computeGeometricAreas(alpha,af,cutter);
}

void acousticInterface::computeGeometricAreas
(
    volScalarField& alpha,surfaceScalarField& af,cutFaceAdvect& cutter
)
{
    reconstructionSchemes& surf=
        mesh_.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");
    surf.reconstruct(true);
    pureTol_=surf.modelDict().lookupOrDefault<scalar>("surfCellTol",1e-8);
    if (!(pureTol_>0 && pureTol_<0.5))
        FatalErrorInFunction << "Invalid reconstruction surfCellTol" << exit(FatalError);
    cells_.clear(); faces_.clear();
    const label nb=mesh_.nBoundaryFaces(), ni=mesh_.nInternalFaces();
    vectorField nc(nb,vector::zero), pc(nb,vector::zero), cc(nb,vector::zero);
    scalarField ac(nb,Zero);
    neighbourGlobal_.setSize(nb,-1);
    label invalid=0;
    forAll(alpha,i)
    {
        CellData c;
        c.centre=mesh_.C()[i]; c.alpha=alpha[i];
        c.normal=surf.normal()[i]; c.planeCentre=surf.centre()[i];
        if (!std::isfinite(c.alpha) || c.alpha < -1e-12 || c.alpha>1+1e-12
            || !finiteVector(c.normal) || !finiteVector(c.planeCentre))
        {
            Pout<< "Invalid PLIC geometry in global cell " << global_.toGlobal(i)
                << " local cell " << i << " faces " << mesh_.cells()[i]
                << " alpha=" << c.alpha << " plane centre " << c.planeCentre
                << " normal " << c.normal << nl;
            ++invalid;
        }
        if(mag(c.normal)>VSMALL) c.normal/=mag(c.normal);
        else if(c.alpha>pureTol_ && c.alpha<1-pureTol_)
        {
            Pout<< "Missing PLIC plane in local cell " << i << " centre "
                << c.centre << " alpha=" << c.alpha
                << " global cell " << global_.toGlobal(i)
                << " faces " << mesh_.cells()[i] << " surfCellTol=" << pureTol_ << nl;
            ++invalid;
        }
        cells_.emplace(global_.toGlobal(i),c);
    }
    if(returnReduce(invalid,sumOp<label>()))
        FatalErrorInFunction << "Invalid or missing acoustic PLIC data" << exit(FatalError);

    for(label b=0;b<nb;++b)
    {
        label own=mesh_.faceOwner()[ni+b];
        const CellData& c=cells_.at(global_.toGlobal(own));
        nc[b]=c.normal; pc[b]=c.planeCentre; cc[b]=c.centre;
        ac[b]=c.alpha; neighbourGlobal_[b]=global_.toGlobal(own);
    }
    syncTools::swapBoundaryFaceList(mesh_,nc);
    syncTools::swapBoundaryFaceList(mesh_,pc);
    syncTools::swapBoundaryFaceList(mesh_,cc);
    syncTools::swapBoundaryFaceList(mesh_,ac);
    syncTools::swapBoundaryFaceList(mesh_,neighbourGlobal_);
    const surfaceScalarField linear(fvc::interpolate(alpha));
    af==linear;
    label coincident=0, two=0, one=0;
    for(label f=0;f<mesh_.nFaces();++f)
    {
        const label own=mesh_.faceOwner()[f];
        const CellData& a=cells_.at(global_.toGlobal(own));
        CellData b;
        label nei=-1, ng=-1, patch=-1;
        if(f<ni)
        {
            nei=mesh_.faceNeighbour()[f]; ng=global_.toGlobal(nei); b=cells_.at(ng);
        }
        else
        {
            patch=mesh_.boundaryMesh().whichPatch(f);
            if(symmetryPatch(mesh_.boundaryMesh()[patch])) continue;
            if(isA<processorPolyPatch>(mesh_.boundaryMesh()[patch]))
            {
                ng=neighbourGlobal_[f-ni]; b.normal=nc[f-ni]; b.planeCentre=pc[f-ni];
                b.centre=cc[f-ni]; b.alpha=ac[f-ni];
            }
        }
        vector sf=mesh_.faceAreas()[f];
        FaceData fd; fd.face=f; fd.owner=own; fd.neighbour=nei; fd.neighbourGlobal=ng;
        for(const CellData* c : std::initializer_list<const CellData*>{&a,&b})
        {
            if(c==&b && ng<0) continue;
            if(mag(c->normal)>VSMALL)
            {
                Candidate ca;
                ca.centre=c->planeCentre; ca.normal=c->normal;
                ca.area=cutArea(mesh_.faces()[f],mesh_.points(),sf,ca.centre,ca.normal,cutter);
                fd.candidates.push_back(ca);
            }
        }
        scalar value=0;
        if(!fd.candidates.empty())
        {
            for(const auto& ca:fd.candidates) value+=ca.area.fraction;
            value/=fd.candidates.size();
            if(fd.candidates.size()==2) ++two; else ++one;
        }
        else if(ng>=0 && ((a.alpha<=pureTol_ && b.alpha>=1-pureTol_)
                        || (b.alpha<=pureTol_ && a.alpha>=1-pureTol_)))
        {
            // On a coincident interface scalar face coverage is ambiguous.
            // Keep the baseline scalar fallback, but use a sharp flux plane.
            value=f<ni ? linear[f] :
                linear.boundaryField()[patch][f-mesh_.boundaryMesh()[patch].start()];
            Candidate ca;
            ca.centre=mesh_.faceCentres()[f];
            ca.normal=sf/mag(sf)*(b.alpha>a.alpha ? 1:-1);
            ca.area.fraction=value; ca.area.liquid=value*sf; ca.area.gas=sf-ca.area.liquid;
            fd.candidates.push_back(ca);
            ++coincident;
        }
        else value=a.alpha>=1-pureTol_ ? 1:0;
        if(f<ni) af[f]=value;
        else af.boundaryFieldRef()[patch][f-mesh_.boundaryMesh()[patch].start()]=value;
        if(!fd.candidates.empty()) faces_.emplace(f,std::move(fd));
    }
    // Both neighbouring planes were exchanged; both copies now compute the
    // same unsigned average. Assert it, never overwrite a disagreement.
    scalarField shared(nb,Zero);
    forAll(mesh_.boundary(),p) forAll(af.boundaryField()[p],i)
        shared[mesh_.boundaryMesh()[p].start()+i-ni]=af.boundaryField()[p][i];
    scalarField received(shared);
    syncTools::swapBoundaryFaceList(mesh_,received);
    scalar mismatch=0;
    forAll(mesh_.boundary(),p) if(isA<processorPolyPatch>(mesh_.boundaryMesh()[p]))
        forAll(af.boundaryField()[p],i)
        {
            label k=mesh_.boundaryMesh()[p].start()+i-ni;
            mismatch=max(mismatch,mag(shared[k]-received[k]));
        }
    reduce(mismatch,maxOp<scalar>());
    if(mismatch>1e-12)
        FatalErrorInFunction << "PLIC processor area disagreement " << mismatch << exit(FatalError);
    Info<< "PLIC areas: onePlane=" << returnReduce(one,sumOp<label>())
        << " twoPlanes=" << returnReduce(two,sumOp<label>())
        << " coincident=" << returnReduce(coincident,sumOp<label>())
        << " processorMismatch=" << mismatch
        << " (processor faces counted on both ranks)" << nl;
}

void acousticInterface::exchangeCells(const labelList& needed)
{
    if(!Pstream::parRun()) return;
    List<DynamicList<label>> requests(Pstream::nProcs());
    for(const label id:needed)
        if(cells_.find(id)==cells_.end())
            requests[global_.whichProcID(id)].append(id);
    PstreamBuffers ask(Pstream::commsTypes::nonBlocking);
    forAll(requests,p) if(p!=Pstream::myProcNo() && !requests[p].empty())
        { UOPstream stream(p,ask); stream<<requests[p]; }
    ask.finishedSends();
    PstreamBuffers answer(Pstream::commsTypes::nonBlocking);
    forAll(requests,p) if(p!=Pstream::myProcNo() && ask.recvDataCount(p))
    {
        labelList ids; UIPstream request(p,ask); request>>ids;
        UOPstream os(p,answer);
        os<<ids.size();
        for(const label id:ids)
        {
            if(!global_.isLocal(id))
                FatalErrorInFunction << "Halo request routed to wrong rank" << exit(FatalError);
            const CellData& c=cells_.at(id);
            os<<id<<c.centre<<c.planeCentre<<c.normal<<c.alpha<<c.pml<<c.adjacent;
        }
    }
    answer.finishedSends();
    forAll(requests,p) if(p!=Pstream::myProcNo() && answer.recvDataCount(p))
    {
        UIPstream is(p,answer); label count; is>>count;
        for(label i=0;i<count;++i)
        {
            label id; CellData c;
            is>>id>>c.centre>>c.planeCentre>>c.normal>>c.alpha>>c.pml>>c.adjacent;
            cells_[id]=std::move(c);
        }
    }
    for(const label id:needed) if(cells_.find(id)==cells_.end())
        FatalErrorInFunction << "Missing acoustic halo cell " << id << exit(FatalError);
}

void acousticInterface::build(const volTensorField& sigma,scalar rhoL,scalar rhoG)
{
    if(!transmission()) return;
    if(mesh_.moving() || mesh_.topoChanging())
        FatalErrorInFunction << "plicTransmission requires a stationary mesh" << exit(FatalError);
    if (!(rhoL>0 && rhoG>0 && std::isfinite(rhoL) && std::isfinite(rhoG)))
        FatalErrorInFunction << "Positive finite phase densities required" << exit(FatalError);
    const label ni=mesh_.nInternalFaces();
    forAll(mesh_.cells(),i)
    {
        CellData& c=cells_.at(global_.toGlobal(i));
        c.pml=mag(sigma[i]);
        std::set<label> neighbours;
        for(const label f:mesh_.cells()[i])
        {
            if(f<ni) neighbours.insert(global_.toGlobal
                (mesh_.faceOwner()[f]==i ? mesh_.faceNeighbour()[f]:mesh_.faceOwner()[f]));
            else if(isA<processorPolyPatch>
                (mesh_.boundaryMesh()[mesh_.boundaryMesh().whichPatch(f)]))
                neighbours.insert(neighbourGlobal_[f-ni]);
        }
        c.adjacent=asList(neighbours);
    }
    std::set<label> wanted;
    std::map<label,std::vector<std::set<label>>> stencils;
    label badBoundary=0;
    for(const auto& entry:faces_)
    {
        const FaceData& f=entry.second;
        if(f.neighbourGlobal<0)
        {
            Pout<< "Interface meets physical patch at face " << f.face
                << " centre " << mesh_.faceCentres()[f.face] << nl;
            ++badBoundary; continue;
        }
        stencils[f.face].push_back({global_.toGlobal(f.owner),f.neighbourGlobal});
        wanted.insert(f.neighbourGlobal);
    }
    if(returnReduce(badBoundary,sumOp<label>()))
        FatalErrorInFunction << "plicTransmission does not support interface contact with "
            << "nonsymmetry physical boundaries" << exit(FatalError);
    exchangeCells(asList(wanted));
    for(label ring=1;ring<=maxRings_;++ring)
    {
        wanted.clear();
        for(auto& entry:stencils)
        {
            std::set<label> expanded=entry.second.back();
            for(const label id:entry.second.back())
            {
                const CellData& c=cells_.at(id);
                expanded.insert(c.adjacent.begin(),c.adjacent.end());
            }
            entry.second.push_back(expanded);
            wanted.insert(expanded.begin(),expanded.end());
        }
        exchangeCells(asList(wanted));
    }
    operators_.clear();
    label invalid=0, maxUsed=0;
    scalar worst=0;
    for(const auto& entry:faces_)
    {
        const FaceData& f=entry.second;
        const label og=global_.toGlobal(f.owner), ng=f.neighbourGlobal;
        if(f.neighbour<0 && og>ng) continue; // one builder per processor face
        if(cells_.at(og).pml>SMALL || cells_.at(ng).pml>SMALL)
        {
            Pout<< "Interface/PML overlap at face " << f.face << nl;
            ++invalid; continue;
        }
        bool success=false;
        for(label ring=1;ring<=maxRings_ && !success;++ring)
        {
            DynamicList<label> ids;
            ids.append(og); ids.append(ng);
            for(const label id:stencils.at(f.face)[ring])
                if(id!=og && id!=ng && cells_.at(id).pml<=SMALL) ids.append(id);
            pointField points(ids.size());
            forAll(points,i) points[i]=cells_.at(ids[i]).centre;
            scalarField weights(ids.size(),Zero);
            scalar condition=0;
            success=true;
            for(const Candidate& c:f.candidates)
            {
                FitResult fit=transmissionFit
                (
                    points,mesh_.faceCentres()[f.face],c.centre,c.normal,
                    c.area.liquid,c.area.gas,basis_,rhoL,rhoG,svdTol_,maxCondition_
                );
                if(!fit.valid) {success=false; break;}
                weights+=fit.weights/scalar(f.candidates.size());
                condition=max(condition,fit.condition);
            }
            if(success)
            {
                FaceOperator op;
                op.face=f.face; op.owner=f.owner; op.neighbour=f.neighbour;
                op.cells=labelList(ids); op.weights=weights;
                op.rings=ring; op.condition=condition;
                operators_.push_back(std::move(op));
                maxUsed=max(maxUsed,ring); worst=max(worst,condition);
            }
        }
        if(!success)
        {
            Pout<< "No valid transmission fit for face " << f.face << " centre "
                << mesh_.faceCentres()[f.face] << " after " << maxRings_ << " rings"
                << " ownerGlobal=" << og << " neighbourGlobal=" << ng
                << " stencil " << asList(stencils.at(f.face).back()) << nl;
            for (const Candidate& c:f.candidates)
                Pout<< "  plane centre " << c.centre << " normal " << c.normal
                    << " liquid area " << c.area.liquid << " gas area " << c.area.gas << nl;
            ++invalid;
        }
    }
    if(returnReduce(invalid,sumOp<label>()))
        FatalErrorInFunction << "Cannot construct plicTransmission operator; no legacy fallback."
            << exit(FatalError);
    shareOperators();
    reduce(maxUsed,maxOp<label>()); reduce(worst,maxOp<scalar>());
    Info<< "PLIC transmission: faces=" << returnReduce(label(operators_.size()),sumOp<label>())
        << " maxRingsUsed=" << maxUsed << " maxCondition=" << worst
        << " (processor faces counted on both ranks)" << nl;
}

void acousticInterface::shareOperators()
{
    if(!Pstream::parRun()) return;
    const auto& patches=mesh_.boundaryMesh();
    std::map<label,const FaceOperator*> byFace;
    for(const auto& op:operators_) byFace[op.face]=&op;
    PstreamBuffers buffers(Pstream::commsTypes::nonBlocking);
    forAll(patches,p) if(isA<processorPolyPatch>(patches[p]))
    {
        const auto& patch=refCast<const processorPolyPatch>(patches[p]);
        UOPstream os(patch.neighbProcNo(),buffers);
        DynamicList<label> faceIds;
        forAll(patch,i) if(byFace.count(patch.start()+i)) faceIds.append(i);
        os<<faceIds;
        for(const label i:faceIds)
        {
            const auto& op=*byFace.at(patch.start()+i);
            os<<op.cells<<op.weights<<op.rings<<op.condition;
        }
    }
    buffers.finishedSends();
    forAll(patches,p) if(isA<processorPolyPatch>(patches[p]))
    {
        const auto& patch=refCast<const processorPolyPatch>(patches[p]);
        UIPstream is(patch.neighbProcNo(),buffers);
        labelList faceIds; is>>faceIds;
        for(const label i:faceIds)
        {
            FaceOperator op; op.face=patch.start()+i;
            op.owner=mesh_.faceOwner()[op.face];
            is>>op.cells>>op.weights>>op.rings>>op.condition;
            op.weights*=-1;
            if(byFace.count(op.face))
                FatalErrorInFunction << "Duplicate processor transmission operator" << exit(FatalError);
            operators_.push_back(std::move(op));
        }
    }
    std::sort(operators_.begin(),operators_.end(),
        [](const FaceOperator& a,const FaceOperator& b){return a.face<b.face;});
}

void acousticInterface::mask(surfaceScalarField& coefficient) const
{
    for(const auto& op:operators_)
    {
        if(op.face<mesh_.nInternalFaces()) coefficient[op.face]=0;
        else
        {
            label p=mesh_.boundaryMesh().whichPatch(op.face);
            coefficient.boundaryFieldRef()[p][op.face-mesh_.boundaryMesh()[p].start()]=0;
        }
    }
}

void acousticInterface::write
(
    const surfaceScalarField& af,const volScalarField& pre,const volScalarField& pim
) const
{
    if(!diagnostics_) return;
    af.write();
    const fileName dir=mesh_.time().path()/"postProcessing"/"acousticInterface"/mesh_.time().timeName();
    mkDir(dir);
    OFstream geometry(dir/"geometry.tsv");
    geometry.precision(17);
    geometry<<"face\tcandidate\tliquidFraction\tliquidArea\tgasArea\tplaneCentre\tplaneNormal\n";
    for(const auto& entry:faces_)
    {
        label candidate=0;
        for(const auto& c:entry.second.candidates)
            geometry<<entry.first<<'\t'<<candidate++<<'\t'<<c.area.fraction<<'\t'
                <<c.area.liquid<<'\t'<<c.area.gas<<'\t'<<c.centre<<'\t'<<c.normal<<'\n';
    }
    OFstream os(dir/"operators.tsv");
    os.precision(17);
    os<<"face\townerGlobal\tneighbourGlobal\trings\tcondition\tcolumnGlobal\tweight\n";
    for(const auto& op:operators_) forAll(op.cells,k)
        os<<op.face<<'\t'<<global_.toGlobal(op.owner)<<'\t'
            <<(op.neighbour>=0 ? global_.toGlobal(op.neighbour):
                neighbourGlobal_[op.face-mesh_.nInternalFaces()])
            <<'\t'<<op.rings<<'\t'<<op.condition<<'\t'<<op.cells[k]<<'\t'<<op.weights[k]<<'\n';
    // Geometry is fixed. Only pressure values are exchanged at output time.
    std::set<label> remote;
    for(const auto& op:operators_) for(const label id:op.cells)
        if(!global_.isLocal(id)) remote.insert(id);
    List<DynamicList<label>> requests(Pstream::nProcs());
    for(label id:remote) requests[global_.whichProcID(id)].append(id);
    std::map<label,Vector2D<scalar>> pressure;
    forAll(pre,i) pressure[global_.toGlobal(i)]={pre[i],pim[i]};
    PstreamBuffers ask(Pstream::commsTypes::nonBlocking);
    forAll(requests,p) if(p!=Pstream::myProcNo() && !requests[p].empty())
        { UOPstream stream(p,ask); stream<<requests[p]; }
    ask.finishedSends();
    PstreamBuffers answer(Pstream::commsTypes::nonBlocking);
    forAll(requests,p) if(p!=Pstream::myProcNo() && ask.recvDataCount(p))
    {
        labelList ids; UIPstream request(p,ask); request>>ids;
        UOPstream stream(p,answer); stream<<ids.size();
        for(label id:ids) stream<<id<<pressure.at(id).x()<<pressure.at(id).y();
    }
    answer.finishedSends();
    forAll(requests,p) if(p!=Pstream::myProcNo() && answer.recvDataCount(p))
    {
        UIPstream stream(p,answer); label count; stream>>count;
        for(label i=0;i<count;++i)
        {
            label id; scalar a,b; stream>>id>>a>>b; pressure[id]={a,b};
        }
    }
    OFstream flux(dir/"flux.tsv"); flux.precision(17);
    flux<<"face\tPreFlux\tPimFlux\n";
    for(const auto& op:operators_)
    {
        scalar a=0,b=0;
        forAll(op.cells,k)
        {
            a+=op.weights[k]*pressure.at(op.cells[k]).x();
            b+=op.weights[k]*pressure.at(op.cells[k]).y();
        }
        flux<<op.face<<'\t'<<a<<'\t'<<b<<'\n';
    }
}
}}
