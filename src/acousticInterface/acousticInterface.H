#ifndef acousticInterface_H
#define acousticInterface_H
#include "fvCFD.H"
#include "globalIndex.H"
#include "transmissionFit.H"
#include <map>
#include <vector>

namespace Foam { namespace acoustic
{
struct FaceOperator
{
    label face=-1, owner=-1, neighbour=-1;
    labelList cells;
    scalarField weights;
    label rings=0;
    scalar condition=0;
};
class acousticInterface
{
    struct CellData
    {
        point centre=point::zero, planeCentre=point::zero;
        vector normal=vector::zero;
        scalar alpha=0, pml=0;
        labelList adjacent;
    };
    struct Candidate
    {
        point centre;
        vector normal;
        CutArea area;
    };
    struct FaceData
    {
        label face=-1, owner=-1, neighbour=-1, neighbourGlobal=-1;
        std::vector<Candidate> candidates;
    };
    const fvMesh& mesh_;
    const globalIndex& global_;
    word areaName_, fluxName_;
    bool diagnostics_;
    label maxRings_;
    scalar svdTol_, maxCondition_, pureTol_=1e-8;
    std::map<label,CellData> cells_;
    std::map<label,FaceData> faces_;
    std::vector<FaceOperator> operators_;
    vectorField basis_;
    labelList neighbourGlobal_;
    void exchangeCells(const labelList&);
    void computeGeometricAreas(volScalarField&,surfaceScalarField&,cutFaceAdvect&);
    void shareOperators();
public:
    acousticInterface(const fvMesh&,const globalIndex&);
    bool transmission() const { return fluxName_=="plicTransmission"; }
    bool writeDiagnostics() const { return diagnostics_; }
    const std::vector<FaceOperator>& operators() const { return operators_; }
    void computeAreas(volScalarField&,surfaceScalarField&,const surfaceScalarField&,cutFaceAdvect&);
    void build(const volTensorField& sigma,scalar rhoL,scalar rhoG);
    void mask(surfaceScalarField&) const;
    void write(const surfaceScalarField&,const volScalarField&,const volScalarField&) const;
};
}}
#endif
