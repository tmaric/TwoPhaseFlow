#!/usr/bin/env pvbatch
"""The windward edge of the liquid jet in cross flow (run/rhoVoF/jetCrossFlow): the
alpha.water = 0.5 iso-surface of the case at one time, cut by the symmetry plane
z = 0; the smallest x of the cut at every height y (bins of dy) is the windward
shape, written as a CSV (x, y in metres) for windward_compare.py against the
experiment (CXY_q7.51_calc.csv). ParaView's pvbatch (OpenFOAMReader):

    pvbatch windward_extract.py --case <case dir> --time 0.0071 --out windward.csv [--dy 1e-4] [--decomposed]
"""
import argparse, csv, os, sys
from paraview.simple import *
from paraview import servermanager
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy

ap = argparse.ArgumentParser()
ap.add_argument("--case", required=True)
ap.add_argument("--time", type=float, required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--dy", type=float, default=1e-4)
ap.add_argument("--decomposed", action="store_true", help="read processor*/ (default: reconstructed or auto)")
ap.add_argument("--zmax", type=float, default=0.0, help="half width of the band |z| <= zmax of the cut (0: the plane slice)")
a = ap.parse_args()

foam = os.path.join(a.case, "case.foam")
open(foam, "a").close()
r = OpenFOAMReader(FileName=foam)
r.CaseType = "Decomposed Case" if a.decomposed else "Reconstructed Case"
r.MeshRegions = ["internalMesh"]
r.CellArrays = ["alpha.water"]
r.UpdatePipelineInformation()
times = list(r.TimestepValues) if r.TimestepValues else []
t = min(times, key=lambda x: abs(x - a.time)) if times else a.time
print("times", len(times), "using", t)
r.UpdatePipeline(t)
c2p = CellDatatoPointData(Input=r)
c2p.CellDataArraytoprocess = ["alpha.water"]
iso = Contour(Input=c2p)
iso.ContourBy = ["POINTS", "alpha.water"]
iso.Isosurfaces = [0.5]
iso.UpdatePipeline(t)
if a.zmax > 0:
    clip = Clip(Input=iso); clip.ClipType = "Box"
    cut = clip
else:
    sl = Slice(Input=iso); sl.SliceType = "Plane"; sl.SliceType.Origin = [0, 0, 0]; sl.SliceType.Normal = [0, 0, 1]
    cut = sl
cut.UpdatePipeline(t)
d = servermanager.Fetch(cut)
pts = vtk_to_numpy(d.GetPoints().GetData()) if d.GetNumberOfPoints() else np.zeros((0, 3))
if hasattr(d, "GetNumberOfBlocks"):
    # multiblock: gather
    from paraview.vtk import vtkCompositeDataIterator
    allp = []
    it = d.NewIterator(); it.InitTraversal()
    while not it.IsDoneWithTraversal():
        b = it.GetCurrentDataObject()
        if b is not None and b.GetNumberOfPoints():
            allp.append(vtk_to_numpy(b.GetPoints().GetData()))
        it.GoToNextItem()
    pts = np.vstack(allp) if allp else np.zeros((0, 3))
print("cut points", len(pts))
rows = []
if len(pts):
    y = pts[:, 1]; x = pts[:, 0]
    bins = np.floor(y/a.dy).astype(int)
    for b in np.unique(bins):
        m = bins == b
        rows.append((x[m].min(), (b + 0.5)*a.dy, x[m].max(), int(m.sum())))
with open(a.out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["x_windward", "y", "x_leeward", "n_points"]); w.writerows(rows)
print("wrote", a.out, len(rows), "heights; jet height (max y of the cut) %.4f m" % (pts[:, 1].max() if len(pts) else float("nan")))
