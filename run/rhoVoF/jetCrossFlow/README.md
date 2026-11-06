# Liquid jet in high-speed gaseous cross flow

The case of section 4.6 of

> J. Liu, T. Tolle, D. Zuzio, J.-L. Estivalèzes, S. Marquez Damian, T. Marić,
> Inconsistencies in unstructured geometric volume-of-fluid methods for two-phase
> flows with high density ratios, Computers & Fluids 281 (2024) 106375,
> doi:10.1016/j.compfluid.2024.106375 (preprint arXiv:2312.10434),

as published on the branch `feature/density-ratio` of this repository (commit
`d6abfb1`) and in the data archive T. Marić, J. Liu, *Inconsistencies in
unstructured volume-of-fluid methods for two-phase flows with high density ratios -
data*, TUdatalib (2023), doi:10.48328/tudatalib-1277 (the code archive:
doi:10.48328/tudatalib-1276). The setup follows

> D. Zuzio, A. Orazzo, J.-L. Estivalèzes, I. Lagrange, A new efficient momentum
> preserving Level-Set/VOF method for high density and momentum ratio
> incompressible two-phase flows, J. Comput. Phys. 410 (2020) 109342,
> doi:10.1016/j.jcp.2020.109342

(the paper's reference [25]); the experiments are ONERA's,

> A. Desclaux, S. Thuillet, D. Zuzio, J.-M. Senoner, D. Sebbane, V. Bodoc,
> P. Gajan, Experimental and numerical characterization of a liquid jet injected
> into air crossflow with acoustic forcing, Flow Turbul. Combust. 105 (2020)
> 1087-1117, doi:10.1007/s10494-020-00126-0 (reference [58]);
> V. Bodoc, A. Desclaux, P. Gajan, F. Simon, G. Illac, Characterization of
> confined liquid jet injected into oscillating air crossflow, Flow Turbul.
> Combust. 104 (2020) 1-18, doi:10.1007/s10494-019-00037-9 (reference [59]);

and the liquid inlet velocity profile follows

> S. Thuillet, Simulation multi-échelle de l'atomisation d'un jet liquide sous
> l'effet d'un écoulement gazeux transverse en présence d'une perturbation
> acoustique, PhD thesis, Université de Toulouse, 2018,
> https://hal.archives-ouvertes.fr/tel-02023561 (reference [57]).

(The citations were taken from the paper's bibliography, arXiv:2312.10434v2, on
2026-11-06.)

## The case

- Domain [-0.01, 0, -0.01] x [0.03, 0.02, 0.01] m (`system/blockMeshDict` gives
  the box; the mesh is cfMesh's `cartesianMesh` of `flowBoxGeo/jetCrossFlowBox.fms`,
  `system/meshDict` = `cartesianMeshDict` with `maxCellSize 0.02/128` for the
  256 x 128 x 128 mesh of the paper's figures 16/17, `cartesianMeshDict.coarse` with
  `0.02/64` for 128 x 64 x 64, the paper's "moderate resolution").
- The nozzle (d = 2 mm) at the origin of the bottom patch: `topoSet` + `createPatch`
  cut the patch `jetInlet` out of `airBed` (`initialCase.sh`: blockMesh,
  cartesianMesh, topoSet, createPatch -overwrite, decomposePar).
- Gas inflow 65 m/s from `inlet` (x = -0.01); the liquid inlet velocity profile
  v_y(r) = -21.434 (r/d)^3 + 15.512 (r/d)^2 - 2.5722 (r/d) + 8.6504 m/s (`0/U`,
  `codedFixedValue`, compiled at run time); water/air
  (`constant/transportProperties`), surface tension 0.072 N/m, gravity -9.81 m/s^2
  in y.
- The published solver setting (`system/fvSolution`, `system/controlDict`):
  `interFlow` (the p-U coupling of this repository; `feature/density-ratio` adds
  only the opt-in `-resNon` option), isoAdvection + plicRDF with `isoFaceTol
  1e-8`, `nAlphaBounds 3`, `snapTol 0`, `clip false`, PIMPLE 2 outer / 3 pressure
  correctors, fixed dt 5e-7 s (CFL about 0.2 at the gas inlet), end time 8 ms,
  writes every 0.1 ms; `decomposeParDict` simple (8 2 2) on 32 ranks
  (`isoAdv.sbatch`). The function objects `FO_plicInterface` (the PLIC interface
  as a vtp surface at the writes) and `FO_visRDF` (the RDF field of plicRDF; a
  coded function object whose include path is the repository's
  `src/VoF/lnInclude`, as `$(FOAM_USER_LIBBIN)/../../../src/VoF/lnInclude` under
  `etc/davof-env.sh`).
- `CXY_q7.51_calc.csv`: the experimental windward shape of the jet (ONERA, the
  references above; the momentum-flux ratio q = 7.51 of the experiment against the
  case's 7.8), columns X, Y, Z in metres, the nozzle centre at the origin.
- `results_figures/`: the published figures and ParaView states (LJCF_*; the
  comparison figures overlay the experimental points on the alpha = 0.5
  iso-surface at t = 7.1 ms).

## The comparison with the experiment (`comparison/`)

`comparison/windward_extract.py` (ParaView's `pvbatch`) reads the case (decomposed
or reconstructed) at a time, builds the alpha.water = 0.5 iso-surface, cuts it with
the symmetry plane z = 0 and writes the windward edge (the smallest x of the cut at
every height y) to `comparison/windward_<time>.csv`; `comparison/windward_compare.py`
(matplotlib) plots the windward edges of the runs against `CXY_q7.51_calc.csv` and
prints, per run, the largest horizontal deviation from the experiment over the
experiment's heights (the paper: "around 1 mm, which is 5 % of the jet height" at
t = 7.1 ms on the fine mesh). The paper's figures were made in ParaView
(`results_figures/LJCF_state*.pvsm`); the scripts here reproduce the comparison
without a ParaView session.

```
module load paraview/5.13.3-rzfn          # Lichtenberg
pvbatch comparison/windward_extract.py --case <case> --time 0.0071 --out comparison/windward_RDF.csv
python3 comparison/windward_compare.py --exp CXY_q7.51_calc.csv RDF=comparison/windward_RDF.csv RING=... --out comparison/windward.png
```

## The verification of the isoAdvector package (branch `davof-isoAdvector`, 2026-11-06)

The rows plicRDF with `positioning vofCutCell` (the published setting), plicRDF
with `positioning ccs`, and the `davof` reconstruction with the ring-area normal
and `positioning ccs` on the same mesh, dt, decomposition and ranks; the mesh
check and the static normal test first (leia `docs/davof/REPORT-swept.md`). The
`davof` row needs the face fractions `alphaf.water` at t = 0 (uniform 0: the
domain holds gas; the nozzle patch's values follow the transported state, the
face update takes a boundary face's fraction from the owner's plane, not from the
boundary condition).
