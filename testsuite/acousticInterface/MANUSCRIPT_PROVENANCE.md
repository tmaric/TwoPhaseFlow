# Manuscript data provenance and geometric-averaging evidence

The revised methods paper selects `areaFraction plicAverage; flux legacy;`. All heterogeneous results now included in the paper have explicit geometric-method evidence. Historical single-fluid wave, piston and rigid-sphere results are retained and identified as such: both area models give zero liquid face fraction in these domains. They are not presented as newly rerun geometric cases.

This audit covers **19 figures and 7 tables** in `helmHoltzPML-revised.tex`. The original manuscript and all 23 original external assets are unchanged. Updated assets are confined to `revised-assets/` in the article workspace.

## Changes made from the audit

- Added Section 4.4 and Table 4 comparing all four Cartesian sphere resolutions, pressure errors and serial/eight-rank differences. On 48³ cells, pressure error is 0.02974595 versus 0.02975350, while the decomposition difference is 2.6395e-6 versus 4.9163e-15. The supported benefit is reproducibility with retained accuracy.
- Replaced the stale homogeneous table and plot with the archived results for the mesh that retains orthogonal boundary layers. Used the explicitly volume-weighted study; finest velocity orders are 1.67 and 1.68.
- Replaced both piston tables and their plots with the archived fixed-piston-edge refinements. The far-field table now reports pressure-magnitude error, as defined in the text.
- Rebuilt all layered plots and both sensitivity tables from selected-method fields. Added the three intermediate damping comparisons: 12 runs on 1/4 ranks, all passed. There are now 44 local plus 32 cluster layered runs. Maximum relative pressure difference from the frozen original serial executable remains 1.95e-12.
- Corrected the finest PML real-pressure component error from the stale displayed 3.071e-5 to 3.074e-5; the whole complex-pressure error remains 2.657e-4.
- Traced the rigid-sphere table and all three force plots to archived force records (3 mesh levels, 7 radii, 11 positions). Checked the quoted 1.485% radius maximum and 1.288% position RMS against those records.
- Removed the untraceable 0.0602 fN finer-mesh zero-crossing result and the causal interpretation based on it. The recorded 0.0633 fN residual remains, with a qualified interpretation.

## Inventory

| Item | LaTeX label | Provenance / status |
|---|---|---|
| Figure 1 | `fig:standardFVMControlVolume` | Original illustration, hash preserved; visually checked against the described geometry. Not solver-run evidence. |
| Figure 2 | `fig:unstructuredFVMSchematic` | Native LaTeX/TikZ schematic; checked against equations, geometry and boundary conditions. No simulation data. |
| Figure 3 | `fig:homogeneousBaselineSchematic` | Native LaTeX/TikZ schematic; checked against equations, geometry and boundary conditions. No simulation data. |
| Figure 4 | `fig:homogeneousMeshes` | Original illustration, hash preserved; visually checked against the described geometry. Not solver-run evidence. |
| Figure 5 | `fig:homogeneousConvergence` | Archived Lichtenberg job 54501269: 24 cases; cell-volume-weighted metrics. Updated table and plot; retained single-phase results. |
| Figure 6 | `fig:combinedLayeredSchematic` | Native LaTeX/TikZ schematic; checked against equations, geometry and boundary conditions. No simulation data. |
| Figure 7 | `fig:combinedLayeredPressure` | Selected plicAverage/legacy fields: all 11 distinct physical configurations checked locally; eight also checked in job 54506190. 1200 valid samples per case. |
| Figure 8 | `fig:combinedLayeredAmplitude` | Selected plicAverage/legacy fields: all 11 distinct physical configurations checked locally; eight also checked in job 54506190. 1200 valid samples per case. |
| Figure 9 | `fig:combinedLayeredWaterAir` | Selected plicAverage/legacy fields: all 11 distinct physical configurations checked locally; eight also checked in job 54506190. 1200 valid samples per case. |
| Figure 10 | `fig:pmlSigmaSensitivity` | Selected plicAverage/legacy fields: all 11 distinct physical configurations checked locally; eight also checked in job 54506190. 1200 valid samples per case. |
| Figure 11 | `fig:pmlMeshSensitivity` | Selected plicAverage/legacy fields: all 11 distinct physical configurations checked locally; eight also checked in job 54506190. 1200 valid samples per case. |
| Figure 12 | `fig:pistonConfiguration` | Native LaTeX/TikZ schematic; checked against equations, geometry and boundary conditions. No simulation data. |
| Figure 13 | `fig:pistonNearAmpMethod` | Archived Lichtenberg job 54501043: five meshes with fixed piston edge; source comparison CSVs and metrics. Updated table/plot. |
| Figure 14 | `fig:pistonFarAmpMethod` | Archived Lichtenberg job 54501043: five meshes with fixed piston edge; source comparison CSVs and metrics. Updated table/plot. |
| Figure 15 | `fig:pistonFarFieldPattern` | Archived Lichtenberg job 54501043: five meshes with fixed piston edge; source comparison CSVs and metrics. Updated table/plot. |
| Figure 16 | `fig:gorkovDomainSchematic` | Original illustration, hash preserved; visually checked against the described geometry. Not solver-run evidence. |
| Figure 17 | `fig:gorkovMeshConvergence` | Archived force TSVs in the repository; reprocessed in job 54501044. No new force simulation claimed. |
| Figure 18 | `fig:gorkovRadiusSweep` | Archived force TSVs in the repository; reprocessed in job 54501044. No new force simulation claimed. |
| Figure 19 | `fig:gorkovPositionSweep` | Archived force TSVs in the repository; reprocessed in job 54501044. No new force simulation claimed. |
| Table 1 | `tab:homogeneousConvergence` | Archived Lichtenberg job 54501269: 24 cases; cell-volume-weighted metrics. Updated table and plot; retained single-phase results. |
| Table 2 | `tab:pmlSigmaSensitivity` | Selected plicAverage/legacy fields: all 11 distinct physical configurations checked locally; eight also checked in job 54506190. 1200 valid samples per case. |
| Table 3 | `tab:pmlMeshSensitivity` | Selected plicAverage/legacy fields: all 11 distinct physical configurations checked locally; eight also checked in job 54506190. 1200 valid samples per case. |
| Table 4 | `tab:geometricComparison` | Lichtenberg job 54503635_2: N=16,24,32,48; both area methods on 1/8 ranks. Original recorded study rows. |
| Table 5 | `tab:pistonConvergence` | Archived Lichtenberg job 54501043: five meshes with fixed piston edge; source comparison CSVs and metrics. Updated table/plot. |
| Table 6 | `tab:pistonFarFieldConvergence` | Archived Lichtenberg job 54501043: five meshes with fixed piston edge; source comparison CSVs and metrics. Updated table/plot. |
| Table 7 | `tab:gorkovMeshConvergence` | Archived force TSVs in the repository; reprocessed in job 54501044. No new force simulation claimed. |

## Verification and limits

- 16 homogeneous table rows match volume-weighted source; both pressure and velocity orders checked
- 10 piston table rows match the geometry-aligned near-field and pressure-magnitude far-field source
- All 10 PML sensitivity rows match geometric-method fields sampled at 1200 valid points
- 3 rigid-sphere mesh table rows match archived integrated forces, face counts and relative errors
- All 16 entries in the new sphere comparison table match the original completed study rows

The full source-to-asset mapping, SHA-256 hashes and numeric checks are in [results/manuscript-provenance.json](results/manuscript-provenance.json). Compact numerical records are in [results/manuscript-evidence.json](results/manuscript-evidence.json). Raw historical evidence was retrieved from `/work/scratch/tm83tomy/TwoPhaseFlowAcoustics/.review-testing/method-review-corrections/` and is retained locally under `.review-testing/submission-evidence-audit/historical/`. New damping fields and logs are under `.review-testing/submission-evidence-audit/remaining-damping/`.

The archive does not contain a complete original field/binary record for the rigid-sphere force sweeps. Their provenance reaches the committed integrated-force TSVs and reproducible postprocessing, not a new full solver rerun. Original mesh illustrations are retained as illustrations; image hashes alone are not mesh/run provenance. No stronger claim is made in the paper.

The archived homogeneous postprocessor can fall back to uniform weights if cell volumes are absent. The revised manuscript specifically uses the `volume-54501269` study, which wrote cell volumes before computing errors. Future reproduction must retain that step.

Run the remaining damping check in a fresh directory after loading the documented isolated environment:

```sh
python3 testsuite/acousticInterface/submission_check.py \
  --output /absolute/new/remaining-damping \
  --baseline /absolute/path/acousticHelmholtzFoam-baseline \
  --remaining-damping --ranks 1,4
```

`ROADMAP.md` remains locally ignored and untracked. The experimental transmission flux is not used by the manuscript comparisons.

## Stable publication data location

The audited secondary data have now been gathered in [data/publications/helmholtz-pml](../../data/publications/helmholtz-pml/README.md). [The root provenance file](../../PROVENANCE.md) provides a repository-relative source link for every manuscript figure and table. The paths above remain historical audit origins; they are no longer required to use the publication dataset.
