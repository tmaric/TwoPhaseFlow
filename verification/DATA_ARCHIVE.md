# Data archive: cyclic changes of frame in interface advection

This archive accompanies

> D. Bothe and T. Marić, *Cyclic changes of frame for designing non-retracing
> verification cases for numerical interface tracking*, arXiv preprint, 2026.

Archive DOI: <https://doi.org/10.6084/m9.figshare.35042791> (Figshare).
Code: TwoPhaseFlow, tag `arxiv-nonretracing-advection-v1`,
<https://github.com/tmaric/TwoPhaseFlow/tree/arxiv-nonretracing-advection-v1>.
`MANIFEST.txt` names the commit and lists a sha256 checksum of every file.

It holds the code, the workflows, the run data behind every table and figure,
and the scripts of the two illustrations that use no simulation. Figure and
table captions in the manuscript name the file below that each was produced from.

## Layout

| Path | Contents |
|---|---|
| `src/VoF/movingFrameFlow/`, `src/VoF/movingFrameFlow3D/` | laboratory-frame flux assembly (stream function, vector potentials) |
| `apps/benchmark/advectorVoF/` | the advection solver |
| `verification/reversed-vortex/` | Section 5: 2D reversed vortex |
| `verification/deformation-sphere3d/` | Section 6: 3D shear flow |
| `verification/leveque3d/` | Section 7: 3D LeVeque deformation |
| `verification/<study>/results/` | every figure, table and check the manuscript quotes |
| `verification/<study>/runs/<case>/` | per run: the error time series `volumeFractionError.dat` and the case's `system/` |
| `verification/manuscript-illustrations/` | Figure 1 and the Section 2 checks (no simulation) |
| `verification/postprocess.py` | regenerates tables and convergence figures from `runs/` |
| `verification/requirements.txt` | Python packages, with the versions used |
| `MANIFEST.txt` | commit, tag, build date and a checksum of every file |

Case names are `<scheme>_<frame>_N<resolution>`: scheme `plicRDF` or `isoAlpha`;
frame `none` (the original, retracing test), `frameA` (one revolution
superposed), `frameId` (the identity frame) or `frameOnly` (the frame without
the base field). `viz_*` cases write more often and supply the interface figures
only. Solver logs and volume-fraction fields are not included (up to
$10^2$ MB per run); rerunning a case regenerates them.

## Results per study

| File | Manuscript |
|---|---|
| `results/convergence.csv`, `convergence_table.tex` | convergence tables |
| `results/convergence_shape_error.{pdf,png}` | convergence figures |
| `results/interface*_moving_frame.png`, `*_non_moving_frame.png` | interface snapshots |
| `results/timestep_bound.json` | the bound max‖w‖₁ that fixes the time step |
| `results/courant_check.txt` | per run: Courant number, discrete divergence, volume and boundedness errors reached |
| `reversed-vortex/results/frame_identity_check.txt` | identity frame against the original code path |
| `reversed-vortex/results/frame_only_check.txt` | the frame alone (α ≡ 0) |
| `reversed-vortex/results/interface_evolution.{pdf,png,gif}` | mirror-pair snapshots against the exact interface; animation |
| `<3D study>/results/potential_check.txt` | curl of the vector potential against the base field |
| `<3D study>/results/support_extent.txt` | how far the material reaches: the domain bound |

## Reproducing

Three levels, from cheapest to complete.

1. **Check a number.** Every value quoted in the manuscript is in
   `verification/<study>/results/` (CSV, `.txt`, `.json`).
2. **Regenerate tables and convergence figures from the run data** (minutes,
   no OpenFOAM):

   ```bash
   python3 -m venv .venv && . .venv/bin/activate
   pip install -r verification/requirements.txt
   python3 verification/postprocess.py            # writes results-regenerated/
   diff verification/reversed-vortex/results/convergence.csv \
        results-regenerated/reversed-vortex/convergence.csv
   ```

   The two illustrations: `cd verification/manuscript-illustrations` and
   `python3 cancellation_1d.py`, `python3 two_motions_storyboard.py OUT`.
3. **Rerun the studies.** Build TwoPhaseFlow against OpenFOAM v2512, then, in
   each study directory,

   ```bash
   ./run.sh            # SLURM if sbatch is available, otherwise local
   ```

   See `verification/<study>/README.md` for the protocol and requirements
   (OpenFOAM v2512, Snakemake 9). The published runs used the Lichtenberg
   cluster at TU Darmstadt; the finest 2D and 3D levels need an HPC system.

## Licence

The code is distributed under the GNU General Public License v3 (`LICENSE`),
as TwoPhaseFlow and OpenFOAM are. The results and run data are released under
the licence stated on the archive's Figshare page.
