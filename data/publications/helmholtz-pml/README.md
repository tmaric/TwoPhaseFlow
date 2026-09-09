# Secondary data for the Helmholtz/PML methods paper

This is the portable secondary-data deposit for **An unstructured finite-volume Helmholtz method with perfectly matched layers for heterogeneous two-phase acoustics**. It covers all **19 figures and 7 tables** in the included `helmHoltzPML-revised.tex` snapshot. No access to the original WSL or Lichtenberg working directories is needed to read the data, check the tables, or reproduce the numerical plots.

The data support choosing geometric PLIC face averaging with the previous pressure flux. On the 48-cubed penetrable-sphere mesh, the pressure error changes from 0.02974595 to 0.02975350 while the serial/eight-rank difference falls from 2.6395e-6 to 4.9163e-15. This supports reproducibility with retained pressure accuracy, not improved force accuracy.

## Contents

- `secondary/`: every plotted sampled pressure profile, the convergence and error tables, numerical and analytical radiation-force data, and supporting comparison statistics.
- `methods.json`: physical parameters, numerical selections, normalization, sampling, simulation status and recorded cluster job IDs.
- `DATA_DICTIONARY.md` and `data-dictionary.json`: columns, units, missing-value semantics, row counts and types.
- `provenance.json`: one entry per figure/table, with paths relative to this deposit. `PROVENANCE.md` provides a readable mapping.
- `software-provenance.json`: source-file and binary hashes for each measured two-phase study and the final implementation.
- `origins.json`: collection lineage and original source hashes. Historical `repository:` locations describe origins only; they are not dependencies of this deposit.
- `manuscript/`: frozen LaTeX identification snapshot, exactly the assets used by it, seven original table files, bibliography and native schematic excerpts.
- `scripts/reproduce.py`: regenerates all 12 numerical figures and 7 tables using only the included CSV files. Static schematic/mesh illustrations are retained in their original form.
- `scripts/validate.py`: checks hashes, paths, schemas, counts, finiteness, numerical error measures and every displayed table row.
- `scripts/original/`: the original postprocessing sources, including the archived homogeneous/piston/Gorkov versions; field-to-sample stages in those historical scripts require the original solver environment.
- `metadata-draft.json`: TUdatalib metadata prepared for author review; no DOI, license, collection or dataset-creator assignment has been invented.
- `RIGHTS.md`: distinguishes data and software licenses, including the CC-BY-4.0 license of the original TikZ Figure 1.

## Use without OpenFOAM

Python 3.10 or newer with NumPy and Matplotlib is sufficient for the portable tools:

```sh
python3 scripts/validate.py
python3 scripts/reproduce.py --output ../helmholtz-pml-reproduced
```

The output directory must be new and outside the deposit, so the original payload and checksums remain unchanged. The reproduced plots have the same numerical content; pixel-level identity with historical Matplotlib renderings is not promised. The manuscript assets preserve the exact original appearance.

## What is included and excluded

The 33 CSV files contain all pressure samples and aggregate results used in the manuscript. The homogeneous and piston benchmarks use archived single-phase results, for which both area methods give the same coefficients. The rigid-sphere force records are archived integrated forces; a complete original raw-field/binary archive has not been recovered. The layered and penetrable-sphere comparisons explicitly exercise geometric averaging. These distinctions are stated in `methods.json` and the manuscript.

Raw volume meshes, OpenFOAM cell fields, processor partitions and solver binaries are excluded. The unselected experimental transmission-flux results are outside this manuscript deposit. The included manuscript is a frozen identification snapshot, not a claim that every historical study was rerun. Supporting text-only data include the boundary-skew diagnostic, 76 layered method-comparison runs, sphere matrix entry counts, and the archived force-amplitude scaling check.

## TUdatalib preparation

The title, description, data year, methods and file documentation are prepared. Final dataset creators, collection, DFG subject selection, license and any related article DOI remain deposit metadata choices. The current data tree has more than 20 files and depends on its directory structure, so a ZIP export is provided alongside the uncompressed tree. TUdatalib's guidance allows archives for preserving folder structure or larger file collections. [TUdatalib FAQ](https://tudatalib.ulb.tu-darmstadt.de/docs/en/faq/)

The metadata draft follows the fields described in the [TUdatalib user guide](https://tudatalib.ulb.tu-darmstadt.de/docs/en/nutzer_leitfaden/). No upload, DOI registration or publication is performed by these tools.

## Version identity

`checksums.json` fixes the exact payload. `provenance.json` identifies the manuscript source/PDF hashes, the repository branch and the baseline revision. The implementation commit is recorded separately from its baseline; recorded source and binary hashes distinguish the measured historical versions. This deposit is version `1.0-draft`, assembled on 2026-09-09.

To create another checked ZIP in a new location, run `python3 scripts/make_release.py --output ../helmholtz-pml-secondary-data.zip`. The archive has a fixed file order and timestamps for reproducible transfer.
