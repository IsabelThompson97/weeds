# WEeDS — Weighted Ensemble equilibrium-mode Dispersion Seeding

A pipeline for estimating slow modes of RNA conformational dynamics by combining bidirectionally-seeded weighted ensemble (WE) simulations with a perturbation-mode dispersion strategy. Trajectory segments from many short WE runs are pooled, featurized at the walker start/end pair level, projected with tICA, and reweighted with RiteWeight to recover equilibrium spectral estimates.

This repository holds everything from basis-state preparation through the WE simulation drivers that produce the segment corpus consumed by the downstream tICA + RiteWeight analysis. The full theoretical and methodological background is in `WEeDS.md` (and the rendered PDF alongside it).

<img width="1040" height="780" alt="pdist_evolution" src="https://github.com/user-attachments/assets/ef789e61-17c0-47ea-8e2d-d89be7971708" />


---

## Scientific objective

The object of interest is not any individual trajectory but the time-evolving probability distribution of the weighted ensemble itself. Concretely, we want to estimate the leading non-trivial eigenpairs of the transfer (Koopman) operator on a featurized state space — the slow modes that govern folding/unfolding and intermediate metastability. Two facts about WE shape the pipeline:

1. A single WE run is non-equilibrium: walkers carry a history of resampling events whose marginal distribution depends on the seeding choice and the progress coordinate. The corpus is *importance-weighted*, not Boltzmann-distributed.
2. Slow modes are only identifiable in the window of observability set by the lag time, the segment length, and the coverage of `supp(π)`. A single direction of seeding tends to over-resolve modes aligned with the progress coordinate and attenuate orthogonal ones.

The strategy in this repo addresses both: bidirectional seeding (folded → unfolded *and* unfolded → folded) for coverage, and short perturbation runs initialized from a dispersed set of intermediates to excite modes orthogonal to the progress coordinate. The resulting multi-corpus is then pooled and reweighted.

---

## Pipeline structure

```
Source trajectories (control MD or prior WE)
        │
        ▼
[bstates/]            Prepare basis states
                      → solvated, HRM-corrected, equilibrated struct_X/ directories
                      → bstates.txt indexed by progress coordinate
        │
        ▼
[we_bidirectional/]   WE setup + run scripts for bidirectional seeding   (planned)
                      → folded→unfolded and unfolded→folded WE runs
                      → west.cfg, system.py, runseg.sh, env.sh per direction
        │
        ▼
[we_perturbation/]    Perturbation-mode short WE / dispersion runs        (planned)
                      → many short runs seeded from dispersed intermediates
                      → designed to excite modes orthogonal to the PC
        │
        ▼
[features/]           Per-segment featurization                           (planned)
                      → one feature vector per walker (start, end) pair
                      → pooled across both seeding directions and all perturbation runs
                      → walker weights tracked for downstream reweighting
        │
        ▼
[tica/]               tICA dimensionality reduction                       (planned)
                      → fit on pooled, weighted segment features
                      → produces slow collective variables
        │
        ▼
[riteweight/]         RiteWeight reweighting                              (planned)
                      → recovers equilibrium spectral estimates from the
                        non-equilibrium WE corpus
                      → two configurations (see Background §3.3)
```

Each stage writes its outputs into a versioned subdirectory and a corresponding `PROJECT_LOG.md` that records the iteration count, walker counts, and any parameter changes made between runs.

---

## Current contents

| Path | Status | Description |
|------|--------|-------------|
| `bstates/` | **Implemented** | Basis-state preparation workflow. Two entry paths (control-MD frame extraction or mining prior `west.h5` for intermediates), shared solvation + HRM correction, staged AMBER minimization/equilibration via SGE array jobs, and final assembly into a WESTPA-ready `bstates/` folder. See `bstates/README.md`. |
| `WEeDS_Background_reorganized.md` / `.pdf` | **Reference** | Theoretical foundations, pipeline architecture, mode-excitation analysis, identifiability conditions, and the rationale for bidirectional + perturbation seeding and pooled-corpus reweighting. |
| `WEeDS_Background_audit_report.md` | **Reference** | Adversarial audit notes on the background document. |
| `052326_WE_BidirectionalSeedingBackground.md` | **Reference** | Earlier-dated background notes on the bidirectional seeding rationale. |

## Planned contents

| Path | Purpose |
|------|---------|
| `we_bidirectional/` | WESTPA project directories for the folded→unfolded and unfolded→folded simulations. Each will hold `west.cfg`, `system.py` (progress-coordinate definition, bin mapper, target state), `env.sh`, `run.sh`, `runseg.sh`, and the `bstates.txt` produced by the `bstates/` workflow. Drivers will support the WEED (equilibrium-mode) recycling configuration described in Background §4.5. |
| `we_perturbation/` | Setup and SGE submission scripts for the perturbation-mode dispersion runs: many short WE simulations initialized from a dispersed set of intermediate structures, designed to excite modes orthogonal to the primary progress coordinate. |
| `features/` | Per-segment featurization scripts. For each walker, the start-frame and end-frame structures are featurized (RNA-specific internal coordinates — distances, torsions, base-pairing descriptors — to be finalized), producing one (start, end) feature pair per segment, tagged with the segment's WE weight and provenance (run ID, iteration, walker ID). |
| `tica/` | tICA fitting on the pooled, weighted segment corpus, with lag time(s) and feature set chosen per the identifiability discussion in Background §4.12. |
| `riteweight/` | RiteWeight reweighting of the pooled corpus to recover equilibrium spectral estimates. Both configurations (Background §3.3) will be supported. |
| `analysis/` | Downstream notebooks/figures: implied timescales, mode visualization, coverage diagnostics, destructive-interference checks (Background §4.14). |

---

## Conventions

- **HPC**: Notre Dame CRC, SGE scheduler. GPU queues for `pmemd.cuda`; CPU queues for setup, equilibration assembly, and analysis.
- **Filesystem**: large run outputs (solvated topologies, `traj_segs/`, `seg_logs/`, equilibration outputs) live under `bstates/0[1-3]_*/` and are git-ignored. Only scripts, configs, and small text outputs are tracked.
- **Software stack**: AMBER 24.0 (`pmemd.cuda`, `cpptraj`, `parmed`, `tleap`), WESTPA 2022.10, Python 3 with `h5py`/`numpy`/`mdtraj` for analysis. Load with `module load amber/24.0` (and the WESTPA module before any WE-stage script).
- **Logging**: every stage maintains a `PROJECT_LOG.md` that records iteration progress, parameter changes, and any non-default decisions (bin remapping, recycling toggles, lag-time choices).

---

## Where to start

- To prepare basis states for a new system: `bstates/README.md`.
- To understand *why* the pipeline is structured around bidirectional + perturbation seeding and pooled reweighting: `WEeDS_Background_reorganized.pdf`, §3 (Pipeline Architecture) and §4.6, §4.11, §4.13.
- To see what changes between iterations of a given stage: the `PROJECT_LOG.md` inside that stage's directory.
