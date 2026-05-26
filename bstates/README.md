# prep_bstates

Workflow for preparing WESTPA basis states (bstates) from existing MD trajectories. Structures are extracted from production runs or a prior weighted ensemble simulation, solvated to a uniform water count, force-field-corrected with the Hydrogen Reparametrization Method (HRM), and equilibrated through a staged AMBER minimization and equilibration protocol using SGE array jobs on the CRC cluster.

---

## Overview

A basis state for a WESTPA simulation must be a fully solvated, equilibrated structure with a known progress coordinate value. This workflow takes raw structures from one of two sources and produces a set of `struct_X/` directories that drop directly into a WESTPA `bstates/` folder.

**Two entry paths:**

- **Path A — Control MD trajectories:** Use `00_extract_bstate_frames.sh` to sample random frames from long folded and/or unfolded production trajectories.
- **Path B — Prior WE simulation:** Use `00_extract_WE_intermediates.sh` to mine `west.h5` for segment endpoints that fall within a target region of phase space (e.g., intermediate MinDist values), useful for seeding new simulations with already-sampled intermediates.

Both paths converge at the same solvation and equilibration steps.

**Full pipeline:**

```
Source trajectories
        |
        v
[00] Extract frames / intermediates   →  00_bstate_pdbs/  (stripped RNA PDBs)
        |
        v
[01] Solvate + HRM (solvate.sh)       →  01_bstate_setup/ (solvated _HRM topologies)
        |
        v
[02] Minimize + equilibrate           →  02_bstate_equilibration/ (min/eq outputs)
        |
        v
[02] Extract final frames             →  03_bstates_final/ (struct_X/, bstates.txt)
        |
        v
[04] Analyze progress coordinates     →  04_analyze/  (pcoord table, figures,
        |                                              seed-set recommendations)
        v
   recommended bstate subsets for the next WE run(s)
```

---

## Prerequisites

| Tool | Version | Notes |
|------|---------|-------|
| AMBER | 24.0 | `cpptraj`, `pmemd`, `pmemd.cuda`, `parmed`, `tleap` |
| Python | 3.x | `h5py` for Path B; `numpy` + `matplotlib` for Step 4 analysis (both ship with the AMBER miniconda) |
| WESTPA | 2022.10 | Needed only if Path B reads a live simulation |
| `solvate.sh` | — | Dan Roe's iterative solvation script; must be present as `./solvate.sh` |
| SGE scheduler | — | Notre Dame CRC; `qsub`, `qstat`, `qdel` |

Load the environment before running any script:
```bash
module load amber/24.0
```

---

## Step 0A — Extract frames from control MD trajectories

**Script:** `00_extract_bstate_frames.sh`

Use this when your bstates come from long standard MD runs of the folded and/or unfolded states.

**Edit the configuration block at the top of the script:**

```bash
FOLDED_PRMTOP   # path to folded topology (.prmtop)
FOLDED_TRAJ     # path to folded production trajectory (.nc)
FOLDED_TOTAL_FRAMES=1000000

UNFOLDED_PRMTOP # path to unfolded topology
UNFOLDED_TRAJ   # path to unfolded trajectory
UNFOLDED_TOTAL_FRAMES=1000000

RNA_MASK=":1-14"   # cpptraj mask for RNA residues only
SKIP_FRAC=0.2      # skip first 20% of trajectory as equilibration buffer
```

**Run:**
```bash
bash 00_extract_bstate_frames.sh
```

**What it does:**

1. Uses Python to randomly select 10 frame indices from each trajectory, skipping the first 20% as equilibration buffer. Writes frame lists to `00_bstate_pdbs/folded_frames.txt` and `00_bstate_pdbs/unfolded_frames.txt`.
2. Builds and runs a `cpptraj` input that loads all 10 frames in a single pass, strips everything outside `RNA_MASK`, and writes one PDB per frame.

**Output:**
```
00_bstate_pdbs/
├── folded_frames.txt
├── unfolded_frames.txt
├── extract_folded.cpptraj
├── extract_unfolded.cpptraj
├── folded/
│   ├── folded_01_frame<N>.pdb
│   ├── folded_02_frame<N>.pdb
│   └── ... (10 PDBs)
└── unfolded/
    ├── unfolded_01_frame<N>.pdb
    └── ... (10 PDBs)
```

Each PDB contains only the stripped RNA molecule (no water, no ions). These feed into `01_prep_all_bstates.sh`.

---

## Step 0B — Extract intermediates from a prior WE simulation

**Script:** `00_extract_WE_intermediates.sh`

Use this when you want to seed a new simulation with structures sampled from a prior weighted ensemble run — for example, to populate an intermediate region of the progress coordinate that is otherwise hard to reach from the folded or unfolded endpoints.

**Usage:**
```bash
bash 00_extract_WE_intermediates.sh \
    --we-dir /path/to/westpa/simulation \
    --n 50 \
    --lo 5.0 \
    --hi 20.0 \
    --seed 42
```

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--we-dir` | *(required)* | Path to WESTPA simulation root (must contain `west.h5` and `traj_segs/`) |
| `--n` | 50 | Number of structures to extract |
| `--lo` | 5.0 | MinDist lower bound (Å) |
| `--hi` | 20.0 | MinDist upper bound (Å) |
| `--seed` | random | Integer seed for reproducibility |
| `--out` | `extracted_bstates_<dirname>` | Override output directory name |

**What it does:**

1. **Phase 1 — Scan `west.h5`:** Reads every iteration's segment endpoint pcoord values with `h5py`. Collects all segments where `MinDist` (pcoord dimension 1) falls in `[lo, hi]`. Randomly samples N candidates.
2. **Phase 2 — Extract files:** For each selected segment, dereferences the `struct.prmtop` symlink in `traj_segs/<iter>/<seg>/` to copy the real topology, copies `seg.rst` as the endpoint coordinate, and writes `pcoordreturn.dat` with the `rmsd<TAB>mindist` values.
3. **Phase 3 — Strip and write PDB:** Runs cpptraj on each extracted structure to autoimage and strip solvent, producing a bare RNA PDB.
4. **Phase 4 — Write manifests:** Generates `manifest.txt` (tracks iter/seg/pcoord for each structure) and `bstates.txt` (WESTPA basis state list with equal weights).

**Output:**
```
extracted_bstates_<WE_dirname>/
├── manifest.txt              (idx  iter  seg  rmsd  mindist)
├── bstates.txt               (WESTPA bstates format, equal weights)
└── StructureFiles/
    ├── struct_01/
    │   ├── struct.prmtop     (full solvated topology, symlink dereferenced)
    │   ├── struct.rst        (segment endpoint coordinates)
    │   ├── struct.pdb        (stripped RNA, residues 1-14)
    │   └── pcoordreturn.dat  (rmsd<TAB>mindist, one line)
    ├── struct_02/
    └── ...
```

**After this step,** copy the stripped PDBs into `00_bstate_pdbs/intermediate/` for the next step:
```bash
OUT=extracted_bstates_<WE_dirname>
mkdir -p 00_bstate_pdbs/intermediate
for d in $OUT/StructureFiles/struct_*/; do
    n=$(basename $d | sed "s/struct_//")
    cp $d/struct.pdb 00_bstate_pdbs/intermediate/intermediate_${n}.pdb
done
```

Then in `01_prep_all_bstates.sh`, set `RUN_INTERMEDIATE=true` and `RUN_FOLDED=false RUN_UNFOLDED=false`.

---

## Step 1 — Solvate and apply HRM

**Script:** `01_prep_all_bstates.sh`
**Dependency:** `solvate.sh` must be present in this directory. `module load amber/24.0` must have been run.

This script loops over all PDBs in `00_bstate_pdbs/` and for each one:
1. Calls `solvate.sh` (Dan Roe's iterative solvation tool) to build a truncated octahedral OPC water box with exactly `TARGET_WATERS` waters, then adds neutralizing K⁺ and 150 mM excess KCl.
2. Applies the Hydrogen Reparametrization Method (HRM) using ParmEd to adjust Lennard-Jones cross-terms for RNA aromatic hydrogens against backbone oxygens.

**Configuration (edit at top of script):**

```bash
RUN_FOLDED=true            # set false to skip
RUN_UNFOLDED=true
RUN_INTERMEDIATE=false

TARGET_WATERS=23000        # target water count (before ion swap)
INITIAL_BUFFER_FOLDED=38.0     # starting buffer guess (Å) for folded
INITIAL_BUFFER_UNFOLDED=12.0   # starting buffer guess for unfolded
INITIAL_BUFFER_INTERMEDIATE=20.0

TOLERANCE=5                # max waters over target before trimming

N_NEUTRALIZING=13          # K+ to neutralize RNA charge
N_EXCESS_IONS=61           # excess K+ and Cl- each (150 mM with 23000 waters)

PDB_DIR="00_bstate_pdbs"
OUT_DIR="01_bstate_setup"
SOLVATE_SH="./solvate.sh"
```

> **Buffer guesses:** `solvate.sh` iteratively adjusts the buffer and re-runs tleap until the water count is within `TOLERANCE`. A good initial guess saves iterations. Folded RNA is more compact and needs a larger buffer; unfolded is more extended and needs less. If the script takes many iterations (>20), adjust `INITIAL_BUFFER_*` toward the converged value reported in the summary table.

**Run:**
```bash
# This script has SGE headers — submit as a job on the cluster:
qsub 01_prep_all_bstates.sh

# Or run interactively on a GPU node:
bash 01_prep_all_bstates.sh
```

**What solvate.sh does internally:**

`solvate.sh` is a binary-search loop around tleap's `solvateOct`. It reads a config file (`solvate.in`) and iteratively adjusts the buffer value, running tleap and counting solvent molecules after each attempt. When it lands within tolerance, it trims or accepts the box. The shared `ionsin.ions` file is appended to the tleap session after solvation to add ions.

**HRM ParmEd modifications applied:**

The HRM modifies Lennard-Jones parameters to correct overestimated stacking interactions in OL3:
- Pyrimidine H6 × backbone O (OS/O2): r_min = 2.7 Å, ε = 0.050498 kcal/mol
- Purine H8 × backbone O: same
- Pyrimidine H5 × backbone O: same
- Sugar H1' × backbone O: r_min = 2.7 Å, ε = 0.051662 kcal/mol

**Output per bstate:**
```
01_bstate_setup/<bstate_name>/
├── solvate.in               (solvate.sh config)
├── solvate.log              (solvate.sh iteration log)
├── leapin.ff                (force field load commands for tleap)
├── ionsin.ions              (ion addition commands)
├── <bstate_name>_OL3.prmtop (solvated, OL3+OPC, with ions)
├── <bstate_name>_OL3.rst7
├── parmed_hrm.in            (parmed input for HRM)
├── parmed.log
├── <bstate_name>_HRM.parm7  (HRM-corrected topology, alias)
├── <bstate_name>_HRM.inpcrd (HRM-corrected coordinates, alias)
├── <bstate_name>_HRM.prmtop (HRM topology — used by step 02)
├── <bstate_name>_HRM.rst7   (HRM coordinates — used by step 02)
└── <bstate_name>_HRM.pdb
```

**Verification — read the summary table printed at the end:**
- All water counts should equal `TARGET_WATERS`
- All box lines should show three identical lengths and identical angles (truncated octahedron invariant: a = b = c and α = β = γ = 109.47°)
- All HRM column entries should read `OK`

---

## Step 2 — Minimize and equilibrate

**Script:** `02_03_equilibrate_extract_bstates.sh`

Runs a 10-step AMBER minimization and equilibration protocol on all bstates in parallel using SGE array jobs. Each step is one array job; each bstate is one array task. Steps are chained with `-hold_jid` so they execute sequentially but all bstates within a step run simultaneously.

**Protocol:**

| Step | Type | Length | Restraints | Notes |
|------|------|--------|------------|-------|
| min1 | Minimization | 1000 steps | 500 kcal/mol/Å² on RNA | Solvent relaxes |
| min2 | Minimization | 2500 steps | None | Full system |
| md1 | NVT heating | 100 ps | 25 kcal RNA | 0 → 300 K |
| md2a | NPT | 50 ps | 25 kcal RNA | Restraint release begins |
| md2b | NPT | 50 ps | 20 kcal RNA | |
| md2c | NPT | 50 ps | 15 kcal RNA | |
| md2d | NPT | 50 ps | 10 kcal RNA | |
| md2e | NPT | 50 ps | 5 kcal RNA | |
| md3 | NPT | 200 ps | None | Unrestrained; box volume correction runs after |
| md4 | NVT | 1 ns | None | Production-like; final frame extracted for bstates |

Minimization steps use CPU `pmemd` (CUDA does not support `imin=1`). All MD steps use `pmemd.cuda`.

**The postmd3 box correction:** After md3, `postmd3_calcboxlength.py` reads the average `VOLUME` from `md3.out` and the current box lengths from `md3.rst`. For a truncated octahedron, V = (1/2)(4/3)^(3/2) L³. The script computes the corrected L from the time-averaged volume and writes `md3_NewVolume.rst`, which becomes the input to md4. This ensures md4 starts from a box consistent with the ensemble-averaged density rather than an instantaneous snapshot.

**Configuration:**

```bash
BSTATE_SETUP_DIR="./01_bstate_setup"   # output from step 01
OUTPUT_DIR="./02_bstate_equilibration"
RESTRAINT_MASK=':1-14'                  # AMBER mask for RNA solute
START_STEP="min1"   # resume from a later step if earlier steps already ran
DRY_RUN="false"     # set true to generate scripts without submitting
```

**Set the bstate list:**

```bash
BSTATES_TO_RUN=(
    "folded_01" "folded_02" ... "folded_10"
    "unfolded_01" ... "unfolded_10"
    # "intermediate_01" ... "intermediate_20"
)
```

Only names listed here will be processed. Names must match directory names under `BSTATE_SETUP_DIR/`.

**Run:**
```bash
bash 02_03_equilibrate_extract_bstates.sh
```

> Run this on a login node — it only submits jobs, does not run MD itself. The script generates all AMBER `.in` files and SGE scripts locally, then calls `qsub` once per step.

**Dry run first:**

```bash
# In the script, set DRY_RUN="true", then:
bash 02_03_equilibrate_extract_bstates.sh
# Inspect generated scripts in 02_bstate_equilibration/run_*.sh
# When satisfied, set DRY_RUN="false" and re-run
```

**What the script does in three phases:**

1. **Phase 1 (local):** Creates `02_bstate_equilibration/<bstate_id>/` for each bstate. Copies the `_HRM.prmtop` and `_HRM.rst7`. Generates all 10 AMBER `.in` files and `postmd3_calcboxlength.py`. Writes `bstate_map.txt` (maps SGE task IDs to bstate directories).

2. **Phase 2 (submitted):** Submits 10 SGE array jobs, one per step, chained with `-hold_jid`. Each array task reads `bstate_map.txt` using `$SGE_TASK_ID` to find its working directory and runs the appropriate AMBER command.

3. **Phase 3 (submitted):** Submits a final extraction array that holds on the md4 array. Each task runs cpptraj on `md4.nc` to extract the last frame as `struct.ncrst`, `struct.rst`, `struct.rst7`, and `struct.pdb`. Writes `03_bstates_final/bstates.txt` and `03_bstates_final/struct_mapping.log`.

**Monitor jobs:**
```bash
qstat -u $USER
```

**Kill all submitted jobs if something goes wrong:**
```bash
qdel $(qstat -u $USER | grep eq_ | awk '{print $1}' | sort -u)
```

**Resuming from a later step:** If steps through md2e already completed but md3 failed, set `START_STEP="md3"` and re-run. The script will skip generating and submitting everything before md3.

---

## Step 4 — Analyze basis-state progress coordinates and recommend seed sets

**Scripts (run from the `bstates/` repo root):** `04_compute_raw_distances.sh`, then `04_analyze_bstate_pcoords.py`
**Outputs:** `04_analyze/raw/` (per-bstate raw distances) and `04_analyze/results/` (table, figures, log, recommendations)

Once the final basis states exist, this stage measures every progress coordinate for each `struct_X`, characterizes how diverse the seed set is, and recommends concrete bstate subsets for the next WE run(s). The pcoord definitions mirror `04_analyze/raw/` cpptraj inputs (`1_we_pcoord_distances.in` / `get_pcoord.cpptraj`) and `compute_Q.1.py` so the bstate values sit on the same footing as the live WE simulation's pcoords.

### Step 4a — Compute raw distances

```bash
module load amber
bash 04_compute_raw_distances.sh
```

For each `03_bstates_final/StructureFiles/struct_X/`, this runs `cpptraj` against the 2KOC NMR reference (`04_analyze/raw/2KOCFolded_NMR.{prmtop,rst7}`) and writes, into `04_analyze/raw/struct_X/`:

| File | Contents |
|------|----------|
| `pcoord_candidates.dat` | RMSD to NMR (global / stem / loop), G9 χ dihedral, Rg, end-to-end distance |
| `loop_hbonds_raw.dat` | 6 UUCG tetraloop H-bond distances |
| `stem_hbonds_raw.dat` | 14 Watson–Crick stem H-bond distances |
| `mindist.dat` | G1–C14 native-contact minimum distance (WE pcoord dimension 2) |

### Step 4b — Build the table, plots, and recommendations

```bash
python3 04_analyze_bstate_pcoords.py
```

This assembles a one-row-per-bstate table (RMSDs, MinDist, Q_stem/Q_loop via the `compute_Q.1.py` sigmoid + hard cutoffs, stem base-pair count, Rg, end-to-end, χ_G9), then assigns each bstate to a basin and selects seed subsets.

**Selection is driven only by the measured progress coordinates, never by the `folded_*` / `intermediate_*` / `unfolded_*` directory names.** Basins come from k-means(2) in standardized pcoord space; the transition (barrier) set is defined as a partially-formed stem (`1 ≤ n_stem_bp ≤ 4`) held at an intermediate G1–C14 MinDist. Diverse subsets within each basin are chosen by farthest-point sampling, since within-basin spread is the lever on higher-mode excitation.

**Recommended schemes** (grounded in `README` + `WEeDS_Background_reorganized.md` §4.6, §4.10, §4.11, §4.14, §4.15.1):

- **A — bidirectional 50/50:** diverse folded-side + unfolded-side seeds (the default coverage run).
- **B — bidirectional + transition seeds:** A plus partially-formed-stem seeds to sharpen ψ₂'s node and the barrier-region FES.
- **C — unidirectional unfolded perturbation:** diverse far-from-π unfolded seeds; strengthens ψ₂ in the pooled corpus. Pool with A/B and check the break-even criterion η_i > 1/√N_runs.
- **D — intermediate/barrier-only run:** analyzed, not blindly recommended — it sits on the ψ₂ node (c₂≈0), so it dilutes ψ₂ if pooled equal-weight, but is useful as a down-weighted higher-mode / barrier-resolution probe alongside a both-basins run. See the log for the full excitation-coefficient / pooling derivation.

**Outputs (`04_analyze/results/`):**

| File | Contents |
|------|----------|
| `bstate_pcoords.csv` / `.dat` | full pcoord table, one row per bstate, with data-driven basin and fold fraction |
| `recommendations.csv` | recommended `struct_X` sets per scheme (with role, basin, fold fraction) |
| `analysis.log` | summary statistics, label-vs-pcoord cross-tab, seed-cloud diversity, and the scheme recommendations with their mathematical grounding |
| `fig_pcoord_RMSD_MinDist.png` | the WE 2D pcoord space, colored by data-driven basin |
| `fig_Qstem_Qloop.png` | native-contact space (loop stays formed while stem melts) |
| `fig_pcoord_distributions.png` | per-coordinate histograms across all bstates |
| `fig_recommended_seeds.png` | Scheme B seeds marked on the WE pcoord |
| `fig_Rg_e2e.png` | global compaction, sized by stem base-pair count |

**Verification:**
```bash
ls 04_analyze/results/bstate_pcoords.csv     # 70 data rows + header
column -t -s, 04_analyze/results/recommendations.csv | head
less 04_analyze/results/analysis.log
```

---

## Output directory tree (complete finished workflow)

```
prep_bstates/
├── 00_bstate_pdbs/                        # Step 0A/0B: stripped RNA PDBs
│   ├── folded_frames.txt                  # randomly selected frame numbers
│   ├── unfolded_frames.txt
│   ├── extract_folded.cpptraj             # cpptraj input (generated)
│   ├── extract_unfolded.cpptraj
│   ├── folded/
│   │   ├── folded_01_frame<N>.pdb         # RNA-only, no solvent
│   │   ├── folded_02_frame<N>.pdb
│   │   └── ... (10 PDBs)
│   ├── unfolded/
│   │   ├── unfolded_01_frame<N>.pdb
│   │   └── ... (10 PDBs)
│   └── intermediate/                      # only if using Path B
│       ├── intermediate_01.pdb
│       └── ... (N PDBs)
│
├── extracted_bstates_<WE_dirname>/        # Step 0B only: WE intermediate extraction
│   ├── manifest.txt                       # idx  iter  seg  rmsd  mindist
│   ├── bstates.txt                        # WESTPA basis state list (pre-equilibration)
│   └── StructureFiles/
│       ├── struct_01/
│       │   ├── struct.prmtop              # full solvated topology (dereferenced)
│       │   ├── struct.rst                 # segment endpoint coordinates
│       │   ├── struct.pdb                 # stripped RNA PDB
│       │   └── pcoordreturn.dat           # rmsd<TAB>mindist
│       └── ... (N struct_XX directories)
│
├── 01_bstate_setup/                       # Step 1: solvated + HRM topologies
│   ├── leapin.ff                          # shared force field load commands
│   ├── ionsin.ions                        # shared ion addition commands
│   ├── folded_01/
│   │   ├── solvate.in
│   │   ├── solvate.log
│   │   ├── folded_01_OL3.prmtop
│   │   ├── folded_01_OL3.rst7
│   │   ├── parmed_hrm.in
│   │   ├── parmed.log
│   │   ├── folded_01_HRM.prmtop           # → input to step 02
│   │   ├── folded_01_HRM.rst7             # → input to step 02
│   │   └── folded_01_HRM.pdb
│   ├── folded_02/
│   │   └── ...
│   ├── ... (one directory per bstate)
│   ├── unfolded_01/
│   │   └── ...
│   └── intermediate_01/                   # only if RUN_INTERMEDIATE=true
│       └── ...
│
├── 02_bstate_equilibration/               # Step 2: min/eq working directories
│   ├── bstate_map.txt                     # task_id  bstate_id  workdir (SGE mapping)
│   ├── 02_equilibrate_bstates.log
│   ├── run_min1.sh                        # SGE array script — min1
│   ├── run_min2.sh
│   ├── run_md1.sh
│   ├── run_md2a.sh
│   ├── run_md2b.sh
│   ├── run_md2c.sh
│   ├── run_md2d.sh
│   ├── run_md2e.sh
│   ├── run_md3.sh
│   ├── run_md4.sh
│   ├── run_extract_final.sh               # SGE array script — final frame extraction
│   ├── logs/
│   │   ├── min1_1.log                     # stdout for task 1 of min1 array
│   │   ├── min1_2.log
│   │   ├── ... (one log per step per bstate)
│   │   └── extract_20.log
│   ├── folded_01/
│   │   ├── folded_01_HRM.prmtop           # copied from 01_bstate_setup
│   │   ├── folded_01_HRM.rst7
│   │   ├── min1.in                        # generated AMBER input
│   │   ├── min1.out                       # AMBER output
│   │   ├── min1.rst
│   │   ├── min2.in
│   │   ├── min2.out
│   │   ├── min2.rst
│   │   ├── md1.in
│   │   ├── md1.out
│   │   ├── md1.rst
│   │   ├── md1.nc
│   │   ├── md2a.in  (through md2e.in)
│   │   ├── md2a.out (through md2e.out)
│   │   ├── md2a.rst (through md2e.rst)
│   │   ├── md2a.nc  (through md2e.nc)
│   │   ├── md3.in
│   │   ├── md3.out
│   │   ├── md3.rst
│   │   ├── md3.nc
│   │   ├── postmd3_calcboxlength.py       # box volume correction script
│   │   ├── md3_NewVolume.rst              # corrected box dimensions → md4 input
│   │   ├── md4.in
│   │   ├── md4.out
│   │   ├── md4.rst
│   │   └── md4.nc
│   ├── folded_02/
│   │   └── ...
│   └── ... (one directory per bstate)
│
└── 03_bstates_final/                      # Step 2 Phase 3: WESTPA-ready basis states
    ├── bstates.txt                        # WESTPA basis state list (equal weights)
    ├── struct_mapping.log                 # struct index → bstate ID → source paths
    ├── struct_0/
    │   ├── struct.prmtop                  # topology (copied from eq directory)
    │   ├── struct.ncrst                   # last md4 frame — NetCDF restart
    │   ├── struct.rst                     # last md4 frame — ASCII restart
    │   ├── struct.rst7                    # last md4 frame — NetCDF restart (alias)
    │   ├── struct.pdb                     # last md4 frame — PDB
    │   └── extract.cpptraj                # cpptraj script used for extraction
    ├── struct_1/
    │   └── ...
    └── ... (one directory per bstate, 0-indexed)
│
├── 04_compute_raw_distances.sh           # Step 4a: per-bstate cpptraj driver
├── 04_analyze_bstate_pcoords.py          # Step 4b: table + figures + recommendations
└── 04_analyze/                           # Step 4: progress-coordinate analysis
    ├── raw/                              # Step 4a output
    │   ├── 2KOCFolded_NMR.prmtop         # folded-state RMSD reference
    │   ├── 2KOCFolded_NMR.rst7
    │   ├── struct_0/
    │   │   ├── get_pcoord.cpptraj        # generated cpptraj input
    │   │   ├── pcoord_candidates.dat     # RMSDs, chi_G9, Rg, end-to-end
    │   │   ├── loop_hbonds_raw.dat       # 6 loop H-bond distances
    │   │   ├── stem_hbonds_raw.dat       # 14 stem H-bond distances
    │   │   ├── mindist.dat               # G1–C14 native-contact MinDist
    │   │   └── cpptraj.log
    │   └── ... (one directory per bstate)
    └── results/                          # Step 4b output
        ├── bstate_pcoords.csv            # full pcoord table (+ basin, fold_frac)
        ├── bstate_pcoords.dat            # same table, whitespace-aligned
        ├── recommendations.csv           # recommended struct sets per scheme
        ├── analysis.log                  # stats, diversity, recommendations
        ├── fig_pcoord_RMSD_MinDist.png
        ├── fig_Qstem_Qloop.png
        ├── fig_pcoord_distributions.png
        ├── fig_recommended_seeds.png
        └── fig_Rg_e2e.png
```

---

## Verification checklist

After each step, confirm the following before proceeding.

**After step 0A (`00_extract_bstate_frames.sh`):**
```bash
ls 00_bstate_pdbs/folded/ | wc -l    # should be 10
ls 00_bstate_pdbs/unfolded/ | wc -l  # should be 10
# Check PDB files are non-empty RNA-only structures:
head -5 00_bstate_pdbs/folded/folded_01_frame*.pdb
```

**After step 0B (`00_extract_WE_intermediates.sh`):**
```bash
wc -l extracted_bstates_*/manifest.txt   # header + N rows
ls extracted_bstates_*/StructureFiles/struct_*/struct.pdb | wc -l  # should be N
```

**After step 1 (`01_prep_all_bstates.sh`):**

Read the summary table printed at the end of the run (also in `solvate.log`):
- Water count column: all entries should equal `TARGET_WATERS`
- Box column: all three lengths should be equal; angles should be 109.47°
- HRM column: all `OK`

```bash
# Quick check:
for d in 01_bstate_setup/*/; do
    echo -n "$(basename $d): "
    ls $d/*_HRM.prmtop 2>/dev/null | head -1 || echo "MISSING HRM PRMTOP"
done
```

**After step 2 (`02_03_equilibrate_extract_bstates.sh`):**
```bash
# Check md4 completed for all bstates:
ls 02_bstate_equilibration/*/md4.rst | wc -l   # should equal N bstates

# Check final structures exist:
ls 03_bstates_final/struct_*/struct.pdb | wc -l

# Check for errors in logs:
grep -l "ERROR" 02_bstate_equilibration/logs/*.log

# Check minimization convergence:
grep "FINAL RESULTS" 02_bstate_equilibration/*/min1.out | wc -l

# Check bstates.txt was written:
cat 03_bstates_final/bstates.txt
```

**After step 4 (`04_compute_raw_distances.sh` + `04_analyze_bstate_pcoords.py`):**
```bash
# One raw mindist file per bstate (should equal N bstates):
ls 04_analyze/raw/struct_*/mindist.dat | wc -l

# No cpptraj failures:
grep -l "Error" 04_analyze/raw/struct_*/cpptraj.log

# Table has a row per bstate (N + 1 with header) and all results exist:
wc -l 04_analyze/results/bstate_pcoords.csv
ls 04_analyze/results/*.png 04_analyze/results/*.csv 04_analyze/results/analysis.log

# Sanity-check the extremes: folded ~ low RMSD/MinDist, unfolded ~ high:
column -t -s, 04_analyze/results/bstate_pcoords.csv | head
```

---

## Using the final bstates in WESTPA

Copy `03_bstates_final/` into your WESTPA simulation directory as the `bstates/StructureFiles/` tree. The `bstates.txt` file uses relative paths of the form `StructureFiles/struct_X`, which matches the extraction script's output layout.

In `west.cfg`, point to `bstates.txt`:
```yaml
west:
  system:
    basis_states_file: bstates/bstates.txt
```

Each `struct_X/` directory needs at minimum `struct.ncrst` (or `struct.rst`) and `struct.prmtop`. The `runseg.sh` script for your simulation should be written to read these files by those names, consistent with how other bstates in the `Example2D_RMSDMinDist` template are structured.

To seed a *subset* of these bstates (e.g. a bidirectional or perturbation run rather than all 70), use the Step 4 recommendations: `04_analyze/results/recommendations.csv` lists the chosen `struct_X` per scheme, and `04_analyze/results/analysis.log` explains the rationale. Build a `bstates.txt` from the selected `struct_X` rows, renormalizing the weights so they sum to 1.

---

## Notes on solvate.sh

`solvate.sh` (by Daniel R. Roe) performs iterative buffer adjustment around tleap's `solvateOct`. It reads a plain-text config file and uses a binary-search-like update rule to converge on the target water count. Key parameters passed via `solvate.in`:

| Parameter | Meaning |
|-----------|---------|
| `target` | Exact water molecule count to reach |
| `buffer` | Initial box buffer guess in Å |
| `tol` | Waters over target that can be removed by trimming (not re-solvating) |
| `mode 0` | Truncated octahedral box (`solvateOct`) |
| `solventunit OPCBOX` | OPC water model |
| `molname m` | Name of the solute unit in tleap (must match `ionsin` file) |

Ion addition is handled by the shared `ionsin.ions` file, which is appended to the tleap session after solvation. The `addions m K+ 0` command neutralizes the system first, then `addionsrand` places excess salt at random positions.cd bs