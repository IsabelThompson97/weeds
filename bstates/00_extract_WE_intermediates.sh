#!/bin/bash
# =============================================================================
# 00_extract_WE_intermediates.sh
#
# Extract N random trajectory segment endpoints from a WESTPA west.h5 that
# fall within a MinDist window, package them as stripped RNA PDBs + full
# solvated prmtop/rst, and write a manifest and bstates.txt.
#
# Run this from the bstates/ directory. Output lands in a folder named after
# the WE simulation so you always know where structures came from.
#
# OUTPUT: extracted_bstates_<WE_dirname>/
#   ├── manifest.txt               (idx  iter  seg  rmsd  mindist)
#   ├── bstates.txt                (WESTPA bstates manifest, equal weights)
#   └── StructureFiles/
#       └── struct_NN/
#           ├── struct.prmtop      (real file — symlink chain dereferenced)
#           ├── struct.rst         (seg.rst endpoint coordinates)
#           ├── struct.pdb         (stripped RNA :1-14, no solvent)
#           └── pcoordreturn.dat   (rmsd<TAB>mindist — single line)
#
# USAGE:
#   bash 00_extract_WE_intermediates.sh --we-dir /path/to/WE/sim [options]
#
# OPTIONS:
#   --we-dir  PATH   Path to WESTPA simulation root (required)
#   --n       N      Number of structures to extract (default: 50)
#   --lo      F      MinDist lower bound in Å (default: 5.0)
#   --hi      F      MinDist upper bound in Å (default: 20.0)
#   --seed    N      Random seed for reproducibility (default: random)
#   --out     PATH   Override output directory name
#
# NEXT STEPS after this script finishes:
#   Copy the stripped PDBs into 00_bstate_pdbs/<subfolder>/:
#     OUT=extracted_bstates_<WE_dirname>
#     mkdir -p 00_bstate_pdbs/intermediate
#     for d in $OUT/StructureFiles/struct_*/; do
#         n=$(basename $d | sed "s/struct_//")   # gives 01, 02, ...
#         cp $d/pdb4amber_struct.pdb 00_bstate_pdbs/intermediate/intermediate_${n}.pdb
#     done
#   Then set RUN_INTERMEDIATE=true in 01_prep_all_bstates.sh and run it.
# =============================================================================

# module load amber

set -euo pipefail

# ---------- defaults ----------------------------------------------------------
WE_SIM_DIR=""
N_STRUCTS=50
MINDIST_LO=5.0
MINDIST_HI=20.0
SEED=""
OUT_DIR=""   # set after parsing args, once WE_SIM_DIR is known

# ---------- argument parsing --------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --we-dir) WE_SIM_DIR="$2";  shift 2 ;;
        --n)      N_STRUCTS="$2";   shift 2 ;;
        --lo)     MINDIST_LO="$2";  shift 2 ;;
        --hi)     MINDIST_HI="$2";  shift 2 ;;
        --seed)   SEED="$2";        shift 2 ;;
        --out)    OUT_DIR="$2";     shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$WE_SIM_DIR" ]]; then
    echo "ERROR: --we-dir is required." >&2
    echo "  Usage: bash 00_extract_WE_intermediates.sh --we-dir /path/to/WE/sim" >&2
    exit 1
fi

WE_SIM_DIR=$(realpath "$WE_SIM_DIR")
WE_NAME=$(basename "$WE_SIM_DIR")
H5FILE="$WE_SIM_DIR/west.h5"
TRAJ_SEGS="$WE_SIM_DIR/traj_segs"

if [[ -z "$OUT_DIR" ]]; then
    OUT_DIR="extracted_bstates_${WE_NAME}"
fi

# ---------- environment checks ------------------------------------------------
if [[ -z "${AMBERHOME:-}" ]]; then
    module load westpa/2022.10 cuda/10.2
    module load mpich/4.3.2/gcc/11.5.0
    export AMBERHOME=/opt/crc/a/amber/24.0
    export PATH=$AMBERHOME/bin:$PATH
    export LD_LIBRARY_PATH=$AMBERHOME/lib:$LD_LIBRARY_PATH
fi
CPPTRAJ_BIN="$AMBERHOME/bin/cpptraj"

if [[ ! -f "$CPPTRAJ_BIN" ]]; then
    echo "ERROR: cpptraj not found at $CPPTRAJ_BIN" >&2
    exit 1
fi

if [[ ! -f "$H5FILE" ]]; then
    echo "ERROR: west.h5 not found at $H5FILE" >&2
    exit 1
fi

if [[ -d "$OUT_DIR" ]]; then
    echo "WARNING: $OUT_DIR already exists. Contents may be overwritten."
fi

echo "=== WE simulation : $WE_SIM_DIR"
echo "=== Output        : $OUT_DIR"
echo "=== MinDist range : [${MINDIST_LO}, ${MINDIST_HI}] Å  |  N=${N_STRUCTS}"
echo ""

# ---------- Phase 1: scan west.h5 for candidates via Python -------------------
echo "=== Phase 1: Scanning west.h5 ==="

CANDIDATE_FILE="/tmp/we_candidates_$$.txt"

python3 << PYEOF
import h5py
import random
import sys

h5file    = "$H5FILE"
n_structs = $N_STRUCTS
lo        = $MINDIST_LO
hi        = $MINDIST_HI
seed_str  = "$SEED"
outfile   = "$CANDIDATE_FILE"

candidates = []

with h5py.File(h5file, "r") as h5:
    iter_group = h5["iterations"]
    iter_keys  = sorted(iter_group.keys())
    print(f"  Scanning {len(iter_keys)} iterations...", flush=True)

    for ik in iter_keys:
        pcoord = iter_group[ik]["pcoord"][:]  # (n_segs, n_frames, n_dim)
        endpoint = pcoord[:, -1, :]           # segment endpoint
        rmsd_vals    = endpoint[:, 0]
        mindist_vals = endpoint[:, 1]

        iter_num = int(ik.split("_")[1])
        for seg_idx, (rmsd, mindist) in enumerate(zip(rmsd_vals, mindist_vals)):
            if lo <= mindist <= hi:
                candidates.append((iter_num, seg_idx, float(rmsd), float(mindist)))

print(f"  Found {len(candidates)} candidate segments.", flush=True)

if len(candidates) < n_structs:
    print(f"ERROR: only {len(candidates)} candidates — fewer than requested {n_structs}.",
          file=sys.stderr)
    sys.exit(1)

if seed_str:
    random.seed(int(seed_str))
else:
    import os
    random.seed(int.from_bytes(os.urandom(4), "little"))

selected = random.sample(candidates, n_structs)
selected.sort(key=lambda x: (x[0], x[1]))

with open(outfile, "w") as f:
    for i, (itr, seg, rmsd, mindist) in enumerate(selected):
        f.write(f"{i}\t{itr}\t{seg}\t{rmsd:.4f}\t{mindist:.4f}\n")

print(f"  Selected {n_structs} structures.", flush=True)
PYEOF

echo ""

# ---------- Phase 2 & 3: extract files and generate stripped PDBs -------------
echo "=== Phase 2+3: Extracting structures and generating stripped PDBs ==="

mkdir -p "$OUT_DIR/StructureFiles"
MANIFEST="$OUT_DIR/manifest.txt"
echo -e "idx\titer\tseg\trmsd\tmindist" > "$MANIFEST"

while IFS=$'\t' read -r IDX ITER SEG RMSD MINDIST; do
    ITER_PAD=$(printf '%06d' "$ITER")
    SEG_PAD=$(printf '%06d' "$SEG")
    STRUCT_NAME=$(printf 'struct_%02d' "$(( IDX + 1 ))")
    SEG_DIR="${TRAJ_SEGS}/${ITER_PAD}/${SEG_PAD}"
    STRUCT_DIR="$OUT_DIR/StructureFiles/${STRUCT_NAME}"

    echo "  [$IDX] iter=$ITER_PAD seg=$SEG_PAD  RMSD=${RMSD} MinDist=${MINDIST}"

    if [[ ! -f "$SEG_DIR/seg.rst" ]]; then
        echo "    WARNING: $SEG_DIR/seg.rst not found — skipping." >&2
        continue
    fi
    if [[ ! -e "$SEG_DIR/struct.prmtop" ]]; then
        echo "    WARNING: $SEG_DIR/struct.prmtop not found — skipping." >&2
        continue
    fi

    mkdir -p "$STRUCT_DIR"

    # Dereference symlink chain to get the real prmtop
    cp -L "$SEG_DIR/struct.prmtop" "$STRUCT_DIR/struct.prmtop"
    cp    "$SEG_DIR/seg.rst"        "$STRUCT_DIR/struct.rst"

    printf "%.4f\t%.4f\n" "$RMSD" "$MINDIST" > "$STRUCT_DIR/pcoordreturn.dat"

    # Strip solvent and write PDB
    CPPTRAJ_IN=$(mktemp /tmp/strip_XXXXXX.cpptraj)
    cat > "$CPPTRAJ_IN" << CPPEOF
parm ${STRUCT_DIR}/struct.prmtop
trajin ${STRUCT_DIR}/struct.rst
autoimage :1-14
strip !:1-14
trajout ${STRUCT_DIR}/struct.pdb pdb
go
quit
CPPEOF

    CPPTRAJ_LOG="$STRUCT_DIR/strip.log"
    if "$CPPTRAJ_BIN" -i "$CPPTRAJ_IN" > "$CPPTRAJ_LOG" 2>&1 && [[ -s "$STRUCT_DIR/struct.pdb" ]]; then
        echo -e "${IDX}\t${ITER}\t${SEG}\t${RMSD}\t${MINDIST}" >> "$MANIFEST"
        rm -f "$CPPTRAJ_IN"
    else
        echo "    ERROR: cpptraj failed for $STRUCT_NAME — see $CPPTRAJ_LOG" >&2
        rm -f "$CPPTRAJ_IN"
        continue
    fi

done < "$CANDIDATE_FILE"

# ---------- Phase 4: generate bstates.txt ------------------------------------
echo ""
echo "=== Phase 4: Writing bstates.txt ==="

N_WRITTEN=$(( $(wc -l < "$MANIFEST") - 1 ))
PROB=$(python3 -c "print(f'{1.0/$N_WRITTEN:.20e}')")

BSTATES_TXT="$OUT_DIR/bstates.txt"
> "$BSTATES_TXT"
IDX_OUT=0
while IFS=$'\t' read -r IDX ITER SEG RMSD MINDIST; do
    STRUCT_NAME=$(printf 'struct_%02d' "$(( IDX + 1 ))")
    printf "%d    %s    StructureFiles/%s\n" "$IDX_OUT" "$PROB" "$STRUCT_NAME" >> "$BSTATES_TXT"
    (( IDX_OUT++ )) || true
done < <(tail -n +2 "$MANIFEST")

rm -f "$CANDIDATE_FILE"

# for d in ${OUT_DIR}/StructureFiles/struct_*/; do
#     pdb4amber -i ${d}/struct.pdb -o ${d}/pdb4amber_struct.pdb 2>&1 | grep -v "^$"
# done

# ---------- summary -----------------------------------------------------------
echo ""
echo "=== Done ==="
echo "  WE simulation    : $WE_SIM_DIR"
echo "  Output directory : $OUT_DIR/"
echo "  Structures written: $N_WRITTEN"
echo "  Manifest         : $OUT_DIR/manifest.txt"
echo "  bstates.txt      : $OUT_DIR/bstates.txt"
echo ""
echo "Next — copy stripped PDBs into 00_bstate_pdbs/:"
echo "  mkdir -p 00_bstate_pdbs/intermediate"
echo "  for d in ${OUT_DIR}/StructureFiles/struct_*/; do"
echo '      n=$(basename $d | sed "s/struct_/")'
echo "      cp \$d/struct.pdb 00_bstate_pdbs/intermediate/intermediate_\${n}.pdb"
echo "  done"
echo ""
echo "Then in 01_prep_all_bstates.sh set:"
echo "  RUN_FOLDED=false  RUN_UNFOLDED=false  RUN_INTERMEDIATE=true"
echo "  INITIAL_BUFFER_INTERMEDIATE=20.0"
echo "  bash 01_prep_all_bstates.sh"
