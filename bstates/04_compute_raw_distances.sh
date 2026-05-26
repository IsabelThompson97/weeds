#!/bin/bash
# =============================================================================
# 01_compute_raw_distances.sh
# -----------------------------------------------------------------------------
# Compute all raw progress-coordinate distances for every final basis state
# in 03_bstates_final/StructureFiles/struct_*/ .
#
# Run from the bstates repo root. For each struct_X it runs cpptraj with the
# 2KOC NMR reference and writes per-bstate raw distance files into
# 04_analyze/raw/struct_X/:
#     pcoord_candidates.dat   rmsd_global, rmsd_stem, rmsd_loop, chi_G9, Rg, d_e2e
#     loop_hbonds_raw.dat     6 UUCG loop H-bond distances
#     stem_hbonds_raw.dat     14 WC stem H-bond distances
#     mindist.dat             nativecontacts min distance between G1 and C14
#
# The distance/mask definitions mirror 1_we_pcoord_distances.in and
# get_pcoord.cpptraj exactly so the bstate pcoords are on the same footing
# as the WE simulation pcoords.
#
# Post-process with 02_analyze_bstate_pcoords.py.
#
# Usage:  module load amber   (or amber/24.0);  bash 01_compute_raw_distances.sh
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"     # .../bstates (script lives at repo root)
BSTATE_DIR="$HERE/03_bstates_final/StructureFiles"
REF_PRM="$HERE/04_analyze/raw/2KOCFolded_NMR.prmtop"
REF_RST="$HERE/04_analyze/raw/2KOCFolded_NMR.rst7"
OUT_ROOT="$HERE/04_analyze/raw"

command -v cpptraj >/dev/null 2>&1 || { echo "ERROR: cpptraj not found. Run 'module load amber' first." >&2; exit 1; }
[[ -f "$REF_PRM" && -f "$REF_RST" ]] || { echo "ERROR: NMR reference files missing in $HERE" >&2; exit 1; }

mkdir -p "$OUT_ROOT"

shopt -s nullglob
structs=("$BSTATE_DIR"/struct_*)
echo "Found ${#structs[@]} bstate directories under $BSTATE_DIR"

n=0
for sdir in "${structs[@]}"; do
    sid="$(basename "$sdir")"
    prm="$sdir/struct.prmtop"
    rst="$sdir/struct.rst7"
    [[ -f "$prm" && -f "$rst" ]] || { echo "  skip $sid (missing struct.prmtop/struct.rst7)"; continue; }

    odir="$OUT_ROOT/$sid"
    mkdir -p "$odir"
    cin="$odir/get_pcoord.cpptraj"

    cat > "$cin" <<EOF
parm $prm
trajin $rst

parm $REF_PRM [NMR]
reference $REF_RST parm [NMR]

autoimage :1-14

# --- RMSD (to NMR), G9 chi, Rg, end-to-end -> pcoord_candidates.dat ---
rmsd rmsd_global :1-14&!@H=      reference :1-14&!@H=      out $odir/pcoord_candidates.dat
rmsd rmsd_stem   :1-5,10-14&!@H= reference :1-5,10-14&!@H= out $odir/pcoord_candidates.dat
rmsd rmsd_loop   :6-9&!@H=       reference :6-9&!@H=       out $odir/pcoord_candidates.dat
dihedral chi_G9  :9@O4' :9@C1' :9@N9 :9@C4                 out $odir/pcoord_candidates.dat
radgyr rog :1-14 nomax                                    out $odir/pcoord_candidates.dat
distance d_e2e   :1@O6 :14@O3'                            out $odir/pcoord_candidates.dat

# --- 6 loop H-bond distances -> loop_hbonds_raw.dat ---
distance d1_G9N1_U6O2   :9@N1  :6@O2   out $odir/loop_hbonds_raw.dat
distance d2_G9N2_U6O2   :9@N2  :6@O2   out $odir/loop_hbonds_raw.dat
distance d3_U6O2p_G9O6  :6@O2' :9@O6   out $odir/loop_hbonds_raw.dat
distance d4_U7O2p_G9N7  :7@O2' :9@N7   out $odir/loop_hbonds_raw.dat
distance d5_U7O2p_G9O6  :7@O2' :9@O6   out $odir/loop_hbonds_raw.dat
distance d6_C8N4_U7OP2  :8@N4  :7@OP2  out $odir/loop_hbonds_raw.dat

# --- 14 stem WC H-bond distances -> stem_hbonds_raw.dat ---
distance s1_G1N1_C14N3  :1@N1  :14@N3  out $odir/stem_hbonds_raw.dat
distance s2_G1N2_C14O2  :1@N2  :14@O2  out $odir/stem_hbonds_raw.dat
distance s3_G1O6_C14N4  :1@O6  :14@N4  out $odir/stem_hbonds_raw.dat
distance s4_G2N1_C13N3  :2@N1  :13@N3  out $odir/stem_hbonds_raw.dat
distance s5_G2N2_C13O2  :2@N2  :13@O2  out $odir/stem_hbonds_raw.dat
distance s6_G2O6_C13N4  :2@O6  :13@N4  out $odir/stem_hbonds_raw.dat
distance s7_C3N3_G12N1  :3@N3  :12@N1  out $odir/stem_hbonds_raw.dat
distance s8_C3O2_G12N2  :3@O2  :12@N2  out $odir/stem_hbonds_raw.dat
distance s9_C3N4_G12O6  :3@N4  :12@O6  out $odir/stem_hbonds_raw.dat
distance s10_A4N1_U11N3 :4@N1  :11@N3  out $odir/stem_hbonds_raw.dat
distance s11_A4N6_U11O4 :4@N6  :11@O4  out $odir/stem_hbonds_raw.dat
distance s12_C5N3_G10N1 :5@N3  :10@N1  out $odir/stem_hbonds_raw.dat
distance s13_C5O2_G10N2 :5@O2  :10@N2  out $odir/stem_hbonds_raw.dat
distance s14_C5N4_G10O6 :5@N4  :10@O6  out $odir/stem_hbonds_raw.dat

# --- native-contact min distance between terminal pair G1 / C14 (WE pcoord dim) ---
nativecontacts :1 :14 mindist name nc
run
write $odir/mindist.dat nc[mindist]
quit
EOF

    cpptraj -i "$cin" > "$odir/cpptraj.log" 2>&1 || { echo "  ERROR cpptraj failed for $sid (see $odir/cpptraj.log)"; continue; }
    n=$((n+1))
done

echo "Done: wrote raw distances for $n bstates into $OUT_ROOT"
