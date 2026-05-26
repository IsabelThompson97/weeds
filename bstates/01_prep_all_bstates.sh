#!/bin/bash
#$ -M ithomps3@nd.edu
#$ -m abe
#$ -pe smp 8
#$ -q gpu
#$ -l gpu_card=1
#$ -N solvate_bstate
# =============================================================================
# prep_all_bstates.sh
#
# Wrapper around Dan Roe's Solvate.sh to prepare 20 WE bstates with
# identical water counts in truncated octahedral boxes.
#
# PREREQUISITES:
#   - Solvate.sh in your PATH or in the current directory
#   - Stripped RNA PDBs in bstate_pdbs/folded/ and bstate_pdbs/unfolded/
#   - AMBER (tleap, cpptraj, parmed) loaded
#
# USAGE: Edit configuration below, then: bash prep_all_bstates.sh
# =============================================================================

# ---- CONFIGURATION ----

# Which categories to process (set false to skip already-completed categories)
RUN_FOLDED=true
RUN_UNFOLDED=true
RUN_INTERMEDIATE=true

TARGET_WATERS=23135           # Target water count (before ion addition)
INITIAL_BUFFER_FOLDED=38.0       # Starting buffer guess for folded structures
INITIAL_BUFFER_UNFOLDED=12.0     # Starting buffer guess for unfolded structures
INITIAL_BUFFER_INTERMEDIATE=20.0 # Starting buffer guess for intermediate structures
TOLERANCE=5                      # How many waters over target is acceptable

# Ion counts (from SLTCAP at 150 mM with ~23000 waters)
N_NEUTRALIZING=13           # K+ to neutralize RNA (-13 charge)
N_EXCESS_IONS=61            # Excess K+ and Cl- each (SLTCAP calculation for 150 mM salt with 23000 waters, 4.5 kDa, 0 charge)

# Paths
PDB_DIR="00_bstate_pdbs"
OUT_DIR="01_bstate_setup"  # Output directory for solvated systems
SOLVATE_SH="./01_solvate.sh"  # Path to Dan Roe's Solvate.sh

module load amber/24.0

# ---- VERIFY SOLVATE.SH EXISTS ----
if [ ! -f "${SOLVATE_SH}" ]; then
    echo "Error: 01_solvate.sh not found at ${SOLVATE_SH}"
    echo "Download from: https://github.com/drroe/Solvate.sh"
    exit 1
fi
chmod +x ${SOLVATE_SH}

mkdir -p ${OUT_DIR}

# =============================================================================
# Create the shared leapin file (force field loading)
# =============================================================================
cat > ${OUT_DIR}/leapin.ff << 'EOF'
source leaprc.RNA.OL3
source leaprc.water.opc
EOF

# =============================================================================
# Empty ionsin file: Solvate.sh runs water-only.
#
# Why: Solvate.sh appends its `remove m m.NNN` water-trimming commands AFTER
# the ionsin contents. When ionsin contains addions/addionsrand, the ions are
# placed before the removes, and the remove-target residue indices computed by
# Solvate.sh (LASTRES = SOLUTERES + solvent_res_added) no longer resolve to
# valid residues in the now-reorganized unit. tleap fails with
# "remove: Argument #2 is of type String".
#
# Fix: leave ionsin empty so Solvate.sh's remove branch operates on a
# contiguous water-only residue list. After Solvate.sh converges, a second
# tleap pass (below) inserts addions/addionsrand BEFORE saveamberparm in the
# converged temp.leap.in and re-runs. Final composition is identical to the
# intended setup; only the order of (remove, addions) in tleap is swapped.
# =============================================================================
cat > ${OUT_DIR}/ionsin.ions << 'EOF'
# Intentionally empty — ions are added in a separate tleap pass.
EOF

# =============================================================================
# Function to process one bstate
# =============================================================================
process_bstate() {
    local PDB_FILE=$1
    local BSTATE_NAME=$2
    local INIT_BUFFER=$3
    local WORKDIR="${OUT_DIR}/${BSTATE_NAME}"

    echo ""
    echo "============================================================"
    echo "  Processing: ${BSTATE_NAME}"
    echo "  PDB: ${PDB_FILE}"
    echo "  Initial buffer guess: ${INIT_BUFFER} A"
    echo "============================================================"

    mkdir -p ${WORKDIR}
    cd ${WORKDIR}

    # Copy shared files
    cp ../leapin.ff .
    cp ../ionsin.ions .

    # ----- Create Solvate.sh input config -----
    cat > solvate.in << EOF
target ${TARGET_WATERS}
buffer ${INIT_BUFFER}
pdb ${PDB_FILE}
top ${BSTATE_NAME}_OL3.prmtop
crd ${BSTATE_NAME}_OL3.rst7
leapin leapin.ff
ionsin ionsin.ions
mode 0
solventunit OPCBOX
molname m
tol ${TOLERANCE}
EOF

    # ----- Run Solvate.sh -----
    echo "  Running Solvate.sh..."
    bash "../../${SOLVATE_SH}" solvate.in > solvate.log 2>&1
    SOLVATE_STATUS=$?

    if [ ${SOLVATE_STATUS} -ne 0 ]; then
        echo "  *** ERROR: Solvate.sh failed. Check solvate.log ***"
        cd - > /dev/null
        return 1
    fi

    # Verify water-only output files exist
    if [ ! -s "${BSTATE_NAME}_OL3.prmtop" ] || [ ! -s "${BSTATE_NAME}_OL3.rst7" ]; then
        echo "  *** ERROR: Solvate.sh output files not created. Check solvate.log ***"
        cd - > /dev/null
        return 1
    fi

    # ----- Add ions in a second tleap pass -----
    # Solvate.sh's temp.leap.in contains the converged solvateoct + (optional)
    # remove commands needed to hit exactly TARGET_WATERS. We insert
    # addions/addionsrand BEFORE saveamberparm and re-run tleap, so removes
    # operate on pure water (correct indices) and ions are added afterward.
    echo "  Adding ions (second tleap pass)..."
    awk -v n_ions="${N_EXCESS_IONS}" '
        /^saveamberparm/ {
            print "addions m K+ 0"
            print "addionsrand m K+ " n_ions " Cl- " n_ions
        }
        { print }
    ' temp.leap.in > add_ions.leap.in

    tleap -f add_ions.leap.in > add_ions.leap.out 2>&1
    TLEAP_STATUS=$?

    if [ ${TLEAP_STATUS} -ne 0 ] || ! grep -q "Errors = 0" add_ions.leap.out; then
        echo "  *** ERROR: tleap ion addition failed. Check add_ions.leap.out ***"
        cd - > /dev/null
        return 1
    fi

    if [ ! -s "${BSTATE_NAME}_OL3.prmtop" ] || [ ! -s "${BSTATE_NAME}_OL3.rst7" ]; then
        echo "  *** ERROR: Final ionized topology not created. Check add_ions.leap.out ***"
        cd - > /dev/null
        return 1
    fi

    # Get water count from final (ionized) topology
    FINAL_WATERS=$(echo "parminfo" | cpptraj -p ${BSTATE_NAME}_OL3.prmtop 2>/dev/null | awk '$2 == "solvent" && $3 == "molecules." {print $1}')
    echo "  Final water count (post-ions): ${FINAL_WATERS}"

    # Get box info from rst7
    BOX_LINE=$(tail -1 ${BSTATE_NAME}_OL3.rst7)
    echo "  Box: ${BOX_LINE}"

    # ----- Run parmed HRM -----
    echo "  Applying HRM..."
    cat > parmed_hrm.in << PARMED_EOF
addljtype :C,C5,C3,U,U5,U3@H6
changeljpair :C,C5,C3,U,U5,U3@H6 @%OS 2.7 0.050498
changeljpair :C,C5,C3,U,U5,U3@H6 @%O2 2.7 0.050498

addljtype :A,A5,A3,G,G5,G3@H8
changeljpair :A,A5,A3,G,G5,G3@H8 @%OS 2.7 0.050498
changeljpair :A,A5,A3,G,G5,G3@H8 @%O2 2.7 0.050498

addljtype :C,C5,C3,U,U5,U3@H5
changeljpair :C,C5,C3,U,U5,U3@H5 @%OS 2.7 0.050498
changeljpair :C,C5,C3,U,U5,U3@H5 @%O2 2.7 0.050498

changeljpair @%H1 @%OS 2.7 0.051662
changeljpair @%H1 @%O2 2.7 0.051662

outparm ${BSTATE_NAME}_HRM.parm7 ${BSTATE_NAME}_HRM.inpcrd
outparm ${BSTATE_NAME}_HRM.prmtop ${BSTATE_NAME}_HRM.rst7
outpdb ${BSTATE_NAME}_HRM.pdb
PARMED_EOF

    parmed -i parmed_hrm.in \
           -p ${BSTATE_NAME}_OL3.prmtop \
           -c ${BSTATE_NAME}_OL3.rst7 \
           -O > parmed.log 2>&1

    if [ $? -ne 0 ]; then
        echo "  *** ERROR: parmed failed. Check parmed.log ***"
        cd - > /dev/null
        return 1
    fi

    echo "  Done: ${BSTATE_NAME}"
    cd - > /dev/null
    return 0
}

# =============================================================================
# Process all folded bstates
# =============================================================================
if [ "${RUN_FOLDED}" = "true" ]; then
    echo ""
    echo "########################################"
    echo "#       PROCESSING FOLDED BSTATES      #"
    echo "########################################"

    for PDB in ${PDB_DIR}/folded/folded_*.pdb; do
        BASENAME=$(basename ${PDB} .pdb)
        BSTATE_NAME=$(echo ${BASENAME} | sed 's/_frame[0-9]*//')
        ABS_PDB=$(realpath ${PDB})
        process_bstate "${ABS_PDB}" "${BSTATE_NAME}" "${INITIAL_BUFFER_FOLDED}"
    done
else
    echo ""
    echo "  [skipping folded — RUN_FOLDED=false]"
fi

# =============================================================================
# Process all unfolded bstates
# =============================================================================
if [ "${RUN_UNFOLDED}" = "true" ]; then
    echo ""
    echo "########################################"
    echo "#      PROCESSING UNFOLDED BSTATES     #"
    echo "########################################"

    for PDB in ${PDB_DIR}/unfolded/unfolded_*.pdb; do
        BASENAME=$(basename ${PDB} .pdb)
        BSTATE_NAME=$(echo ${BASENAME} | sed 's/_frame[0-9]*//')
        ABS_PDB=$(realpath ${PDB})
        process_bstate "${ABS_PDB}" "${BSTATE_NAME}" "${INITIAL_BUFFER_UNFOLDED}"
    done
else
    echo ""
    echo "  [skipping unfolded — RUN_UNFOLDED=false]"
fi

# =============================================================================
# Process all intermediate bstates
# =============================================================================
if [ "${RUN_INTERMEDIATE}" = "true" ]; then
    echo ""
    echo "########################################"
    echo "#    PROCESSING INTERMEDIATE BSTATES   #"
    echo "########################################"

    for PDB in ${PDB_DIR}/intermediate/intermediate_*.pdb; do
        BASENAME=$(basename ${PDB} .pdb)
        ABS_PDB=$(realpath ${PDB})
        process_bstate "${ABS_PDB}" "${BASENAME}" "${INITIAL_BUFFER_INTERMEDIATE}"
    done
else
    echo ""
    echo "  [skipping intermediate — RUN_INTERMEDIATE=false]"
fi

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "########################################"
echo "#            SUMMARY REPORT            #"
echo "########################################"
echo ""
echo "Target: ${TARGET_WATERS} waters + ${N_NEUTRALIZING} neutralizing K+ + ${N_EXCESS_IONS} K+/${N_EXCESS_IONS} Cl-"
echo ""
printf "  %-20s  %8s  %8s  %4s  %8s  %s\n" "Bstate" "Waters" "Iters" "HRM" "Density" "Box (rst7 last line)"
printf "  %-20s  %8s  %8s  %4s  %8s  %s\n" "--------------------" "--------" "--------" "----" "--------" "------------------------------------"

for DIR in ${OUT_DIR}/folded_* ${OUT_DIR}/unfolded_* ${OUT_DIR}/intermediate_*; do
    [ -d "${DIR}" ] || continue
    NAME=$(basename ${DIR})

    # Water count
    PRMTOP="${DIR}/${NAME}_OL3.prmtop"
    WATERS=""
    if [ -f "${PRMTOP}" ]; then
        WATERS=$(echo "parminfo" | cpptraj -p ${PRMTOP} 2>/dev/null | awk '$2 == "solvent" && $3 == "molecules." {print $1}')
    fi
    [ -z "${WATERS}" ] && WATERS="N/A"

    # Number of Solvate.sh iterations
    ITERS=$(grep -c "^[0-9]*) Buffer:" ${DIR}/solvate.log 2>/dev/null || echo "N/A")

    # HRM check
    if ls ${DIR}/*_HRM.prmtop 1>/dev/null 2>&1; then
        HRM="OK"
    else
        HRM="FAIL"
    fi

    # Box
    RST7="${DIR}/${NAME}_OL3.rst7"
    BOX=""
    if [ -f "${RST7}" ]; then
        BOX=$(tail -1 ${RST7})
    fi

    DENS=$(grep "Density" ${DIR}/leap.log 2>/dev/null | tail -1 | awk '{print $6}')
    [ -z "${DENS}" ] && DENS="N/A"
    
    printf "  %-20s  %8s  %8s  %4s  %8s  %s\n" "${NAME}" "${WATERS}" "${ITERS}" "${HRM}" "${DENS}" "${BOX}"

done

echo ""
echo "VERIFY:"
echo "  1. All water counts (post-ions) should be: $((TARGET_WATERS - N_NEUTRALIZING - 2 * N_EXCESS_IONS))"
echo "     (= ${TARGET_WATERS} pre-ion waters − ${N_NEUTRALIZING} neutralizing K+ − $((2 * N_EXCESS_IONS)) excess K+/Cl⁻)"
echo "  2. Box lines should show three EQUAL lengths and angles"
echo "  3. All HRM = OK"
echo ""
echo "Next: minimize → heat → NPT equilibrate each _HRM system"
