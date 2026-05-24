#!/bin/bash
# =============================================================================
# 02_equilibrate_bstates.sh
#
# Minimization and equilibration pipeline for WE bstates using SGE array jobs.
#
# Instead of submitting 200 individual jobs (10 steps × 20 bstates), this
# script submits 10 array jobs, each with N tasks (one per bstate). The
# scheduler sees 10 job entries. Each array holds on the previous array
# via -hold_jid, so the full pipeline is:
#
#   array_min1 -> array_min2 -> array_md1 -> ... -> array_md3 -> array_md4
#
# Within each array, all bstates run in parallel (subject to queue limits).
#
# WORKFLOW (10 steps per bstate):
#   min1  ->  min2  ->  md1  ->  md2a -> md2b -> md2c -> md2d -> md2e -> md3
#   -> [postmd3 box recalculation] -> md4
#
# INPUT STRUCTURE (from 01_prep_all_bstates.sh):
#   BSTATE_SETUP_DIR/BSTATE_ID/BSTATE_ID_HRM.{prmtop,rst7}
#
# OUTPUT STRUCTURE:
#   02_bstates_equilibration/BSTATE_ID/{min1.out, min1.rst, ..., md4.rst, md4.nc}
#
# USAGE:
#   1. Edit CONFIGURATION section below
#   2. Set BSTATES_TO_RUN to the list of bstate IDs you want to process
#   3. Run:  bash 02_equilibrate_bstates.sh
#
# NOTES:
#   - Minimization steps use pmemd (CPU) since pmemd.cuda doesn't support
#     imin=1. They still request a GPU slot so the array chain is simple,
#     but they finish in seconds.
#   - The postmd3 box recalculation runs at the tail end of the md3 job.
#   - No email notifications on array tasks (per CRC guidelines).
# =============================================================================

set -u

# =============================================================================
#                           CONFIGURATION
# =============================================================================

# ---- AMBER module ----
AMBER_MODULE="amber"

# ---- Input bstate directory (output of 01_prep_all_bstates.sh) ----
# Structure: BSTATE_SETUP_DIR/BSTATE_ID/BSTATE_ID_HRM.{prmtop,rst7}
BSTATE_SETUP_DIR="./01_bstate_setup"

# ---- Output directory ----
OUTPUT_DIR="./02_bstate_equilibration"

# ---- Restraint mask for RNA solute ----
# 14 nucleotides for UUCG tetraloop (PDB: 2KOC)
RESTRAINT_MASK=':1-14'

# ---- Topology file suffix (from prep script) ----
PRMTOP_SUFFIX="_HRM.prmtop"
RST7_SUFFIX="_HRM.rst7"

# ---- SGE queue settings ----
SGE_QUEUE="gpu"
SGE_GPU_CARDS=1
SGE_PE="smp 8"

# ---- Which step to start from ----
# Options: min1, min2, md1, md2a, md2b, md2c, md2d, md2e, md3, md4
# Use this to resume from a specific step if earlier steps already completed.
# The script will assume all prior restart files exist in the output directory.
START_STEP="min1"

# ---- DRY RUN ----
# Set to "true" to generate all files and scripts without submitting to SGE.
# Useful for inspecting the .in files and run scripts before committing.
DRY_RUN="false"

# =============================================================================
#                    BSTATES TO PROCESS
# =============================================================================
# List the bstate IDs you want to run. These must match the directory names
# under BSTATE_SETUP_DIR/
#
# Examples:
#   BSTATES_TO_RUN=("folded_01" "folded_02" "unfolded_01")
#   BSTATES_TO_RUN=("folded_"{01..10} "unfolded_"{01..10})  # all 20
# =============================================================================

BSTATES_TO_RUN=(
    "folded_01"
    "folded_02"
    "folded_03"
    "folded_04"
    "folded_05"
    "folded_06"
    "folded_07"
    "folded_08"
    "folded_09"
    "folded_10"
    "intermediate_01"
    "intermediate_02"
    "intermediate_03"
    "intermediate_04"
    "intermediate_05"
    "intermediate_06"
    "intermediate_07"
    "intermediate_08"
    "intermediate_09"
    "intermediate_10"
    "intermediate_11"
    "intermediate_12"
    "intermediate_13"
    "intermediate_14"
    "intermediate_15"
    "intermediate_16"
    "intermediate_17"
    "intermediate_18"
    "intermediate_19"
    "intermediate_20"
    "intermediate_21"
    "intermediate_22"
    "intermediate_23"
    "intermediate_24"
    "intermediate_25"
    "intermediate_26"
    "intermediate_27"
    "intermediate_28"
    "intermediate_29"
    "intermediate_30"
    "intermediate_31"
    "intermediate_32"
    "intermediate_33"
    "intermediate_34"
    "intermediate_35"
    "intermediate_36"
    "intermediate_37"
    "intermediate_38"
    "intermediate_39"
    "intermediate_40"
    "intermediate_41"
    "intermediate_42"
    "intermediate_43"
    "intermediate_44"
    "intermediate_45"
    "intermediate_46"
    "intermediate_47"
    "intermediate_48"
    "intermediate_49"
    "intermediate_50"
    "unfolded_01"
    "unfolded_02"
    "unfolded_03"
    "unfolded_04"
    "unfolded_05"
    "unfolded_06"
    "unfolded_07"
    "unfolded_08"
    "unfolded_09"
    "unfolded_10"
)

# =============================================================================
#           STEP DEFINITIONS (do not edit unless changing protocol)
# =============================================================================
ALL_STEPS=("min1" "min2" "md1" "md2a" "md2b" "md2c" "md2d" "md2e" "md3" "md4")

# =============================================================================
#                    FUNCTION: Locate bstate source directory
# =============================================================================
find_bstate_source() {
    local BSTATE_ID=$1
    if [ -d "${BSTATE_SETUP_DIR}/${BSTATE_ID}" ]; then
        echo "${BSTATE_SETUP_DIR}/${BSTATE_ID}"
    else
        echo ""
    fi
}

# =============================================================================
#                    FUNCTION: Generate AMBER .in files
# =============================================================================
generate_input_files() {
    local WORKDIR=$1

    # ---- min1: Restrained minimization ----
    cat > ${WORKDIR}/min1.in << EOF
min1 minimization - solvent relaxation, solute restraints
 &cntrl
  imin=1,
  irest=0,
  ntx=1,
  maxcyc=1000,
  ncyc=500,
  ntr=1,
  restraint_wt=500.0,
  restraintmask='${RESTRAINT_MASK}',
  cut=10.0,
  ntpr=100,
  ntwx=0,
 /
EOF

    # ---- min2: Unrestrained minimization ----
    cat > ${WORKDIR}/min2.in << EOF
min2 minimization - full relaxation
 &cntrl
  imin=1,
  irest=0,
  ntx=1,
  maxcyc=2500,
  ncyc=1000,
  cut=10.0,
  ntpr=100,
  ntwx=0,
 /
EOF

    # ---- md1: Heating 0->300K, restrained, NVT, 100ps ----
    cat > ${WORKDIR}/md1.in << EOF
md1 heating 0 to 300K with restraints on RNA 100ps NVT
 &cntrl
  imin=0,
  irest=0, ntx=1,
  ntb=1,
  ntr=1, restraint_wt=25.0, restraintmask='${RESTRAINT_MASK}',
  cut=10.0,
  ntc=2, ntf=2,
  ntt=3, gamma_ln=1.0,
  tempi=0, temp0=300,
  nstlim=50000, dt=0.002,
  ntpr=100,
  ntwx=100,
  ntwr=1000,
 /
EOF

    # ---- md2a-e: NPT restraint release (25->20->15->10->5 kcal), 50ps each ----
    local WEIGHTS=("25.0" "20.0" "15.0" "10.0" "5.0")
    local LABELS=("a" "b" "c" "d" "e")

    for i in 0 1 2 3 4; do
        cat > ${WORKDIR}/md2${LABELS[$i]}.in << EOF
md2${LABELS[$i]} NPT restraint ${WEIGHTS[$i]} kcal 50ps 300K
 &cntrl
  imin=0,
  irest=1, ntx=5,
  ntb=2,
  ntp=1, barostat=2,
  taup=1.0, pres0=1.0,
  ntc=2, ntf=2,
  cut=10.0,
  ntr=1, restraint_wt=${WEIGHTS[$i]}, restraintmask='${RESTRAINT_MASK}',
  ntt=3, gamma_ln=1.0,
  tempi=300.0, temp0=300.0,
  nstlim=25000, dt=0.002,
  ntpr=100, ntwx=1000, ntwr=1000,
 /
EOF
    done

    # ---- md3: Unrestrained NPT equilibration, 200ps ----
    cat > ${WORKDIR}/md3.in << EOF
md3 equilibration 200ps 300K NPT
 &cntrl
  imin=0,
  irest=1, ntx=5,
  ntb=2,
  ntp=1, barostat=2,
  taup=1.0, pres0=1.0,
  cut=10.0,
  ntc=2, ntf=2,
  ntt=3, gamma_ln=1.0,
  tempi=300.0, temp0=300.0,
  ntxo=1,
  nstlim=100000, dt=0.002,
  ntpr=100, ntwx=1000, ntwr=1000,
 /
EOF

    # ---- md4: NVT equilibration, 1ns ----
    cat > ${WORKDIR}/md4.in << EOF
md4 equilibration 1ns 300K NVT
 &cntrl
  imin=0,
  irest=1, ntx=5,
  ntb=1,
  cut=10.0,
  ntc=2, ntf=2,
  ntt=3, gamma_ln=1.0,
  tempi=300.0, temp0=300.0,
  nstlim=500000, dt=0.002,
  ntpr=500, ntwx=1000, ntwr=100000,
 /
EOF
}

# =============================================================================
#    FUNCTION: Generate the postmd3 box recalculation Python script
# =============================================================================
generate_postmd3_script() {
    local WORKDIR=$1

    cat > ${WORKDIR}/postmd3_calcboxlength.py << 'PYEOF'
#!/usr/bin/env python3
"""
postmd3_calcboxlength.py — Automated box volume correction for truncated octahedron.

Reads the box length from the last line of md3.rst, extracts the average VOLUME
from md3.out, computes the corrected box length, and writes md3_NewVolume.rst
with the updated box dimensions.

For a truncated octahedron: V = (1/2) * (4/3)^(3/2) * L^3
"""
import sys
import re
import shutil

def get_box_length_from_rst(rst_file):
    """Read the three box lengths from the last line of an AMBER rst file."""
    with open(rst_file, 'r') as f:
        lines = f.readlines()
    last_line = lines[-1].strip()
    vals = last_line.split()
    if len(vals) < 6:
        print(f"ERROR: Expected 6 values on last line of {rst_file}, got {len(vals)}")
        sys.exit(1)
    return float(vals[0]), float(vals[1]), float(vals[2])

def get_avg_volume_from_mdout(mdout_file):
    """Extract the average volume from an AMBER mdout file.
    Looks for 'VOLUME' in the A V E R A G E S section."""
    with open(mdout_file, 'r') as f:
        content = f.read()

    avg_match = re.search(r'A V E R A G E S.*?R M S', content, re.DOTALL)
    if not avg_match:
        print("WARNING: Could not find AVERAGES section, averaging all VOLUME lines")
        volumes = re.findall(r'VOLUME\s*=\s*([\d.]+)', content)
        if not volumes:
            print(f"ERROR: No VOLUME entries found in {mdout_file}")
            sys.exit(1)
        vols = [float(v) for v in volumes]
        return sum(vols) / len(vols)

    avg_section = avg_match.group(0)
    vol_match = re.search(r'VOLUME\s*=\s*([\d.]+)', avg_section)
    if not vol_match:
        print(f"ERROR: No VOLUME in AVERAGES section of {mdout_file}")
        sys.exit(1)
    return float(vol_match.group(1))

def main():
    rst_file = "md3.rst"
    mdout_file = "md3.out"
    output_file = "md3_NewVolume.rst"

    L1, L2, L3 = get_box_length_from_rst(rst_file)
    print(f"Current box lengths: {L1:.7f}  {L2:.7f}  {L3:.7f}")
    length = L1

    predicted_volume = 0.5 * (4.0/3.0)**1.5 * length**3
    print(f"Predicted volume from rst: {predicted_volume:.4f}")

    avg_volume = get_avg_volume_from_mdout(mdout_file)
    print(f"Average volume from md3.out: {avg_volume:.4f}")

    scale = (avg_volume / predicted_volume) ** (1.0/3.0)
    new_length = scale * length
    print(f"Scale factor: {scale:.10f}")
    print(f"New box length: {new_length:.7f}")

    shutil.copy2(rst_file, output_file)
    with open(output_file, 'r') as f:
        lines = f.readlines()

    old_vals = lines[-1].split()
    new_last_line = f"{new_length:12.7f}{new_length:12.7f}{new_length:12.7f}"
    new_last_line += f"{float(old_vals[3]):12.7f}{float(old_vals[4]):12.7f}{float(old_vals[5]):12.7f}\n"
    lines[-1] = new_last_line

    with open(output_file, 'w') as f:
        f.writelines(lines)

    print(f"Wrote {output_file} with corrected box dimensions.")

if __name__ == "__main__":
    main()
PYEOF

    chmod +x ${WORKDIR}/postmd3_calcboxlength.py
}

# =============================================================================
#    FUNCTION: Get step index
# =============================================================================
get_step_index() {
    local STEP=$1
    for i in "${!ALL_STEPS[@]}"; do
        if [ "${ALL_STEPS[$i]}" == "${STEP}" ]; then
            echo $i
            return
        fi
    done
    echo -1
}

# =============================================================================
#    Helper: Return the input restart expression for a step.
#    Uses ${PRMTOP_BASE} which resolves at runtime inside the SGE script.
# =============================================================================
get_input_restart_expr() {
    local STEP=$1
    case ${STEP} in
        min1) echo '${PRMTOP_BASE}.rst7' ;;
        min2) echo "min1.rst" ;;
        md1)  echo "min2.rst" ;;
        md2a) echo "md1.rst" ;;
        md2b) echo "md2a.rst" ;;
        md2c) echo "md2b.rst" ;;
        md2d) echo "md2c.rst" ;;
        md2e) echo "md2d.rst" ;;
        md3)  echo "md2e.rst" ;;
        md4)  echo "md3_NewVolume.rst" ;;
    esac
}

get_ref_restart_expr() {
    local STEP=$1
    case ${STEP} in
        min1) echo '${PRMTOP_BASE}.rst7' ;;
        min2) echo "min1.rst" ;;
        md1)  echo "min2.rst" ;;
        md2a) echo "md1.rst" ;;
        md2b) echo "md2a.rst" ;;
        md2c) echo "md2b.rst" ;;
        md2d) echo "md2c.rst" ;;
        md2e) echo "md2d.rst" ;;
        md3)  echo "md2e.rst" ;;
        md4)  echo "md3.rst" ;;
    esac
}

# =============================================================================
#    FUNCTION: Generate one SGE array script for a given step
#
#    Each array task reads bstate_map.txt to find its bstate ID and working
#    directory based on $SGE_TASK_ID, then runs the appropriate AMBER command.
# =============================================================================
generate_array_script() {
    local STEP=$1
    local ABS_OUTPUT_DIR=$2
    local N_TASKS=$3
    local SCRIPT="${ABS_OUTPUT_DIR}/run_${STEP}.sh"

    # ---- Get restart file expressions for this step ----
    local INPUT_RST_EXPR
    local REF_RST_EXPR
    INPUT_RST_EXPR=$(get_input_restart_expr "${STEP}")
    REF_RST_EXPR=$(get_ref_restart_expr "${STEP}")

    # ---- Write the script header (SGE directives) ----
    cat > ${SCRIPT} << SGEEOF
#!/bin/bash
#$ -pe ${SGE_PE}
#$ -q ${SGE_QUEUE}
#$ -l gpu_card=${SGE_GPU_CARDS}
#$ -N eq_${STEP}
#$ -t 1-${N_TASKS}
#$ -j y
#$ -o ${ABS_OUTPUT_DIR}/logs/${STEP}_\$TASK_ID.log

# ---- Load AMBER ----
module load ${AMBER_MODULE}

# ---- Read bstate mapping: SGE_TASK_ID -> bstate working directory ----
MAPFILE="${ABS_OUTPUT_DIR}/bstate_map.txt"
BSTATE_ID=\$(awk -v id=\$SGE_TASK_ID 'NR==id {print \$1}' \${MAPFILE})
WORKDIR=\$(awk -v id=\$SGE_TASK_ID 'NR==id {print \$2}' \${MAPFILE})

if [ -z "\${WORKDIR}" ] || [ ! -d "\${WORKDIR}" ]; then
    echo "ERROR: Task \${SGE_TASK_ID} — cannot find workdir from \${MAPFILE}"
    exit 1
fi

cd \${WORKDIR}

# ---- Identify topology ----
PRMTOP=\$(ls *_HRM.prmtop 2>/dev/null | head -1)
if [ -z "\${PRMTOP}" ]; then
    echo "ERROR: No *_HRM.prmtop found in \${WORKDIR}"
    exit 1
fi
PRMTOP_BASE="\${PRMTOP%.prmtop}"

# ---- Input and reference restart files ----
INPUT_RST="${INPUT_RST_EXPR}"
REF_RST="${REF_RST_EXPR}"

echo "[\$(date)] Starting ${STEP} for \${BSTATE_ID} (task \${SGE_TASK_ID})"
echo "  Topology: \${PRMTOP}"
echo "  Input:    \${INPUT_RST}"
echo "  Ref:      \${REF_RST}"
SGEEOF

    # ---- Step-specific AMBER command ----
    if [[ "${STEP}" == min* ]]; then
        # Minimization: CPU pmemd (pmemd.cuda does not support imin=1)
        cat >> ${SCRIPT} << SGEEOF

\$AMBERHOME/bin/pmemd -O \\
    -i ${STEP}.in \\
    -o ${STEP}.out \\
    -p \${PRMTOP} \\
    -c \${INPUT_RST} \\
    -r ${STEP}.rst \\
    -ref \${REF_RST}
STATUS=\$?
echo "[\$(date)] ${STEP} finished with exit code \${STATUS}"
exit \${STATUS}
SGEEOF

    elif [[ "${STEP}" == "md3" ]]; then
        # md3 + postmd3 box recalculation
        cat >> ${SCRIPT} << SGEEOF

\$AMBERHOME/bin/pmemd.cuda -O \\
    -i ${STEP}.in \\
    -o ${STEP}.out \\
    -p \${PRMTOP} \\
    -c \${INPUT_RST} \\
    -r ${STEP}.rst \\
    -x ${STEP}.nc \\
    -ref \${REF_RST}
STATUS=\$?

if [ \${STATUS} -ne 0 ]; then
    echo "[\$(date)] ERROR: ${STEP} failed with exit code \${STATUS}"
    exit \${STATUS}
fi

echo "[\$(date)] ${STEP} completed. Running box recalculation..."
python3 postmd3_calcboxlength.py
BOX_STATUS=\$?

if [ \${BOX_STATUS} -ne 0 ]; then
    echo "[\$(date)] ERROR: postmd3 box recalculation failed"
    exit \${BOX_STATUS}
fi

echo "[\$(date)] Box recalculation complete. md3_NewVolume.rst written."
exit 0
SGEEOF

    else
        # All other MD steps: pmemd.cuda
        cat >> ${SCRIPT} << SGEEOF

\$AMBERHOME/bin/pmemd.cuda -O \\
    -i ${STEP}.in \\
    -o ${STEP}.out \\
    -p \${PRMTOP} \\
    -c \${INPUT_RST} \\
    -r ${STEP}.rst \\
    -x ${STEP}.nc \\
    -ref \${REF_RST}
STATUS=\$?
echo "[\$(date)] ${STEP} finished with exit code \${STATUS}"
exit \${STATUS}
SGEEOF
    fi

    chmod +x ${SCRIPT}
}

# =============================================================================
#                              MAIN
# =============================================================================

# ---- Validate START_STEP ----
START_IDX=$(get_step_index "${START_STEP}")
if [ "${START_IDX}" -eq -1 ]; then
    log "ERROR: Invalid START_STEP '${START_STEP}'"
    log "Valid options: ${ALL_STEPS[*]}"
    exit 1
fi

# ---- Create output directory structure ----
ABS_OUTPUT_DIR="$(mkdir -p "${OUTPUT_DIR}" && realpath "${OUTPUT_DIR}")"
mkdir -p "${ABS_OUTPUT_DIR}/logs"

# ---- Log file ----
SCRIPT_LOG="${ABS_OUTPUT_DIR}/02_equilibrate_bstates.log"
> "${SCRIPT_LOG}"
log() {
    echo "$@" | tee -a "${SCRIPT_LOG}"
}

log ""
log "============================================================"
log "  AMBER Min/Eq Pipeline for WE Bstates (Array Jobs)"
log "  $(date)"
log "============================================================"
log ""
log "  Setup dir:      ${BSTATE_SETUP_DIR}"
log "  Output dir:     ${OUTPUT_DIR}"
log "  Restraint mask: ${RESTRAINT_MASK}"
log "  Start step:     ${START_STEP}"
log "  Dry run:        ${DRY_RUN}"
log "  Bstates:        ${#BSTATES_TO_RUN[@]}"
log "  Log:            ${SCRIPT_LOG}"
log ""

# =============================================================================
# PHASE 1: Prepare all bstate directories (runs locally, no job submission)
# =============================================================================
log "============================================================"
log "  Phase 1: Preparing bstate directories"
log "============================================================"
log ""

N_BSTATES=${#BSTATES_TO_RUN[@]}
PREP_SUCCESS=0
PREP_FAIL=0

# Write bstate mapping file: one line per task, "BSTATE_ID  /abs/path/to/workdir"
MAPFILE="${ABS_OUTPUT_DIR}/bstate_map.txt"
> ${MAPFILE}

for idx in "${!BSTATES_TO_RUN[@]}"; do
    BSTATE_ID="${BSTATES_TO_RUN[$idx]}"
    TASK_NUM=$((idx + 1))

    log "  [${TASK_NUM}/${N_BSTATES}] ${BSTATE_ID}..."

    # ---- Locate source ----
    SRC_DIR=$(find_bstate_source "${BSTATE_ID}")
    if [ -z "${SRC_DIR}" ]; then
        log "    *** ERROR: Cannot find ${BSTATE_ID} under ${BSTATE_SETUP_DIR}/"
        PREP_FAIL=$((PREP_FAIL + 1))
        continue
    fi

    # ---- Verify source files ----
    PRMTOP_FILE="${SRC_DIR}/${BSTATE_ID}${PRMTOP_SUFFIX}"
    RST7_FILE="${SRC_DIR}/${BSTATE_ID}${RST7_SUFFIX}"

    if [ ! -f "${PRMTOP_FILE}" ]; then
        log "    *** ERROR: Topology not found: ${PRMTOP_FILE}"
        PREP_FAIL=$((PREP_FAIL + 1))
        continue
    fi
    if [ ! -f "${RST7_FILE}" ]; then
        log "    *** ERROR: Coordinates not found: ${RST7_FILE}"
        PREP_FAIL=$((PREP_FAIL + 1))
        continue
    fi

    # ---- Create bstate workdir ----
    WORKDIR="${ABS_OUTPUT_DIR}/${BSTATE_ID}"
    mkdir -p "${WORKDIR}"

    # ---- Copy topology and coordinates (each bstate has a unique topology) ----
    PRMTOP_BASE="${BSTATE_ID}_HRM"
    cp "${PRMTOP_FILE}" "${WORKDIR}/${PRMTOP_BASE}.prmtop"
    cp "${RST7_FILE}" "${WORKDIR}/${PRMTOP_BASE}.rst7"

    # ---- Generate AMBER .in files ----
    generate_input_files "${WORKDIR}"

    # ---- Generate postmd3 Python script ----
    generate_postmd3_script "${WORKDIR}"

    # ---- Write to mapping file ----
    echo "${BSTATE_ID}  ${WORKDIR}" >> ${MAPFILE}

    PREP_SUCCESS=$((PREP_SUCCESS + 1))
done

log ""
log "  Prepared: ${PREP_SUCCESS} / ${N_BSTATES}"
if [ ${PREP_FAIL} -gt 0 ]; then
    log "  Failed:   ${PREP_FAIL} / ${N_BSTATES}"
    log ""
    log "  *** Fix errors above before continuing. ***"
    log "  *** Mapping file may be incomplete — do not submit. ***"
    exit 1
fi

# Verify mapping file has the right number of lines
MAP_LINES=$(wc -l < ${MAPFILE})
if [ "${MAP_LINES}" -ne "${N_BSTATES}" ]; then
    log "  *** ERROR: Mapping file has ${MAP_LINES} lines, expected ${N_BSTATES}"
    exit 1
fi

log ""
log "  Mapping file: ${MAPFILE}"
log "  Logs:         ${ABS_OUTPUT_DIR}/logs/"
log ""

# =============================================================================
# PHASE 2: Generate and submit array jobs
# =============================================================================
log "============================================================"
log "  Phase 2: Submitting array jobs"
log "============================================================"
log ""

PREV_JOBID=""
FIRST_JOBID=""
LAST_JOBID=""

for i in "${!ALL_STEPS[@]}"; do
    STEP="${ALL_STEPS[$i]}"

    # Skip steps before START_STEP
    if [ "$i" -lt "${START_IDX}" ]; then
        continue
    fi

    # ---- Generate the array script ----
    generate_array_script "${STEP}" "${ABS_OUTPUT_DIR}" "${N_BSTATES}"

    if [ "${DRY_RUN}" == "true" ]; then
        log "  ${STEP}: array script generated (dry run, not submitted)"
        PREV_JOBID="DRY"
        continue
    fi

    # ---- Build qsub command with dependency on previous array ----
    QSUB_CMD="qsub"
    if [ -n "${PREV_JOBID}" ]; then
        QSUB_CMD="qsub -hold_jid ${PREV_JOBID}"
    fi

    # ---- Submit ----
    SUBMIT_OUT=$(${QSUB_CMD} "${ABS_OUTPUT_DIR}/run_${STEP}.sh" 2>&1)
    SUBMIT_STATUS=$?

    if [ ${SUBMIT_STATUS} -ne 0 ]; then
        log "  *** ERROR submitting ${STEP}: ${SUBMIT_OUT}"
        log "  *** Stopping. Previously submitted arrays will still run."
        exit 1
    fi

    # Extract job ID from: "Your job-array 123456.1-20:1 ("eq_min1") has been submitted"
    PREV_JOBID=$(echo "${SUBMIT_OUT}" | awk '{print $3}' | cut -d. -f1)
    if [ -z "${FIRST_JOBID}" ]; then
        FIRST_JOBID="${PREV_JOBID}"
    fi
    LAST_JOBID="${PREV_JOBID}"
    log "  ${STEP} -> array job ${PREV_JOBID} (${N_BSTATES} tasks)"

done

# =============================================================================
# PHASE 3: Generate final structure extraction job
#
# After md4 completes for all bstates, extract the last frame from each
# md4 trajectory into 03_bstates_final/struct_X/ as struct.ncrst, struct.rst,
# and struct.pdb. Then generate bstates.txt for WESTPA.
# =============================================================================
log "============================================================"
log "  Phase 3: Final structure extraction"
log "============================================================"
log ""

FINAL_DIR="$(mkdir -p "./03_bstates_final" && realpath "./03_bstates_final")"
mkdir -p "${FINAL_DIR}/StructureFiles"

# ---- Generate the extraction array script ----
cat > ${ABS_OUTPUT_DIR}/run_extract_final.sh << SGEEOF
#!/bin/bash
#$ -pe ${SGE_PE}
#$ -q ${SGE_QUEUE}
#$ -l gpu_card=${SGE_GPU_CARDS}
#$ -N eq_extract
#$ -t 1-${N_BSTATES}
#$ -j y
#$ -o ${ABS_OUTPUT_DIR}/logs/extract_\$TASK_ID.log

module load ${AMBER_MODULE}

# ---- Read bstate mapping ----
MAPFILE="${ABS_OUTPUT_DIR}/bstate_map.txt"
BSTATE_ID=\$(awk -v id=\$SGE_TASK_ID 'NR==id {print \$1}' \${MAPFILE})
WORKDIR=\$(awk -v id=\$SGE_TASK_ID 'NR==id {print \$2}' \${MAPFILE})

if [ -z "\${WORKDIR}" ] || [ ! -d "\${WORKDIR}" ]; then
    echo "ERROR: Task \${SGE_TASK_ID} — cannot find workdir"
    exit 1
fi

# ---- Identify topology ----
PRMTOP=\$(ls \${WORKDIR}/*_HRM.prmtop 2>/dev/null | head -1)
if [ -z "\${PRMTOP}" ]; then
    echo "ERROR: No topology found in \${WORKDIR}"
    exit 1
fi

# ---- Verify md4 completed ----
if [ ! -f "\${WORKDIR}/md4.rst" ]; then
    echo "ERROR: md4.rst not found in \${WORKDIR} — md4 may not have completed"
    exit 1
fi

# ---- Create output directory: struct_0, struct_1, ... (0-indexed) ----
STRUCT_IDX=\$((\${SGE_TASK_ID} - 1))
STRUCT_DIR="${FINAL_DIR}/StructureFiles/struct_\${STRUCT_IDX}"
mkdir -p \${STRUCT_DIR}

echo "[\$(date)] Extracting final structure for \${BSTATE_ID} -> struct_\${STRUCT_IDX}"

cd \${STRUCT_DIR}

# ---- Extract last frame from md4 as .ncrst, .rst, .rst7, and .pdb ----
ccat > extract.cpptraj << CPPTRAJ_EOF
parm \${PRMTOP}
trajin \${WORKDIR}/md4.nc lastframe
trajout struct.ncrst ncrestart
trajout struct.rst restart
trajout struct.rst7 restart
go
clear trajin
trajin \${WORKDIR}/md4.nc lastframe
autoimage
trajout struct.pdb pdb
go
quit
CPPTRAJ_EOF

cpptraj -i extract.cpptraj > extract.log 2>&1
STATUS=\$?

if [ \${STATUS} -ne 0 ]; then
    log "ERROR: cpptraj extraction failed. Check \${STRUCT_DIR}/extract.log"
    exit \${STATUS}
fi

# ---- Also copy the topology for reference ----
cp \${PRMTOP} \${STRUCT_DIR}/struct.prmtop

log "[\$(date)] Done: \${STRUCT_DIR}/"
ls -la \${STRUCT_DIR}/struct.*
exit 0
SGEEOF

chmod +x ${ABS_OUTPUT_DIR}/run_extract_final.sh

# ---- Generate bstates.txt (can be written now since we know N and ordering) ----
BSTATES_TXT="${FINAL_DIR}/bstates.txt"
> ${BSTATES_TXT}

for idx in $(seq 0 $((N_BSTATES - 1))); do
    WEIGHT=$(python3 -c "print(f'{1.0/${N_BSTATES}:.17e}')")
    echo "${idx}    ${WEIGHT}    StructureFiles/struct_${idx}" >> ${BSTATES_TXT}
done

# ---- Generate mapping log: struct index -> original bstate ID and source ----
MAPPING_LOG="${FINAL_DIR}/struct_mapping.log"
cat > ${MAPPING_LOG} << MAPEOF
# =============================================================================
# Bstate -> Structure Mapping
# Generated: $(date)
#
# This file records the correspondence between the original bstate IDs
# (from 01_prep_all_bstates.sh) and the final struct_X directories used
# by WESTPA. Each structure was extracted as the last frame of md4.
#
# Columns: struct_index  bstate_id  source_directory  equilibration_directory
# =============================================================================
MAPEOF

for idx in "${!BSTATES_TO_RUN[@]}"; do
    BSTATE_ID="${BSTATES_TO_RUN[$idx]}"
    SRC_DIR=$(find_bstate_source "${BSTATE_ID}")
    EQ_DIR="${ABS_OUTPUT_DIR}/${BSTATE_ID}"
    echo "struct_${idx}    ${BSTATE_ID}    ${SRC_DIR}    ${EQ_DIR}" >> ${MAPPING_LOG}
done

log "  Final structure dir: ${FINAL_DIR}"
log "  bstates.txt:         ${BSTATES_TXT}"
log "  Mapping log:         ${MAPPING_LOG}"
log ""
log "  --- bstates.txt ---"
cat ${BSTATES_TXT} | tee -a "${SCRIPT_LOG}"
log ""
log "  --- struct_mapping.log ---"
cat ${MAPPING_LOG} | tee -a "${SCRIPT_LOG}"
log ""

# ---- Submit extraction job (holds on md4 array) ----
if [ "${DRY_RUN}" == "true" ]; then
    log "  extract: array script generated (dry run, not submitted)"
else
    QSUB_CMD="qsub"
    if [ -n "${PREV_JOBID}" ]; then
        QSUB_CMD="qsub -hold_jid ${PREV_JOBID}"
    fi

    SUBMIT_OUT=$(${QSUB_CMD} "${ABS_OUTPUT_DIR}/run_extract_final.sh" 2>&1)
    SUBMIT_STATUS=$?

    if [ ${SUBMIT_STATUS} -ne 0 ]; then
        log "  *** ERROR submitting extract: ${SUBMIT_OUT}"
        exit 1
    fi

    EXTRACT_JOBID=$(echo "${SUBMIT_OUT}" | awk '{print $3}' | cut -d. -f1)
    LAST_JOBID="${EXTRACT_JOBID}"
    log "  extract -> array job ${EXTRACT_JOBID} (${N_BSTATES} tasks)"
fi

log ""

# =============================================================================
#                           SUMMARY
# =============================================================================
log ""
log "########################################"
log "#            SUMMARY                   #"
log "########################################"
log ""
log "  Bstates:    ${N_BSTATES}"
log "  Steps:      ${ALL_STEPS[*]:${START_IDX}} + extract"
log "  Array jobs: $((${#ALL_STEPS[@]} - START_IDX + 1)) (one per step, ${N_BSTATES} tasks each)"
log ""
log "  Equilibration: ${ABS_OUTPUT_DIR}/<BSTATE_ID>/"
log "  Final output:  ${FINAL_DIR}/struct_<X>/"
log "  Logs:          ${ABS_OUTPUT_DIR}/logs/"
log "  bstates.txt:   ${BSTATES_TXT}"
log ""
log "  Bstate mapping (SGE_TASK_ID -> bstate -> struct index):"
awk '{printf "    Task %d -> %-20s -> struct_%d\n", NR, $1, NR-1}' ${MAPFILE} | tee -a "${SCRIPT_LOG}"
log ""
if [ "${DRY_RUN}" == "true" ]; then
    log "  *** DRY RUN — no jobs were submitted ***"
    log "  Inspect the generated files, then set DRY_RUN=\"false\" and re-run."
    log ""
else
    log "  Monitor:    qstat -u $(whoami)"
    log "  Kill all:   qdel ${FIRST_JOBID} ${LAST_JOBID}"
    log "              (or: qdel \$(qstat -u \$(whoami) | grep eq_ | awk '{print \$1}' | sort -u))"
    log ""
fi
log "  After completion, verify:"
log "    1. Each bstate has md4.rst:  ls ${ABS_OUTPUT_DIR}/*/md4.rst"
log "    2. Final structures:         ls ${FINAL_DIR}/StructureFiles/struct_*/struct.pdb"
log "    3. No errors in logs:        grep -l ERROR ${ABS_OUTPUT_DIR}/logs/*.log"
log "    4. Min convergence:          grep 'FINAL RESULTS' ${ABS_OUTPUT_DIR}/*/min*.out"
log ""
