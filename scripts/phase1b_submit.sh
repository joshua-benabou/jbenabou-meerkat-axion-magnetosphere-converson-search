#!/bin/bash
#SBATCH --job-name=meerkat_rfi
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --array=0-85
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/phase1b_subband_%a.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/phase1b_subband_%a.err
#
# Phase 1b: RFI flagging on all subbands.
# Each subband is independent — all 86 run in parallel.
# Applies known MeerKAT L-band RFI masks + tfcrop automatic flagging.
# Run AFTER Phase 1 (split) and BEFORE Phase 2 (contsub).
#

SCRIPT_DIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts"
mkdir -p "${SCRIPT_DIR}/logs"

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

echo "=== Phase 1b: RFI flagging, subband ${SLURM_ARRAY_TASK_ID} ==="
echo "Node: $(hostname)"
echo "Start: $(date)"

python3 "${SCRIPT_DIR}/phase1b_flag_rfi.py"

echo "End: $(date)"
echo "Exit code: $?"
