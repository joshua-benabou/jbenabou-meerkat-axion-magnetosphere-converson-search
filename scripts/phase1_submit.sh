#!/bin/bash
#SBATCH --job-name=meerkat_split
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --array=0-85%8
#SBATCH --requeue
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/phase1_subband_%a.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/phase1_subband_%a.err
#
# Phase 1: Split monolithic MS into 86 subbands with TOPO->LSRK conversion.
#
# Array jobs 0-85, with %8 = max 8 running simultaneously.
# Input MS is striped across 8 Lustre HDD OSTs, so ~8 parallel readers
# is roughly one per OST. Increase if I/O isn't saturated; decrease if
# jobs are too slow from contention.
#
# Total channels: 32768
# Channels per subband: 383 (~10 MHz at 26.123 kHz)
# Number of subbands: 86 (last subband has 32768 - 85*383 = 213 channels)
#
# Estimated time per subband: ~30-60 min (TBD from test results)
# Estimated total wall time at %4: ~10-15 hours
#

SCRIPT_DIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Create log directory
mkdir -p "${PROJECT_DIR}/scripts/logs"

# Activate CASA environment
source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

echo "=== Phase 1: Subband ${SLURM_ARRAY_TASK_ID} ==="
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

python3 "${SCRIPT_DIR}/phase1_split_subbands.py"

echo ""
echo "End: $(date)"
echo "Exit code: $?"
