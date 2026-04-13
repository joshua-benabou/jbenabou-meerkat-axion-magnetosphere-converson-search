#!/bin/bash
#SBATCH --job-name=phase5_rfi_flag
#SBATCH --partition=lr6
#SBATCH --account=pc_heptheory
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=scripts/logs/phase5_rfi_%j.out
#SBATCH --error=scripts/logs/phase5_rfi_%j.err

# Phase 5: RFI Channel Flagging
# Scans all ~32,768 dirty FITS images and produces a channel quality CSV + plot.
# READ-ONLY on all FITS files.

set -euo pipefail

PROJECT_DIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project"

# Create log directory if needed
mkdir -p "${PROJECT_DIR}/scripts/logs"

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

echo "============================================="
echo "Phase 5: RFI Channel Flagging"
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}MB"
echo "============================================="

cd "${PROJECT_DIR}"

python scripts/phase5_rfi_flagging.py \
    --workers ${SLURM_CPUS_PER_TASK} \
    --image-dir "${PROJECT_DIR}/images" \
    --output-csv "${PROJECT_DIR}/rfi_channel_flags.csv" \
    --output-plot "${PROJECT_DIR}/plots/rfi_overview.png"

echo "============================================="
echo "End time: $(date)"
echo "============================================="
