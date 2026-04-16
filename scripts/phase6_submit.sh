#!/bin/bash
#SBATCH --job-name=phase6_sanity
#SBATCH --partition=lr6
#SBATCH --account=pc_heptheory
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=scripts/logs/phase6_sanity_%j.out
#SBATCH --error=scripts/logs/phase6_sanity_%j.err

# Phase 6: Sanity Checks on MeerKAT GC Imaging Data
# READ-ONLY on all FITS files. Produces diagnostic plots in plots/sanity_*.png.

set -euo pipefail

PROJECT_DIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project"

# Create log directory if needed
mkdir -p "${PROJECT_DIR}/scripts/logs"

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

echo "============================================="
echo "Phase 6: Sanity Checks"
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}MB"
echo "============================================="

cd "${PROJECT_DIR}"

python scripts/phase6_sanity_checks.py \
    --workers ${SLURM_CPUS_PER_TASK} \
    --sample-stride 40 \
    --injection-trials 200

echo "============================================="
echo "End time: $(date)"
echo "============================================="
