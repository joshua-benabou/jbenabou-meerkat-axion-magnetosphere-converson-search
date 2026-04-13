#!/bin/bash
#SBATCH --job-name=axion_search
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=240G
#SBATCH --time=04:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/phase7_axion_search.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/phase7_axion_search.err
#
# Phase 7: Axion search via sideband background subtraction.
#
# Runs the search pipeline across all NS targets in the template bank.
# Each target loads ~100 channel FITS images (sideband), so memory usage
# scales as ~100 * 512*512 * 4 bytes ~ 100 MB per target (very modest).
#
# The main bottleneck is I/O: reading ~100 FITS files per target.
# With ~5 targets (dummy bank), this should complete in minutes.
# With a full template bank of hundreds of NS, expect ~1-2 hours.
#

SCRIPT_DIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts"
mkdir -p "${SCRIPT_DIR}/logs"

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

echo "=== Phase 7: Axion Search ==="
echo "Node: $(hostname)"
echo "Start: $(date)"

python3 "${SCRIPT_DIR}/phase7_axion_search.py" \
    --n_sideband 50 \
    --n_guard 5 \
    --threshold 5.0 \
    --bg_method median

echo "End: $(date)"
echo "Exit code: $?"
