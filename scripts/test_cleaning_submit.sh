#!/bin/bash
#SBATCH --job-name=test_clean
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=240G
#SBATCH --time=04:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_cleaning_%j.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_cleaning_%j.err
#
# Test cleaning strategies on subband 30 (mid-band, low RFI).
# Tests 3 approaches x 6 channels + continuum image.
# Expected runtime: ~1-2 hours.
#
# Submit: sbatch test_cleaning_submit.sh
#

SCRIPT_DIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts"
mkdir -p "${SCRIPT_DIR}/logs"

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

echo "=== Test Cleaning Strategy ==="
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Start: $(date)"

python3 "${SCRIPT_DIR}/test_cleaning_strategy.py" 30

echo "End: $(date)"
echo "Exit code: $?"
