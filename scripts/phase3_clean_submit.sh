#!/bin/bash
#SBATCH --job-name=meerkat_clean
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --array=0-85
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=240G
#SBATCH --time=10:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/phase3_clean_%a.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/phase3_clean_%a.err
#
# Phase 3b: Cleaned imaging with 10% peak threshold (no mask).
# Measures peak from existing dirty images, sets threshold = 10% of peak.
# One tclean call per subband, all 86 subbands in parallel.
# Output: images/subband_XXX/cleaned/chan_NNNN_FFFF.FFFMHz.fits
#
# Test single subband first:
#   sbatch --array=30 phase3_clean_submit.sh
#
# Then full run:
#   sbatch phase3_clean_submit.sh
#

SCRIPT_DIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts"
mkdir -p "${SCRIPT_DIR}/logs"

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

echo "=== Phase 3b: Cleaned imaging subband ${SLURM_ARRAY_TASK_ID} ==="
echo "Node: $(hostname)"
echo "Start: $(date)"

python3 "${SCRIPT_DIR}/phase3_clean_cube.py"

echo "End: $(date)"
echo "Exit code: $?"
