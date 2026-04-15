#!/bin/bash
#SBATCH --job-name=debug_clean
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/debug_cleaning_%j.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/debug_cleaning_%j.err
#
# Debug cleaning: test tclean on a single channel with correct noise measurement.
# Tests 7 different threshold strategies on subband 30, channel 100.
# Expected runtime: ~30-60 min (7 tclean runs of ~5-7 min each).
#
# Submit: sbatch debug_cleaning_submit.sh

SCRIPT_DIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts"
mkdir -p "${SCRIPT_DIR}/logs"

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

echo "=== Debug Cleaning ==="
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Start: $(date)"

python3 "${SCRIPT_DIR}/debug_cleaning.py" 30 100

echo "End: $(date)"
echo "Exit code: $?"
