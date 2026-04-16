#!/bin/bash
#SBATCH --job-name=debug_v2
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=240G
#SBATCH --time=04:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/debug_v2_%j.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/debug_v2_%j.err
#
# Debug cleaning v2: test multiscale + masking strategies on single channel.
# Tests 7 configurations: dirty, hogbom, multiscale, with/without masks.
# Expected runtime: ~1-2 hours.
#
# Submit: sbatch debug_cleaning_v2_submit.sh

SCRIPT_DIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts"
mkdir -p "${SCRIPT_DIR}/logs"

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

echo "Start: $(date)"
python3 "${SCRIPT_DIR}/debug_cleaning_v2.py" 30 100
echo "End: $(date)"
echo "Exit code: $?"
