#!/bin/bash
#SBATCH --job-name=meerkat_image
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --array=0-85
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=90G
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/phase3_subband_%a.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/phase3_subband_%a.err
#
# Phase 3: Channel imaging for all subbands.
# Each subband is independent so all 86 can run on 86 nodes in parallel.
# Within each node, ~30 worker processes image channels concurrently.
# 86 nodes x 30 workers = ~2580 channels being imaged simultaneously.
#
# With 383 channels per subband and ~30 parallel workers per node,
# each subband should finish in ~13x the time of a single channel image.
#

SCRIPT_DIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts"
mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "$(dirname "$SCRIPT_DIR")/images"

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

echo "=== Phase 3: Imaging, subband ${SLURM_ARRAY_TASK_ID} ==="
echo "Node: $(hostname)"
echo "Start: $(date)"

python3 "${SCRIPT_DIR}/phase3_image_channels.py"

echo "End: $(date)"
echo "Exit code: $?"
