#!/bin/bash
#SBATCH --job-name=meerkat_img
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --array=0-85
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=240G
#SBATCH --time=10:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/phase3_cube_%a.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/phase3_cube_%a.err
#
# Phase 3: Image all channels in each subband using tclean cube mode (dirty, niter=0).
# One tclean call per subband, all 86 subbands in parallel on 86 nodes.
# Each produces 383 individual FITS files (512x512, 2arcsec, Briggs r=0.5).
# Jobs that find a missing subband MS will skip gracefully.
#

SCRIPT_DIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts"
mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "$(dirname "$SCRIPT_DIR")/images"

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

echo "=== Phase 3: Imaging subband ${SLURM_ARRAY_TASK_ID} ==="
echo "Node: $(hostname)"
echo "Start: $(date)"

python3 "${SCRIPT_DIR}/phase3_image_cube.py"

echo "End: $(date)"
echo "Exit code: $?"
