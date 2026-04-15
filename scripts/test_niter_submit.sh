#!/bin/bash
#SBATCH --job-name=test_niter
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=240G
#SBATCH --time=04:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_niter_%j.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_niter_%j.err

SCRIPT_DIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts"
mkdir -p "${SCRIPT_DIR}/logs"

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

echo "Start: $(date)"
python3 "${SCRIPT_DIR}/test_niter.py"
echo "End: $(date)"
