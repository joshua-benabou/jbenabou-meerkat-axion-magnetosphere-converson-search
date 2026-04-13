#!/bin/bash
#SBATCH --job-name=test_16n_split
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --array=0-15
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_16node_split_%a.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_16node_split_%a.err
#
# 16-node parallel split test to find I/O saturation point.
# Each task splits a different 383-channel subband.
#

TEST_DIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0"
MS="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/SCI-20210212-SS-01/17TBdataset/scratch/kat/1621534878_20231022T19_11_29/1621534878_sdp_l0.ms"

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

IDX=$SLURM_ARRAY_TASK_ID
CHAN_START=$((IDX * 383))
CHAN_END=$((CHAN_START + 382))
OUT="${TEST_DIR}/n16_${IDX}.ms"

echo "=== 16-node split test: task $IDX on $(hostname) ==="
echo "Channels: ${CHAN_START}~${CHAN_END}"
echo "Start: $(date)"

python3 -c "
import time, os, shutil
from casatasks import split
import casatools

out = '${OUT}'
if os.path.exists(out):
    shutil.rmtree(out)

t0 = time.time()
split(vis='${MS}', outputvis=out, field='SgrA*', spw='0:${CHAN_START}~${CHAN_END}', datacolumn='data')
elapsed = time.time() - t0

msmd = casatools.msmetadata()
msmd.open(out)
nc = msmd.nchan(0)
nr = msmd.nrows()
msmd.close()

size_mb = sum(os.path.getsize(os.path.join(dp, f)) for dp, dn, fn in os.walk(out) for f in fn) / 1e6
print(f'DONE. Time: {elapsed:.1f}s ({elapsed/60:.1f} min), Size: {size_mb:.1f} MB, Chans: {nc}, Rows: {nr:.0f}')
"

echo "End: $(date)"
