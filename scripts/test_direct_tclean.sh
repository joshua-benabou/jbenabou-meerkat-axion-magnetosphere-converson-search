#!/bin/bash
#SBATCH --job-name=test_tclean
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_direct_tclean_%a.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_direct_tclean_%a.err
#
# Test: can we skip splitting entirely and image directly from the monolithic MS?
# 4 nodes each image a different single channel via tclean.
# This tests both speed and whether parallel tclean jobs interfere.
#

MS="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/SCI-20210212-SS-01/17TBdataset/scratch/kat/1621534878_20231022T19_11_29/1621534878_sdp_l0.ms"
OUTDIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0/direct_tclean"

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

mkdir -p "$OUTDIR"

IDX=$SLURM_ARRAY_TASK_ID
# Pick 4 channels spread across the band
CHANNELS=(17000 17100 17200 17300)
CHAN=${CHANNELS[$IDX]}

echo "=== Direct tclean test: task $IDX on $(hostname) ==="
echo "Channel: $CHAN"
echo "Start: $(date)"

python3 -c "
import time, os, shutil
from casatasks import tclean, exportfits

ms = '${MS}'
chan = ${CHAN}
outdir = '${OUTDIR}'
imagename = os.path.join(outdir, f'direct_chan_{chan}')
fitsname = imagename + '.fits'

# Clean up any prior run
for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
    p = imagename + suffix
    if os.path.isdir(p):
        shutil.rmtree(p)
if os.path.exists(fitsname):
    os.remove(fitsname)

print(f'Imaging channel {chan} directly from monolithic MS...')
t0 = time.time()

tclean(
    vis=ms,
    imagename=imagename,
    field='SgrA*',
    specmode='mfs',
    spw=f'0:{chan}',
    imsize=[512, 512],
    cell='2arcsec',
    niter=0,
    weighting='briggs',
    robust=0.5,
    savemodel='none',
    pbcor=False,
)

t_tclean = time.time() - t0

# Export FITS
exportfits(imagename=imagename + '.image', fitsimage=fitsname, overwrite=True)

t_total = time.time() - t0

# Cleanup CASA products
for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
    p = imagename + suffix
    if os.path.isdir(p):
        shutil.rmtree(p)

size_kb = os.path.getsize(fitsname) / 1024
print(f'DONE. tclean: {t_tclean:.1f}s, total: {t_total:.1f}s, FITS: {size_kb:.0f} KB')
"

echo "End: \$(date)"
