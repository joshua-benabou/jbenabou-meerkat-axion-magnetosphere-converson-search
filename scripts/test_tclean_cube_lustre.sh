#!/bin/bash
#SBATCH --job-name=test_cube2
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=90G
#SBATCH --time=10:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_tclean_cube_lustre.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_tclean_cube_lustre.err
#
# Test tclean cube mode directly from Lustre (no /tmp copy).
# Uses a completed production subband.
#

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

# Use first completed production subband
SUBBAND="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/subbands/subband_000.ms"
OUTDIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0/cube_results"

echo "=== tclean CUBE MODE from Lustre on $(hostname), CPUs: $(nproc) ==="
echo "Start: $(date)"

python3 << 'PYEOF'
import time, os, shutil, sys
from casatasks import tclean, exportfits
import casatools

SUBBAND = "/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/subbands/subband_000.ms"
OUTDIR = "/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0/cube_results"
os.makedirs(OUTDIR, exist_ok=True)

msmd = casatools.msmetadata()
msmd.open(SUBBAND)
nchan = msmd.nchan(0)
freqs = msmd.chanfreqs(0)
msmd.close()
print(f"Subband: {nchan} channels, {freqs[0]/1e6:.3f}-{freqs[-1]/1e6:.3f} MHz", flush=True)

# Test 1: 20 channels cube
print("\n=== CUBE 20 channels, niter=0 ===", flush=True)
imagename = os.path.join(OUTDIR, "cube20")
for s in ['.image','.model','.pb','.psf','.residual','.sumwt','.mask']:
    if os.path.isdir(imagename+s): shutil.rmtree(imagename+s)
t0 = time.time()
tclean(vis=SUBBAND, imagename=imagename, specmode='cube',
       nchan=20, start=0, width=1, imsize=[512,512], cell='2arcsec',
       niter=0, weighting='briggs', robust=0.5, savemodel='none', pbcor=False)
t20 = time.time() - t0
print(f"  Time: {t20:.1f}s ({t20/60:.1f}min), per chan: {t20/20:.1f}s", flush=True)
for s in ['.image','.model','.pb','.psf','.residual','.sumwt','.mask']:
    if os.path.isdir(imagename+s): shutil.rmtree(imagename+s)

# Test 2: 100 channels cube
print("\n=== CUBE 100 channels, niter=0 ===", flush=True)
imagename = os.path.join(OUTDIR, "cube100")
for s in ['.image','.model','.pb','.psf','.residual','.sumwt','.mask']:
    if os.path.isdir(imagename+s): shutil.rmtree(imagename+s)
t0 = time.time()
tclean(vis=SUBBAND, imagename=imagename, specmode='cube',
       nchan=100, start=0, width=1, imsize=[512,512], cell='2arcsec',
       niter=0, weighting='briggs', robust=0.5, savemodel='none', pbcor=False)
t100 = time.time() - t0
print(f"  Time: {t100:.1f}s ({t100/60:.1f}min), per chan: {t100/100:.1f}s", flush=True)
for s in ['.image','.model','.pb','.psf','.residual','.sumwt','.mask']:
    if os.path.isdir(imagename+s): shutil.rmtree(imagename+s)

# Test 3: all 383 channels cube
print("\n=== CUBE ALL 383 channels, niter=0 ===", flush=True)
imagename = os.path.join(OUTDIR, "cube383")
for s in ['.image','.model','.pb','.psf','.residual','.sumwt','.mask']:
    if os.path.isdir(imagename+s): shutil.rmtree(imagename+s)
t0 = time.time()
tclean(vis=SUBBAND, imagename=imagename, specmode='cube',
       nchan=383, start=0, width=1, imsize=[512,512], cell='2arcsec',
       niter=0, weighting='briggs', robust=0.5, savemodel='none', pbcor=False)
t383 = time.time() - t0
print(f"  Time: {t383:.1f}s ({t383/60:.1f}min), per chan: {t383/383:.1f}s", flush=True)

# Don't clean up cube383 — we'll inspect it
print(f"\n=== SUMMARY ===", flush=True)
print(f"  20 chan: {t20:.0f}s ({t20/20:.1f}s/chan)", flush=True)
print(f"  100 chan: {t100:.0f}s ({t100/100:.1f}s/chan)", flush=True)
print(f"  383 chan: {t383:.0f}s ({t383/383:.1f}s/chan)", flush=True)
print(f"\n=== ALL DONE ===", flush=True)
PYEOF

echo "End: $(date)"
