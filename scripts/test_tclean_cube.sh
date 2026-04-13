#!/bin/bash
#SBATCH --job-name=test_cube
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=90G
#SBATCH --time=10:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_tclean_cube.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_tclean_cube.err
#
# Test tclean cube mode: one call images all 383 channels.
# Also tests parallel=True for multi-core gridding.
# Copies subband to local /tmp first.
#

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

SUBBAND="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0/multinode_0.ms"
LOCAL="/tmp/subband_cube.ms"
OUTDIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0/cube_results"

echo "=== tclean CUBE MODE test on $(hostname), CPUs: $(nproc) ==="
echo "Start: $(date)"

# Copy to local
echo "=== COPYING TO /tmp ==="
t0=$(date +%s)
if [ ! -d "$LOCAL" ]; then
    cp -r "$SUBBAND" "$LOCAL"
fi
echo "Copy time: $(($(date +%s) - t0))s"

python3 << 'PYEOF'
import time, os, shutil, glob, sys
from casatasks import tclean, exportfits
import casatools

LOCAL = "/tmp/subband_cube.ms"
OUTDIR = "/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0/cube_results"
os.makedirs(OUTDIR, exist_ok=True)

msmd = casatools.msmetadata()
msmd.open(LOCAL)
nchan = msmd.nchan(0)
freqs = msmd.chanfreqs(0)
msmd.close()
print(f"Subband: {nchan} channels, {freqs[0]/1e6:.3f}-{freqs[-1]/1e6:.3f} MHz", flush=True)

# Test 1: Cube mode, small subset (20 channels), niter=0
print("\n=== CUBE MODE: 20 channels, niter=0 ===", flush=True)
imagename = os.path.join(OUTDIR, "cube_20chan_dirty")
for s in ['.image','.model','.pb','.psf','.residual','.sumwt','.mask']:
    if os.path.isdir(imagename+s): shutil.rmtree(imagename+s)

t0 = time.time()
tclean(
    vis=LOCAL,
    imagename=imagename,
    specmode='cube',
    nchan=20,
    start=0,
    width=1,
    imsize=[512,512],
    cell='2arcsec',
    niter=0,
    weighting='briggs',
    robust=0.5,
    savemodel='none',
    pbcor=False,
)
t_20 = time.time() - t0
print(f"  Time: {t_20:.1f}s ({t_20/60:.1f} min)", flush=True)
print(f"  Per channel: {t_20/20:.1f}s", flush=True)

# Cleanup
for s in ['.image','.model','.pb','.psf','.residual','.sumwt','.mask']:
    if os.path.isdir(imagename+s): shutil.rmtree(imagename+s)

# Test 2: Cube mode, all 383 channels, niter=0
print("\n=== CUBE MODE: ALL 383 channels, niter=0 ===", flush=True)
imagename = os.path.join(OUTDIR, "cube_383chan_dirty")
for s in ['.image','.model','.pb','.psf','.residual','.sumwt','.mask']:
    if os.path.isdir(imagename+s): shutil.rmtree(imagename+s)

t0 = time.time()
tclean(
    vis=LOCAL,
    imagename=imagename,
    specmode='cube',
    nchan=383,
    start=0,
    width=1,
    imsize=[512,512],
    cell='2arcsec',
    niter=0,
    weighting='briggs',
    robust=0.5,
    savemodel='none',
    pbcor=False,
)
t_383 = time.time() - t0
print(f"  Time: {t_383:.1f}s ({t_383/60:.1f} min)", flush=True)
print(f"  Per channel: {t_383/383:.1f}s", flush=True)

# Check output
cube_image = imagename + '.image'
if os.path.isdir(cube_image):
    from casatools import image as imagetool
    ia = imagetool()
    ia.open(cube_image)
    shape = ia.shape()
    ia.close()
    print(f"  Cube shape: {shape}", flush=True)

# Test 3: Try parallel=True if available
print("\n=== CUBE MODE: 20 channels, niter=0, parallel=True ===", flush=True)
imagename = os.path.join(OUTDIR, "cube_20chan_parallel")
for s in ['.image','.model','.pb','.psf','.residual','.sumwt','.mask']:
    if os.path.isdir(imagename+s): shutil.rmtree(imagename+s)

t0 = time.time()
try:
    tclean(
        vis=LOCAL,
        imagename=imagename,
        specmode='cube',
        nchan=20,
        start=0,
        width=1,
        imsize=[512,512],
        cell='2arcsec',
        niter=0,
        weighting='briggs',
        robust=0.5,
        savemodel='none',
        pbcor=False,
        parallel=True,
    )
    t_par = time.time() - t0
    print(f"  Time: {t_par:.1f}s ({t_par/60:.1f} min)", flush=True)
    print(f"  Speedup vs non-parallel: {t_20/t_par:.1f}x", flush=True)
except Exception as e:
    print(f"  parallel=True failed: {e}", flush=True)
    print(f"  (May need mpicasa for parallelism)", flush=True)

# Cleanup all
for s in ['.image','.model','.pb','.psf','.residual','.sumwt','.mask']:
    for prefix in ['cube_20chan_dirty','cube_383chan_dirty','cube_20chan_parallel']:
        p = os.path.join(OUTDIR, prefix+s)
        if os.path.isdir(p): shutil.rmtree(p)

print(f"\n=== SUMMARY ===", flush=True)
print(f"  20-chan cube: {t_20:.0f}s ({t_20/20:.1f}s/chan)", flush=True)
print(f"  383-chan cube: {t_383:.0f}s ({t_383/383:.1f}s/chan)", flush=True)
print(f"  383 × single-chan MFS (estimated): {383*t_20/20:.0f}s ({383*t_20/20/60:.0f} min)", flush=True)
print(f"  Cube speedup: {383*t_20/20/t_383:.1f}x", flush=True)
print(f"\n=== ALL DONE ===", flush=True)
PYEOF

rm -rf /tmp/subband_cube.ms
echo "End: $(date)"
