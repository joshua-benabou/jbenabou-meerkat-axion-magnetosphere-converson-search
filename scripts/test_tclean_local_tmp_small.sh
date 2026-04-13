#!/bin/bash
#SBATCH --job-name=test_tmp_s
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_tclean_local_tmp_small.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_tclean_local_tmp_small.err
#
# Tests tclean from local /tmp:
# - Sequential dirty (niter=0) vs cleaned (niter=10000)
# - Parallel 4, 8, 14 workers (dirty)
#

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

SUBBAND="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0/multinode_0.ms"
LOCAL_SUBBAND="/tmp/subband_local.ms"
OUTDIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0/tclean_results"

echo "=== tclean test on $(hostname), CPUs: $(nproc) ==="
echo "Start: $(date)"

# Copy to local
echo "=== COPYING TO /tmp ==="
t0=$(date +%s)
if [ ! -d "$LOCAL_SUBBAND" ]; then
    cp -r "$SUBBAND" "$LOCAL_SUBBAND"
fi
echo "Copy time: $(($(date +%s) - t0))s"

python3 << 'PYEOF'
import time, os, shutil, sys
from multiprocessing import Pool
from casatasks import tclean, exportfits
import casatools

LOCAL = "/tmp/subband_local.ms"
OUTDIR = "/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0/tclean_results"
os.makedirs(OUTDIR, exist_ok=True)

msmd = casatools.msmetadata()
msmd.open(LOCAL)
nchan = msmd.nchan(0)
freqs = msmd.chanfreqs(0)
msmd.close()
print(f"Subband: {nchan} channels, {freqs[0]/1e6:.3f}-{freqs[-1]/1e6:.3f} MHz", flush=True)

def run_tclean(vis, imagename, niter):
    kwargs = dict(
        vis=vis, imagename=imagename, specmode='mfs',
        imsize=[512,512], cell='2arcsec', weighting='briggs',
        robust=0.5, savemodel='none', pbcor=False, niter=niter,
    )
    if niter > 0:
        kwargs['usemask'] = 'auto-multithresh'
        kwargs['sidelobethreshold'] = 2.0
        kwargs['noisethreshold'] = 5.0
    tclean(**kwargs)

def image_channel_dirty(chan):
    imagename = os.path.join(OUTDIR, f"dirty_{chan:04d}")
    fitsname = imagename + ".fits"
    for s in ['.image','.model','.pb','.psf','.residual','.sumwt','.mask']:
        if os.path.isdir(imagename+s): shutil.rmtree(imagename+s)
    if os.path.exists(fitsname): os.remove(fitsname)
    t0 = time.time()
    run_tclean(LOCAL, imagename, niter=0)
    t_tcl = time.time() - t0
    exportfits(imagename=imagename+'.image', fitsimage=fitsname, overwrite=True)
    t_tot = time.time() - t0
    for s in ['.image','.model','.pb','.psf','.residual','.sumwt','.mask']:
        if os.path.isdir(imagename+s): shutil.rmtree(imagename+s)
    return (chan, t_tcl, t_tot)

def image_channel_cleaned(chan):
    imagename = os.path.join(OUTDIR, f"clean_{chan:04d}")
    fitsname = imagename + ".fits"
    for s in ['.image','.model','.pb','.psf','.residual','.sumwt','.mask']:
        if os.path.isdir(imagename+s): shutil.rmtree(imagename+s)
    if os.path.exists(fitsname): os.remove(fitsname)
    t0 = time.time()
    run_tclean(LOCAL, imagename, niter=10000)
    t_tcl = time.time() - t0
    exportfits(imagename=imagename+'.image', fitsimage=fitsname, overwrite=True)
    t_tot = time.time() - t0
    for s in ['.image','.model','.pb','.psf','.residual','.sumwt','.mask']:
        if os.path.isdir(imagename+s): shutil.rmtree(imagename+s)
    return (chan, t_tcl, t_tot)

# === 1. Sequential dirty vs cleaned ===
print("\n=== SEQUENTIAL DIRTY (niter=0) ===", flush=True)
for c in [0, 100]:
    ch, t_tcl, t_tot = image_channel_dirty(c)
    print(f"  Chan {ch}: tclean={t_tcl:.1f}s total={t_tot:.1f}s", flush=True)

print("\n=== SEQUENTIAL CLEANED (niter=10000) ===", flush=True)
for c in [0, 100]:
    ch, t_tcl, t_tot = image_channel_cleaned(c)
    print(f"  Chan {ch}: tclean={t_tcl:.1f}s total={t_tot:.1f}s", flush=True)

# === 2. Parallel dirty from local /tmp ===
print("\n=== PARALLEL 4 workers DIRTY ===", flush=True)
t0 = time.time()
with Pool(4) as pool:
    res = pool.map(image_channel_dirty, [10,50,90,130])
wall = time.time() - t0
for ch,t_tcl,t_tot in res:
    print(f"  Chan {ch}: tclean={t_tcl:.1f}s total={t_tot:.1f}s", flush=True)
print(f"  Wall: {wall:.1f}s ({wall/60:.1f}min)", flush=True)

print("\n=== PARALLEL 8 workers DIRTY ===", flush=True)
t0 = time.time()
with Pool(8) as pool:
    res = pool.map(image_channel_dirty, list(range(8)))
wall = time.time() - t0
for ch,t_tcl,t_tot in res:
    print(f"  Chan {ch}: tclean={t_tcl:.1f}s total={t_tot:.1f}s", flush=True)
print(f"  Wall: {wall:.1f}s ({wall/60:.1f}min)", flush=True)

print("\n=== PARALLEL 14 workers DIRTY ===", flush=True)
t0 = time.time()
with Pool(14) as pool:
    res = pool.map(image_channel_dirty, list(range(14)))
wall = time.time() - t0
for ch,t_tcl,t_tot in res:
    print(f"  Chan {ch}: tclean={t_tcl:.1f}s total={t_tot:.1f}s", flush=True)
print(f"  Wall: {wall:.1f}s ({wall/60:.1f}min)", flush=True)

# === 3. Parallel cleaned from local /tmp ===
print("\n=== PARALLEL 4 workers CLEANED ===", flush=True)
t0 = time.time()
with Pool(4) as pool:
    res = pool.map(image_channel_cleaned, [10,50,90,130])
wall = time.time() - t0
for ch,t_tcl,t_tot in res:
    print(f"  Chan {ch}: tclean={t_tcl:.1f}s total={t_tot:.1f}s", flush=True)
print(f"  Wall: {wall:.1f}s ({wall/60:.1f}min)", flush=True)

print("\n=== ALL DONE ===", flush=True)
PYEOF

rm -rf /tmp/subband_local.ms
echo "End: $(date)"
