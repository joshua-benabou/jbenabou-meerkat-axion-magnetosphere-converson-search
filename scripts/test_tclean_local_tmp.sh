#!/bin/bash
#SBATCH --job-name=test_tmp_tcl
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=90G
#SBATCH --time=04:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_tclean_local_tmp.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_tclean_local_tmp.err
#
# Test tclean with subband copied to local /tmp SSD.
# This avoids Lustre I/O contention for parallel workers.
#

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

SUBBAND="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0/multinode_0.ms"
LOCAL_SUBBAND="/tmp/subband_local.ms"
OUTDIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0/subband_tclean_local"

echo "=== tclean from LOCAL /tmp test on $(hostname) ==="
echo "CPUs: $(nproc), /tmp free: $(df -h /tmp | tail -1 | awk '{print $4}')"
echo "Start: $(date)"

# Step 1: Copy subband to local /tmp
echo ""
echo "=== COPYING SUBBAND TO /tmp ==="
t0_copy=$(date +%s)
cp -r "$SUBBAND" "$LOCAL_SUBBAND"
t1_copy=$(date +%s)
copy_time=$((t1_copy - t0_copy))
echo "Copy time: ${copy_time}s ($(echo "scale=1; $copy_time/60" | bc) min)"
echo "/tmp used after copy: $(du -sh /tmp/subband_local.ms | cut -f1)"

python3 << 'PYEOF'
import time, os, shutil
from multiprocessing import Pool
from casatasks import tclean, exportfits
import casatools

LOCAL = "/tmp/subband_local.ms"
OUTDIR = "/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0/subband_tclean_local"
os.makedirs(OUTDIR, exist_ok=True)

msmd = casatools.msmetadata()
msmd.open(LOCAL)
nchan = msmd.nchan(0)
freqs = msmd.chanfreqs(0)
msmd.close()
print(f"Subband: {nchan} channels, {freqs[0]/1e6:.3f}-{freqs[-1]/1e6:.3f} MHz")

def image_one_channel(chan):
    freq_mhz = freqs[chan] / 1e6
    imagename = os.path.join(OUTDIR, f"local_chan_{chan:04d}_{freq_mhz:.3f}MHz")
    fitsname = imagename + ".fits"

    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = imagename + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)
    if os.path.exists(fitsname):
        os.remove(fitsname)

    t0 = time.time()
    tclean(
        vis=LOCAL,
        imagename=imagename,
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

    exportfits(imagename=imagename + '.image', fitsimage=fitsname, overwrite=True)
    t_total = time.time() - t0

    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = imagename + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)

    return (chan, t_tclean, t_total)

# Sequential: 2 channels (quick baseline from local disk)
print("\n=== SEQUENTIAL FROM LOCAL /tmp (2 channels) ===")
for chan in [0, 100]:
    ch, t_tcl, t_tot = image_one_channel(chan)
    print(f"  Chan {ch}: tclean={t_tcl:.1f}s total={t_tot:.1f}s")

# Parallel: 8 workers
print("\n=== PARALLEL 8 workers FROM LOCAL /tmp ===")
t0 = time.time()
with Pool(8) as pool:
    results = pool.map(image_one_channel, [10, 50, 90, 130, 170, 210, 250, 350])
wall_8 = time.time() - t0
for ch, t_tcl, t_tot in results:
    print(f"  Chan {ch}: tclean={t_tcl:.1f}s total={t_tot:.1f}s")
print(f"  Wall time: {wall_8:.1f}s ({wall_8/60:.1f} min)")

# Parallel: 16 workers
print("\n=== PARALLEL 16 workers FROM LOCAL /tmp ===")
t0 = time.time()
with Pool(16) as pool:
    results = pool.map(image_one_channel, list(range(16)))
wall_16 = time.time() - t0
for ch, t_tcl, t_tot in results:
    print(f"  Chan {ch}: tclean={t_tcl:.1f}s total={t_tot:.1f}s")
print(f"  Wall time: {wall_16:.1f}s ({wall_16/60:.1f} min)")

# Parallel: 30 workers
print("\n=== PARALLEL 30 workers FROM LOCAL /tmp ===")
t0 = time.time()
with Pool(30) as pool:
    results = pool.map(image_one_channel, list(range(30)))
wall_30 = time.time() - t0
print(f"  Wall time: {wall_30:.1f}s ({wall_30/60:.1f} min)")
print(f"  Estimated 383 chans at 30 workers: {383/30 * wall_30/30 * 30:.0f}s ({383/30 * wall_30/30 * 30 / 60:.1f} min)")

print(f"\n=== SUMMARY ===")
print(f"  Sequential: ~{(results[0][2] + results[1][2])/2 if len(results) >= 2 else 0:.0f}s/chan from local /tmp")
print(f"  8 parallel wall: {wall_8:.0f}s for 8 chans")
print(f"  16 parallel wall: {wall_16:.0f}s for 16 chans")
print(f"  30 parallel wall: {wall_30:.0f}s for 30 chans")

print(f"\n=== ALL DONE ===")
PYEOF

# Cleanup local copy
rm -rf /tmp/subband_local.ms

echo "End: $(date)"
