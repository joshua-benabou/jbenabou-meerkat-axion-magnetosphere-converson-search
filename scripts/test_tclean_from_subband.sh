#!/bin/bash
#SBATCH --job-name=test_sub_tcl
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=90G
#SBATCH --time=02:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_tclean_from_subband.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_tclean_from_subband.err
#
# Test tclean speed on a split subband (~154 GB) vs the monolithic MS (17 TB).
# Images 4 channels sequentially, then 4 in parallel.
#

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

SUBBAND="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0/multinode_0.ms"
OUTDIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0/subband_tclean"

echo "=== tclean from split subband test on $(hostname) ==="
echo "Subband: $SUBBAND"
echo "CPUs: $(nproc)"
echo "Start: $(date)"

python3 << 'PYEOF'
import time, os, shutil
from multiprocessing import Pool
from casatasks import tclean, exportfits
import casatools

SUBBAND = "/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0/multinode_0.ms"
OUTDIR = "/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0/subband_tclean"
os.makedirs(OUTDIR, exist_ok=True)

# Verify subband
msmd = casatools.msmetadata()
msmd.open(SUBBAND)
nchan = msmd.nchan(0)
freqs = msmd.chanfreqs(0)
msmd.close()
print(f"Subband: {nchan} channels, {freqs[0]/1e6:.3f}-{freqs[-1]/1e6:.3f} MHz")

def image_one_channel(chan):
    freq_mhz = freqs[chan] / 1e6
    imagename = os.path.join(OUTDIR, f"sub_chan_{chan:04d}_{freq_mhz:.3f}MHz")
    fitsname = imagename + ".fits"

    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = imagename + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)
    if os.path.exists(fitsname):
        os.remove(fitsname)

    t0 = time.time()
    tclean(
        vis=SUBBAND,
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

# Sequential: image 4 channels one at a time
print("\n=== SEQUENTIAL (4 channels) ===")
seq_times = []
for chan in [0, 100, 200, 300]:
    chan_idx, t_tcl, t_tot = image_one_channel(chan)
    print(f"  Chan {chan_idx}: tclean={t_tcl:.1f}s total={t_tot:.1f}s")
    seq_times.append(t_tot)
print(f"  Avg: {sum(seq_times)/len(seq_times):.1f}s per channel")

# Parallel: image 8 channels simultaneously
print("\n=== PARALLEL (8 channels, 8 workers) ===")
t0 = time.time()
with Pool(8) as pool:
    results = pool.map(image_one_channel, [10, 50, 90, 130, 170, 210, 250, 350])
parallel_wall = time.time() - t0
for chan, t_tcl, t_tot in results:
    print(f"  Chan {chan}: tclean={t_tcl:.1f}s total={t_tot:.1f}s")
print(f"  Wall time for 8 channels: {parallel_wall:.1f}s ({parallel_wall/60:.1f} min)")
print(f"  Speedup vs 8 sequential: {sum(seq_times)/len(seq_times)*8/parallel_wall:.1f}x")

# Parallel: 30 channels
print("\n=== PARALLEL (30 channels, 30 workers) ===")
t0 = time.time()
with Pool(30) as pool:
    results = pool.map(image_one_channel, list(range(30)))
parallel_wall_30 = time.time() - t0
print(f"  Wall time for 30 channels: {parallel_wall_30:.1f}s ({parallel_wall_30/60:.1f} min)")
avg_seq = sum(seq_times)/len(seq_times)
print(f"  Speedup vs 30 sequential: {avg_seq*30/parallel_wall_30:.1f}x")
print(f"  Estimated time for 383 channels at 30 workers: {383/30*parallel_wall_30/30:.1f}s ({383/30*parallel_wall_30/30/60:.1f} min)")

print(f"\n=== ALL DONE ===")
PYEOF

echo "End: $(date)"
