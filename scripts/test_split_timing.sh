#!/bin/bash
#SBATCH --job-name=test_split
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=90G
#SBATCH --time=04:00:00
#SBATCH --output=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_split_timing.out
#SBATCH --error=/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/test_split_timing.err
#
# Test split timing on a real compute node with 32 cores.
# Runs 1-chan, 10-chan, and 383-chan splits sequentially, then
# a parallel test with 4 simultaneous 383-chan splits.
#

SCRIPT_DIR="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
TEST_DIR="${PROJECT_DIR}/test_phase0"
MS_PATH="/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/SCI-20210212-SS-01/17TBdataset/scratch/kat/1621534878_20231022T19_11_29/1621534878_sdp_l0.ms"

mkdir -p "${SCRIPT_DIR}/logs"

source /global/home/users/jbenabou/opt/miniforge3/etc/profile.d/conda.sh
conda activate /global/home/users/osning/.conda/envs/casa_cookbook

echo "=== Split Timing Tests on $(hostname) ==="
echo "CPUs: $(nproc), RAM: $(free -h | awk '/Mem:/{print $2}')"
echo "Start: $(date)"
echo ""

python3 << 'PYEOF'
import time, os, shutil
from casatasks import split
import casatools

MS = "/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/SCI-20210212-SS-01/17TBdataset/scratch/kat/1621534878_20231022T19_11_29/1621534878_sdp_l0.ms"
TEST = "/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_phase0"

def timed_split(name, spw, nchan_label):
    out = os.path.join(TEST, f"slurm_{name}.ms")
    if os.path.exists(out):
        shutil.rmtree(out)

    print(f"\n--- {name} ({nchan_label}) ---")
    t0 = time.time()
    split(vis=MS, outputvis=out, field='SgrA*', spw=spw, datacolumn='data')
    elapsed = time.time() - t0

    # Verify
    msmd = casatools.msmetadata()
    msmd.open(out)
    nc = msmd.nchan(0)
    nr = msmd.nrows()
    msmd.close()

    size_mb = sum(os.path.getsize(os.path.join(dp, f)) for dp, dn, fn in os.walk(out) for f in fn) / 1e6
    print(f"  Time: {elapsed:.1f}s ({elapsed/60:.1f} min)")
    print(f"  Size: {size_mb:.1f} MB")
    print(f"  Channels: {nc}, Rows: {nr:.0f}")
    return elapsed

# Sequential tests
print("=== SEQUENTIAL TESTS ===")
t1 = timed_split("1chan", "0:17000", "1 channel")
t10 = timed_split("10chan", "0:16997~17006", "10 channels")
t383 = timed_split("383chan", "0:16997~17379", "383 channels")

print(f"\n=== SEQUENTIAL SUMMARY ===")
print(f"  1 chan:   {t1:.1f}s")
print(f"  10 chan:  {t10:.1f}s")
print(f"  383 chan: {t383:.1f}s")
print(f"  Scaling: 10/1={t10/t1:.2f}x, 383/1={t383/t1:.2f}x")

# Parallel test: 4 splits at once
print(f"\n=== PARALLEL TEST (4 x 383-chan splits) ===")
from multiprocessing import Pool

def parallel_split(args):
    idx, spw_start = args
    name = f"parallel_{idx}"
    spw_end = spw_start + 382
    out = os.path.join(TEST, f"slurm_{name}.ms")
    if os.path.exists(out):
        shutil.rmtree(out)
    t0 = time.time()
    split(vis=MS, outputvis=out, field='SgrA*', spw=f'0:{spw_start}~{spw_end}', datacolumn='data')
    return (idx, time.time() - t0)

t0_parallel = time.time()
with Pool(4) as pool:
    results = pool.map(parallel_split, [(0, 0), (1, 383), (2, 766), (3, 1149)])
parallel_wall = time.time() - t0_parallel

print(f"  Wall time: {parallel_wall:.1f}s ({parallel_wall/60:.1f} min)")
for idx, elapsed in results:
    print(f"  Job {idx}: {elapsed:.1f}s")
print(f"  Speedup vs serial: {t383*4/parallel_wall:.2f}x")

print(f"\n=== ALL TESTS COMPLETE ===")
print(f"End: {time.strftime('%Y-%m-%d %H:%M:%S')}")
PYEOF
