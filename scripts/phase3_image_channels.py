#!/usr/bin/env python3
"""
Phase 3: Image every channel in a subband using tclean.

Produces FITS images per channel — both dirty (niter=0) and optionally cleaned
(niter>0). Uses specmode='mfs' with single-channel spw selection to image each
channel independently. Parallelized across channels within each node using
multiprocessing.

Copies the subband MS to local /tmp before imaging to avoid Lustre I/O
contention between parallel workers. Falls back to Lustre if /tmp is too small.

Usage:
    Called by phase3_submit.sh via SLURM array job.
    SLURM_ARRAY_TASK_ID selects which subband to process.
    Set MAKE_CLEANED=1 env var to also produce cleaned images (slower).
"""

import os
import sys
import time
import glob
import shutil
import subprocess
from multiprocessing import Pool, cpu_count
import casatools
from casatasks import tclean, exportfits

# === Configuration ===
SUBBANDS_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/subbands'
)
IMAGES_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/images'
)

# tclean parameters
IMSIZE = [512, 512]       # Start small; scale to 4096 after validation
CELL = '2arcsec'
WEIGHTING = 'briggs'
ROBUST = 0.5
SAVEMODEL = 'none'        # CRITICAL: never save model back to MS

# Cleaning parameters (used when MAKE_CLEANED=1)
CLEAN_NITER = 10000
CLEAN_THRESHOLD = '3sigma'
CLEAN_AUTOMASK = 5.0

# Whether to produce cleaned images in addition to dirty
MAKE_CLEANED = os.environ.get('MAKE_CLEANED', '0') == '1'

# Parallelism: use most of the node's cores, leaving a couple free
NWORKERS = max(1, int(os.environ.get('SLURM_CPUS_PER_TASK', cpu_count())) - 2)


def copy_to_tmp(src_ms):
    """Copy subband MS to /tmp for fast local I/O. Returns local path or src_ms on failure."""
    local_ms = f'/tmp/subband_{os.getpid()}.ms'

    # Check /tmp has enough space (need ~160GB free)
    try:
        stat = os.statvfs('/tmp')
        free_gb = stat.f_bavail * stat.f_frsize / 1e9
        if free_gb < 160:
            print(f"  /tmp has only {free_gb:.0f}GB free, staying on Lustre", flush=True)
            return src_ms, False
    except Exception:
        return src_ms, False

    print(f"  Copying subband to /tmp ({free_gb:.0f}GB free)...", flush=True)
    t0 = time.time()
    try:
        subprocess.run(['cp', '-r', src_ms, local_ms], check=True)
        elapsed = time.time() - t0
        print(f"  Copy done: {elapsed:.0f}s ({elapsed/60:.1f} min)", flush=True)
        return local_ms, True
    except Exception as e:
        print(f"  Copy failed ({e}), staying on Lustre", flush=True)
        if os.path.exists(local_ms):
            shutil.rmtree(local_ms, ignore_errors=True)
        return src_ms, False


def _run_tclean_and_export(input_ms, imagename, fitsname, spw, niter):
    """Run tclean on a single channel and export to FITS, cleaning up CASA products."""
    tclean_kwargs = dict(
        vis=input_ms,
        imagename=imagename,
        specmode='mfs',
        spw=spw,
        imsize=IMSIZE,
        cell=CELL,
        weighting=WEIGHTING,
        robust=ROBUST,
        savemodel=SAVEMODEL,
        pbcor=False,
        niter=niter,
    )
    if niter > 0:
        tclean_kwargs['usemask'] = 'auto-multithresh'
        tclean_kwargs['sidelobethreshold'] = 2.0
        tclean_kwargs['noisethreshold'] = CLEAN_AUTOMASK

    tclean(**tclean_kwargs)

    casa_image = imagename + '.image'
    if os.path.exists(casa_image):
        exportfits(imagename=casa_image, fitsimage=fitsname, overwrite=True)

    # Clean up CASA intermediate products
    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        path = imagename + suffix
        if os.path.isdir(path):
            shutil.rmtree(path)


def image_channel(args):
    """Image a single channel. Designed to run in a worker process."""
    input_ms, output_dir, chan, freq_mhz, make_cleaned = args

    spw = f'0:{chan}'
    results = []

    # --- Dirty image ---
    dirty_name = os.path.join(output_dir, 'dirty', f'chan_{chan:04d}_{freq_mhz:.3f}MHz')
    dirty_fits = dirty_name + '.fits'

    t0 = time.time()
    try:
        if not os.path.exists(dirty_fits):
            os.makedirs(os.path.dirname(dirty_name), exist_ok=True)
            _run_tclean_and_export(input_ms, dirty_name, dirty_fits, spw, niter=0)
        results.append(('dirty', chan, time.time() - t0, 'ok'))
    except Exception as e:
        results.append(('dirty', chan, time.time() - t0, f'ERROR: {e}'))

    # --- Cleaned image ---
    if make_cleaned:
        clean_name = os.path.join(output_dir, 'cleaned', f'chan_{chan:04d}_{freq_mhz:.3f}MHz')
        clean_fits = clean_name + '.fits'

        t0 = time.time()
        try:
            if not os.path.exists(clean_fits):
                os.makedirs(os.path.dirname(clean_name), exist_ok=True)
                _run_tclean_and_export(
                    input_ms, clean_name, clean_fits, spw, niter=CLEAN_NITER,
                )
            results.append(('cleaned', chan, time.time() - t0, 'ok'))
        except Exception as e:
            results.append(('cleaned', chan, time.time() - t0, f'ERROR: {e}'))

    return results


def main():
    subband_idx = int(os.environ.get('SLURM_ARRAY_TASK_ID', sys.argv[1] if len(sys.argv) > 1 else '0'))

    # Use contsub MS if available, otherwise raw subband
    contsub_ms = os.path.join(SUBBANDS_DIR, f'subband_{subband_idx:03d}.ms.contsub')
    raw_ms = os.path.join(SUBBANDS_DIR, f'subband_{subband_idx:03d}.ms')

    if os.path.exists(contsub_ms):
        src_ms = contsub_ms
    elif os.path.exists(raw_ms):
        src_ms = raw_ms
        print(f"  WARNING: Using raw MS (no continuum subtraction)", flush=True)
    else:
        print(f"Subband {subband_idx}: no input MS found, skipping.", flush=True)
        return

    output_dir = os.path.join(IMAGES_DIR, f'subband_{subband_idx:03d}')
    os.makedirs(output_dir, exist_ok=True)

    # Copy to local /tmp for fast I/O
    input_ms, is_local = copy_to_tmp(src_ms)

    # Get number of channels in this subband
    msmd = casatools.msmetadata()
    msmd.open(input_ms)
    nchan = msmd.nchan(0)
    freqs = msmd.chanfreqs(0)
    msmd.close()

    print(f"=== Phase 3: Imaging subband {subband_idx} ===", flush=True)
    print(f"  Input: {input_ms} ({'local /tmp' if is_local else 'Lustre'})", flush=True)
    print(f"  Channels: {nchan}", flush=True)
    print(f"  Freq range: {freqs[0]/1e6:.3f}-{freqs[-1]/1e6:.3f} MHz", flush=True)
    print(f"  Image size: {IMSIZE}, cell: {CELL}", flush=True)
    print(f"  Dirty: always | Cleaned: {'yes' if MAKE_CLEANED else 'no'}", flush=True)
    print(f"  Workers: {NWORKERS}", flush=True)
    print(f"  Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)

    # Build work list
    work = []
    for chan in range(nchan):
        freq_mhz = freqs[chan] / 1e6
        work.append((input_ms, output_dir, chan, freq_mhz, MAKE_CLEANED))

    t0_total = time.time()

    # Process channels in parallel
    completed = 0
    errors = 0
    with Pool(processes=NWORKERS) as pool:
        for result_list in pool.imap_unordered(image_channel, work):
            for img_type, chan, elapsed, status in result_list:
                if status == 'ok':
                    completed += 1
                elif status != 'skipped':
                    errors += 1
                    print(f"  {img_type} chan {chan}: {status}", flush=True)

            total_done = completed + errors
            if total_done % 50 == 0 or total_done >= nchan:
                rate = completed / max(time.time() - t0_total, 1) * 60
                print(f"  Progress: {completed}/{nchan} done, {errors} errors "
                      f"({rate:.1f} images/min)", flush=True)

    # Cleanup local copy
    if is_local and os.path.exists(input_ms):
        print(f"  Cleaning up /tmp copy...", flush=True)
        shutil.rmtree(input_ms, ignore_errors=True)

    total_elapsed = time.time() - t0_total
    n_dirty = len(glob.glob(os.path.join(output_dir, 'dirty', '*.fits')))
    n_clean = len(glob.glob(os.path.join(output_dir, 'cleaned', '*.fits')))
    print(f"  Completed: {n_dirty} dirty, {n_clean} cleaned FITS images, {errors} errors", flush=True)
    print(f"  Total elapsed: {total_elapsed:.1f}s ({total_elapsed/60:.1f} min)", flush=True)
    print(f"  End time: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)


if __name__ == '__main__':
    main()
