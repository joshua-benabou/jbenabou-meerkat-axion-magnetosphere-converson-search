#!/usr/bin/env python3
"""
Phase 3b: Cleaned channel imaging using tclean cube mode.

Strategy: measure peak flux from existing dirty images, then clean with
threshold = 10% of peak. No mask (the GC fills much of the field).

This replaces the broken auto-multithresh approach, which set thresholds
based on full-image RMS (dominated by sidelobes) and produced zero clean
components. The 10% peak threshold was validated in debug_cleaning.py:
- Reduces corner RMS from 40 -> 26 mJy (34% improvement)
- 10,781 model pixels, 14,987 mJy model flux
- 1%, 5%, and 10% peak all converge to the same solution

Usage:
    Called by phase3_clean_submit.sh via SLURM array job.
    SLURM_ARRAY_TASK_ID selects which subband to process.
"""

import os
import sys
import time
import shutil
import glob
import numpy as np
import casatools
from casatasks import tclean, exportfits, imsubimage, imstat

# === Configuration ===
SUBBANDS_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/subbands'
)
IMAGES_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/images'
)

IMSIZE = [512, 512]
CELL = '2arcsec'
NITER = 100000
WEIGHTING = 'briggs'
ROBUST = 0.5
SAVEMODEL = 'none'
THRESHOLD_FRAC = 0.10  # 10% of peak flux


def measure_peak_from_dirty(subband_idx):
    """Measure peak flux from existing dirty FITS files for this subband.

    Samples ~10 channels spread across the subband (avoiding channel 0
    which has a known edge artifact) and returns the median peak.
    """
    dirty_dir = os.path.join(IMAGES_DIR, f'subband_{subband_idx:03d}')
    fits_files = sorted(glob.glob(os.path.join(dirty_dir, 'chan_*.fits')))

    if not fits_files:
        return None

    # Sample channels spread across the subband
    n = len(fits_files)
    # Skip channel 0 (edge artifact), sample ~10 evenly spaced
    indices = np.linspace(max(1, n // 20), n - 1, min(10, n)).astype(int)

    from astropy.io import fits as pyfits
    peaks = []
    for idx in indices:
        with pyfits.open(fits_files[idx]) as hdul:
            data = hdul[0].data
            while data.ndim > 2:
                data = data[0]
            peaks.append(np.nanmax(data))

    median_peak = np.median(peaks)
    print(f"  Peak from dirty images: median={median_peak*1e3:.3f} mJy "
          f"(sampled {len(peaks)} channels)", flush=True)
    return median_peak


def main():
    subband_idx = int(os.environ.get('SLURM_ARRAY_TASK_ID',
                                      sys.argv[1] if len(sys.argv) > 1 else '0'))

    input_ms = os.path.join(SUBBANDS_DIR, f'subband_{subband_idx:03d}.ms')
    output_dir = os.path.join(IMAGES_DIR, f'subband_{subband_idx:03d}', 'cleaned')

    if not os.path.exists(input_ms):
        print(f"Subband {subband_idx}: {input_ms} does not exist, skipping.", flush=True)
        return

    # Get channel info
    msmd = casatools.msmetadata()
    msmd.open(input_ms)
    nchan = msmd.nchan(0)
    freqs = msmd.chanfreqs(0)
    msmd.close()

    # Measure peak from existing dirty images to set threshold
    peak_jy = measure_peak_from_dirty(subband_idx)
    if peak_jy is None or peak_jy <= 0:
        print(f"Subband {subband_idx}: could not measure peak from dirty images. "
              f"Run dirty imaging first.", flush=True)
        return

    threshold_jy = THRESHOLD_FRAC * peak_jy
    threshold_mjy = threshold_jy * 1e3
    print(f"  Cleaning threshold: {THRESHOLD_FRAC*100:.0f}% of peak = "
          f"{threshold_mjy:.3f} mJy", flush=True)

    # Remove old broken cleaned FITS if they exist
    if os.path.isdir(output_dir):
        old_fits = glob.glob(os.path.join(output_dir, '*.fits'))
        if old_fits:
            print(f"  Removing {len(old_fits)} old cleaned FITS from broken run...",
                  flush=True)
            for f in old_fits:
                os.remove(f)

    os.makedirs(output_dir, exist_ok=True)

    print(f"=== Phase 3b: Cleaned imaging subband {subband_idx} ===", flush=True)
    print(f"  Input: {input_ms}", flush=True)
    print(f"  Channels: {nchan}", flush=True)
    print(f"  Freq range: {freqs[0]/1e6:.3f}-{freqs[-1]/1e6:.3f} MHz", flush=True)
    print(f"  niter: {NITER}", flush=True)
    print(f"  threshold: {threshold_mjy:.3f} mJy ({THRESHOLD_FRAC*100:.0f}% of peak)", flush=True)
    print(f"  Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)

    # Step 1: tclean cube mode with cleaning
    imagename = os.path.join(output_dir, f'subband_{subband_idx:03d}_clean_cube')

    # Clean up any prior incomplete CASA products
    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = imagename + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)

    t0 = time.time()
    tclean(
        vis=input_ms,
        imagename=imagename,
        specmode='cube',
        nchan=nchan,
        start=0,
        width=1,
        imsize=IMSIZE,
        cell=CELL,
        niter=NITER,
        weighting=WEIGHTING,
        robust=ROBUST,
        savemodel=SAVEMODEL,
        pbcor=False,
        threshold=f'{threshold_mjy}mJy',
    )
    t_tclean = time.time() - t0
    print(f"  tclean clean cube done: {t_tclean:.1f}s ({t_tclean/60:.1f} min)", flush=True)

    # Quick sanity check: model flux from the cube
    model_image = imagename + '.model'
    if os.path.isdir(model_image):
        model_stats = imstat(model_image)
        print(f"  Model flux (total): {model_stats['sum'][0]*1e3:.3f} mJy", flush=True)
        print(f"  Model max:          {model_stats['max'][0]*1e3:.3f} mJy/beam", flush=True)

    # Step 2: Export each channel to individual FITS
    cube_image = imagename + '.image'
    if not os.path.isdir(cube_image):
        print(f"  ERROR: {cube_image} not found after tclean!", flush=True)
        return

    print(f"  Exporting {nchan} channels to FITS...", flush=True)
    t0_export = time.time()
    n_exported = 0
    for chan in range(nchan):
        freq_mhz = freqs[chan] / 1e6
        fitsname = os.path.join(output_dir, f'chan_{chan:04d}_{freq_mhz:.3f}MHz.fits')
        if os.path.exists(fitsname):
            continue
        chan_image = os.path.join(output_dir, f'_tmp_chan_{chan:04d}.image')
        if os.path.isdir(chan_image):
            shutil.rmtree(chan_image)
        imsubimage(
            imagename=cube_image,
            outfile=chan_image,
            chans=str(chan),
        )
        exportfits(
            imagename=chan_image,
            fitsimage=fitsname,
            overwrite=True,
        )
        shutil.rmtree(chan_image)
        n_exported += 1
        if n_exported % 50 == 0:
            print(f"    Exported {n_exported}/{nchan} (chan {chan})", flush=True)

    t_export = time.time() - t0_export
    print(f"  Export done: {n_exported} new files in "
          f"{t_export:.1f}s ({t_export/60:.1f} min)", flush=True)

    # Verify all FITS files exist before cleaning up CASA products
    n_fits = len([f for f in os.listdir(output_dir) if f.endswith('.fits')])
    if n_fits >= nchan:
        print(f"  All {n_fits} FITS files present, cleaning up CASA products", flush=True)
        for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
            p = imagename + suffix
            if os.path.isdir(p):
                shutil.rmtree(p)
    else:
        print(f"  WARNING: Only {n_fits}/{nchan} FITS -- keeping CASA products "
              f"for recovery", flush=True)

    t_total = time.time() - t0
    print(f"  Total: {t_total:.1f}s ({t_total/60:.1f} min), {n_fits} FITS files", flush=True)
    print(f"  End time: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)


if __name__ == '__main__':
    main()
