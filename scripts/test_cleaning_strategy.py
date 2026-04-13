#!/usr/bin/env python3
"""
Test cleaning strategies on a single subband to find one that works.

The previous attempt (auto-multithresh, noisethreshold=5.0) produced zero clean
components because per-channel SNR (~2.7) is far below the threshold.

This script tests the continuum-mask approach:
1. Make a continuum image (average all channels) → high SNR
2. Threshold it to create a clean mask
3. Clean a few test channels using that mask + a fixed threshold
4. Compare dirty vs cleaned images and save diagnostic plots

Usage:
    python test_cleaning_strategy.py [subband_idx]
    Default: subband 30 (mid-band, low RFI, ~1156 MHz)
"""

import os
import sys
import time
import shutil
import numpy as np
import casatools
from casatasks import tclean, exportfits, imsubimage, imstat, makemask, immath

# === Configuration ===
PROJECT_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project'
)
SUBBANDS_DIR = os.path.join(PROJECT_DIR, 'subbands')
TEST_DIR = os.path.join(PROJECT_DIR, 'test_cleaning')

IMSIZE = [512, 512]
CELL = '2arcsec'
WEIGHTING = 'briggs'
ROBUST = 0.5
SAVEMODEL = 'none'

# Test channels: spread across subband, avoiding edge (chan 0 has artifact)
TEST_CHANNELS = [10, 50, 100, 150, 200, 300]


def make_continuum_image(input_ms, output_dir, nchan):
    """Make a high-SNR continuum image by averaging all channels."""
    imagename = os.path.join(output_dir, 'continuum')

    # Clean up prior products
    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = imagename + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)

    print("  Making continuum image (niter=5000, all channels averaged)...", flush=True)
    t0 = time.time()
    tclean(
        vis=input_ms,
        imagename=imagename,
        specmode='mfs',          # Multi-frequency synthesis = average all channels
        imsize=IMSIZE,
        cell=CELL,
        niter=5000,
        weighting=WEIGHTING,
        robust=ROBUST,
        savemodel=SAVEMODEL,
        pbcor=False,
        threshold='0.0mJy',      # Let it clean deeply
        usemask='auto-multithresh',  # auto-mask works here because SNR is high
        sidelobethreshold=2.0,
        noisethreshold=4.0,
        lownoisethreshold=1.5,
        minbeamfrac=0.3,
    )
    dt = time.time() - t0
    print(f"  Continuum image done: {dt:.1f}s ({dt/60:.1f} min)", flush=True)

    # Get stats
    stats = imstat(imagename + '.image')
    print(f"  Continuum peak: {stats['max'][0]*1e3:.3f} mJy/beam", flush=True)
    print(f"  Continuum RMS:  {stats['rms'][0]*1e3:.3f} mJy/beam", flush=True)
    print(f"  Continuum SNR:  {stats['max'][0]/stats['rms'][0]:.1f}", flush=True)

    # Export continuum to FITS for inspection
    exportfits(
        imagename=imagename + '.image',
        fitsimage=os.path.join(output_dir, 'continuum.fits'),
        overwrite=True,
    )

    return imagename


def make_continuum_mask(continuum_imagename, output_dir, sigma_threshold=5.0):
    """Create a clean mask from the continuum image by thresholding at N*sigma."""
    cont_image = continuum_imagename + '.image'
    mask_image = os.path.join(output_dir, 'continuum_mask.image')

    if os.path.isdir(mask_image):
        shutil.rmtree(mask_image)

    # Get noise level from continuum
    stats = imstat(cont_image)
    rms = stats['rms'][0]
    threshold = sigma_threshold * rms

    print(f"  Creating mask at {sigma_threshold}*sigma = {threshold*1e3:.3f} mJy/beam", flush=True)

    # Create mask: pixels above threshold = 1, below = 0
    immath(
        imagename=cont_image,
        outfile=mask_image,
        expr=f'iif(IM0>{threshold}, 1.0, 0.0)',
    )

    # Count mask pixels
    mask_stats = imstat(mask_image)
    total_pixels = IMSIZE[0] * IMSIZE[1]
    mask_pixels = int(mask_stats['sum'][0])
    print(f"  Mask pixels: {mask_pixels}/{total_pixels} "
          f"({100*mask_pixels/total_pixels:.1f}%)", flush=True)

    # Export mask to FITS
    exportfits(
        imagename=mask_image,
        fitsimage=os.path.join(output_dir, 'continuum_mask.fits'),
        overwrite=True,
    )

    return mask_image


def clean_single_channel(input_ms, output_dir, chan, freq_mhz, mask_image,
                         nchan_total, method_name, niter=10000, threshold_mjy=None):
    """Clean a single channel using a given mask and threshold."""
    imagename = os.path.join(output_dir, f'{method_name}_chan_{chan:04d}')

    # Clean up prior products
    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = imagename + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)

    tclean_args = dict(
        vis=input_ms,
        imagename=imagename,
        specmode='cube',
        nchan=1,
        start=chan,
        width=1,
        imsize=IMSIZE,
        cell=CELL,
        niter=niter,
        weighting=WEIGHTING,
        robust=ROBUST,
        savemodel=SAVEMODEL,
        pbcor=False,
    )

    if mask_image is not None:
        tclean_args['usemask'] = 'user'
        tclean_args['mask'] = mask_image

    if threshold_mjy is not None:
        tclean_args['threshold'] = f'{threshold_mjy}mJy'
    else:
        tclean_args['threshold'] = '0.0mJy'

    t0 = time.time()
    tclean(**tclean_args)
    dt = time.time() - t0

    # Get stats
    stats = imstat(imagename + '.image')
    model_stats = imstat(imagename + '.model')
    resid_stats = imstat(imagename + '.residual')

    peak = stats['max'][0]
    rms = resid_stats['rms'][0]
    model_sum = model_stats['sum'][0]
    n_components = np.count_nonzero(model_stats['max'][0])

    print(f"    {method_name} chan {chan} ({freq_mhz:.1f} MHz): "
          f"peak={peak*1e3:.3f} mJy, residRMS={rms*1e3:.3f} mJy, "
          f"modelFlux={model_sum*1e3:.3f} mJy, time={dt:.1f}s", flush=True)

    # Export cleaned image, residual, and model to FITS
    for suffix, label in [('.image', 'image'), ('.residual', 'residual'), ('.model', 'model')]:
        exportfits(
            imagename=imagename + suffix,
            fitsimage=os.path.join(output_dir, f'{method_name}_chan_{chan:04d}_{label}.fits'),
            overwrite=True,
        )

    # Cleanup CASA products
    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = imagename + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)

    return {
        'method': method_name,
        'channel': chan,
        'freq_mhz': freq_mhz,
        'peak_mJy': peak * 1e3,
        'resid_rms_mJy': rms * 1e3,
        'model_flux_mJy': model_sum * 1e3,
        'time_s': dt,
    }


def make_dirty_channel(input_ms, output_dir, chan, freq_mhz):
    """Make a dirty (niter=0) image for a single channel for comparison."""
    imagename = os.path.join(output_dir, f'dirty_chan_{chan:04d}')

    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = imagename + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)

    tclean(
        vis=input_ms,
        imagename=imagename,
        specmode='cube',
        nchan=1,
        start=chan,
        width=1,
        imsize=IMSIZE,
        cell=CELL,
        niter=0,
        weighting=WEIGHTING,
        robust=ROBUST,
        savemodel=SAVEMODEL,
        pbcor=False,
    )

    # Export
    exportfits(
        imagename=imagename + '.image',
        fitsimage=os.path.join(output_dir, f'dirty_chan_{chan:04d}_image.fits'),
        overwrite=True,
    )

    stats = imstat(imagename + '.image')

    # Cleanup
    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = imagename + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)

    return stats


def main():
    subband_idx = int(sys.argv[1]) if len(sys.argv) > 1 else 30

    input_ms = os.path.join(SUBBANDS_DIR, f'subband_{subband_idx:03d}.ms')
    output_dir = os.path.join(TEST_DIR, f'subband_{subband_idx:03d}')

    if not os.path.exists(input_ms):
        print(f"ERROR: {input_ms} does not exist", flush=True)
        return

    os.makedirs(output_dir, exist_ok=True)

    print(f"=== Test Cleaning Strategy: subband {subband_idx} ===", flush=True)
    print(f"  Input: {input_ms}", flush=True)
    print(f"  Output: {output_dir}", flush=True)
    print(f"  Start: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)

    # Get channel info
    msmd = casatools.msmetadata()
    msmd.open(input_ms)
    nchan = msmd.nchan(0)
    freqs = msmd.chanfreqs(0)
    msmd.close()
    print(f"  Channels: {nchan}, Freq: {freqs[0]/1e6:.3f}-{freqs[-1]/1e6:.3f} MHz", flush=True)

    # Step 1: Continuum image
    print("\n--- Step 1: Continuum Image ---", flush=True)
    cont_imagename = make_continuum_image(input_ms, output_dir, nchan)

    # Step 2: Create mask from continuum
    print("\n--- Step 2: Continuum Mask ---", flush=True)
    mask_5sig = make_continuum_mask(cont_imagename, output_dir, sigma_threshold=5.0)

    # Also try 3-sigma mask
    mask_3sig_path = os.path.join(output_dir, 'continuum_mask_3sig.image')
    if os.path.isdir(mask_3sig_path):
        shutil.rmtree(mask_3sig_path)
    cont_stats = imstat(cont_imagename + '.image')
    rms = cont_stats['rms'][0]
    immath(
        imagename=cont_imagename + '.image',
        outfile=mask_3sig_path,
        expr=f'iif(IM0>{3.0*rms}, 1.0, 0.0)',
    )
    mask_3sig_stats = imstat(mask_3sig_path)
    print(f"  3-sigma mask pixels: {int(mask_3sig_stats['sum'][0])}/{IMSIZE[0]*IMSIZE[1]}", flush=True)
    exportfits(
        imagename=mask_3sig_path,
        fitsimage=os.path.join(output_dir, 'continuum_mask_3sig.fits'),
        overwrite=True,
    )

    # Step 3: Get per-channel noise estimate from a dirty image
    print("\n--- Step 3: Dirty images for comparison ---", flush=True)
    dirty_stats = {}
    for chan in TEST_CHANNELS:
        if chan >= nchan:
            continue
        freq_mhz = freqs[chan] / 1e6
        print(f"  Dirty channel {chan} ({freq_mhz:.1f} MHz)...", flush=True)
        dirty_stats[chan] = make_dirty_channel(input_ms, output_dir, chan, freq_mhz)

    # Compute typical per-channel noise from dirty images
    dirty_rms_values = [dirty_stats[ch]['rms'][0] for ch in dirty_stats]
    median_rms = np.median(dirty_rms_values) * 1e3  # mJy
    print(f"\n  Median per-channel dirty RMS: {median_rms:.3f} mJy/beam", flush=True)
    print(f"  3-sigma threshold: {3*median_rms:.3f} mJy", flush=True)

    # Step 4: Test cleaning strategies
    results = []

    print("\n--- Step 4a: Clean with 5-sigma continuum mask + 3*rms threshold ---", flush=True)
    for chan in TEST_CHANNELS:
        if chan >= nchan:
            continue
        freq_mhz = freqs[chan] / 1e6
        r = clean_single_channel(
            input_ms, output_dir, chan, freq_mhz, mask_5sig,
            nchan, method_name='mask5sig_thr3rms',
            niter=10000, threshold_mjy=3*median_rms
        )
        results.append(r)

    print("\n--- Step 4b: Clean with 3-sigma continuum mask + 3*rms threshold ---", flush=True)
    for chan in TEST_CHANNELS:
        if chan >= nchan:
            continue
        freq_mhz = freqs[chan] / 1e6
        r = clean_single_channel(
            input_ms, output_dir, chan, freq_mhz, mask_3sig_path,
            nchan, method_name='mask3sig_thr3rms',
            niter=10000, threshold_mjy=3*median_rms
        )
        results.append(r)

    print("\n--- Step 4c: Clean with NO mask, just threshold (3*rms) ---", flush=True)
    for chan in TEST_CHANNELS:
        if chan >= nchan:
            continue
        freq_mhz = freqs[chan] / 1e6
        r = clean_single_channel(
            input_ms, output_dir, chan, freq_mhz, None,
            nchan, method_name='nomask_thr3rms',
            niter=10000, threshold_mjy=3*median_rms
        )
        results.append(r)

    # Step 5: Summary
    print("\n\n=== RESULTS SUMMARY ===", flush=True)
    print(f"{'Method':<25} {'Chan':>5} {'Freq MHz':>10} {'Peak mJy':>10} "
          f"{'ResidRMS':>10} {'ModelFlux':>12} {'Time s':>8}", flush=True)
    print("-" * 85, flush=True)
    for r in results:
        print(f"{r['method']:<25} {r['channel']:>5} {r['freq_mhz']:>10.1f} "
              f"{r['peak_mJy']:>10.3f} {r['resid_rms_mJy']:>10.3f} "
              f"{r['model_flux_mJy']:>12.3f} {r['time_s']:>8.1f}", flush=True)

    # Save results to CSV
    import csv
    csv_path = os.path.join(output_dir, 'cleaning_test_results.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print(f"\nResults saved to {csv_path}", flush=True)
    print(f"FITS files saved in {output_dir}/", flush=True)
    print(f"\nEnd: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)


if __name__ == '__main__':
    main()
