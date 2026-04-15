#!/usr/bin/env python3
"""
Debug cleaning: test tclean on a single channel with correct noise measurement.

The previous test failed because imstat RMS over the full image includes
GC emission + sidelobes (~95 mJy), which is NOT the thermal noise.
Setting threshold = 3*95 = 285 mJy > peak (270 mJy) → zero cleaning.

Fix: measure noise in the image corners (off-source), use that for threshold.
Also test several threshold levels to find what works.

Usage:
    python debug_cleaning.py [subband_idx] [channel]
    Default: subband 30, channel 100
"""

import os
import sys
import time
import shutil
import numpy as np
import casatools
from casatasks import tclean, exportfits, imsubimage, imstat

PROJECT_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project'
)
SUBBANDS_DIR = os.path.join(PROJECT_DIR, 'subbands')
DEBUG_DIR = os.path.join(PROJECT_DIR, 'test_cleaning', 'debug')

IMSIZE = [512, 512]
CELL = '2arcsec'
WEIGHTING = 'briggs'
ROBUST = 0.5
SAVEMODEL = 'none'


def measure_off_source_noise(imagename):
    """Measure RMS in the four corners of the image (off-source regions)."""
    nx, ny = IMSIZE
    # Use 64x64 pixel boxes in each corner
    box_size = 64
    corners = [
        f'0,0,{box_size-1},{box_size-1}',                          # bottom-left
        f'{nx-box_size},0,{nx-1},{box_size-1}',                    # bottom-right
        f'0,{ny-box_size},{box_size-1},{ny-1}',                    # top-left
        f'{nx-box_size},{ny-box_size},{nx-1},{ny-1}',              # top-right
    ]

    rms_values = []
    for i, box in enumerate(corners):
        stats = imstat(imagename, box=box)
        rms_values.append(stats['rms'][0])
        print(f"    Corner {i} RMS: {stats['rms'][0]*1e3:.3f} mJy/beam "
              f"(max={stats['max'][0]*1e3:.3f}, min={stats['min'][0]*1e3:.3f})")

    median_rms = np.median(rms_values)
    print(f"    Median corner RMS: {median_rms*1e3:.3f} mJy/beam")
    return median_rms


def run_tclean(input_ms, imagename, chan, niter, threshold_mjy, label):
    """Run tclean on a single channel and report results."""
    # Clean up prior products
    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = imagename + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)

    print(f"\n  [{label}] niter={niter}, threshold={threshold_mjy:.3f} mJy")
    t0 = time.time()
    tclean(
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
        threshold=f'{threshold_mjy}mJy',
    )
    dt = time.time() - t0

    # Stats
    img = imagename + '.image'
    model = imagename + '.model'
    resid = imagename + '.residual'

    img_stats = imstat(img)
    model_stats = imstat(model)
    resid_stats = imstat(resid)

    peak = img_stats['max'][0]
    model_flux = model_stats['sum'][0]
    model_max = model_stats['max'][0]
    resid_rms_full = resid_stats['rms'][0]

    # Off-source residual noise
    print(f"    Residual off-source noise:")
    resid_rms_corner = measure_off_source_noise(resid)

    # Count non-zero model pixels
    # Read model image to count
    ia = casatools.image()
    ia.open(model)
    model_data = ia.getchunk()
    ia.close()
    n_model_pixels = np.count_nonzero(model_data)
    total_pixels = model_data.size

    print(f"    Peak:            {peak*1e3:.3f} mJy/beam")
    print(f"    Model flux:      {model_flux*1e3:.3f} mJy")
    print(f"    Model max:       {model_max*1e3:.3f} mJy/beam")
    print(f"    Model pixels:    {n_model_pixels}/{total_pixels}")
    print(f"    Resid RMS(full): {resid_rms_full*1e3:.3f} mJy/beam")
    print(f"    Resid RMS(corn): {resid_rms_corner*1e3:.3f} mJy/beam")
    print(f"    Time:            {dt:.1f}s")

    # Export cleaned image and residual
    exportfits(imagename=img,
               fitsimage=imagename + '_image.fits', overwrite=True)
    exportfits(imagename=resid,
               fitsimage=imagename + '_residual.fits', overwrite=True)
    exportfits(imagename=model,
               fitsimage=imagename + '_model.fits', overwrite=True)

    # Cleanup CASA products
    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = imagename + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)

    return {
        'label': label,
        'niter': niter,
        'threshold_mjy': threshold_mjy,
        'peak_mJy': peak * 1e3,
        'model_flux_mJy': model_flux * 1e3,
        'model_max_mJy': model_max * 1e3,
        'n_model_pixels': n_model_pixels,
        'resid_rms_full_mJy': resid_rms_full * 1e3,
        'resid_rms_corner_mJy': resid_rms_corner * 1e3,
        'time_s': dt,
    }


def main():
    subband_idx = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    chan = int(sys.argv[2]) if len(sys.argv) > 2 else 100

    input_ms = os.path.join(SUBBANDS_DIR, f'subband_{subband_idx:03d}.ms')
    output_dir = os.path.join(DEBUG_DIR, f'sb{subband_idx:03d}_ch{chan:04d}')

    if not os.path.exists(input_ms):
        print(f"ERROR: {input_ms} does not exist")
        return

    os.makedirs(output_dir, exist_ok=True)

    # Get freq info
    msmd = casatools.msmetadata()
    msmd.open(input_ms)
    nchan = msmd.nchan(0)
    freqs = msmd.chanfreqs(0)
    msmd.close()
    freq_mhz = freqs[chan] / 1e6

    print(f"=== Debug Cleaning ===")
    print(f"  Subband: {subband_idx}, Channel: {chan}/{nchan}")
    print(f"  Frequency: {freq_mhz:.3f} MHz")
    print(f"  MS: {input_ms}")
    print(f"  Output: {output_dir}")
    print(f"  Start: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    # Step 1: Make dirty image and measure CORRECT noise
    print(f"\n--- Step 1: Dirty image + noise measurement ---")
    dirty_name = os.path.join(output_dir, 'dirty')
    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = dirty_name + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)

    tclean(
        vis=input_ms, imagename=dirty_name,
        specmode='cube', nchan=1, start=chan, width=1,
        imsize=IMSIZE, cell=CELL, niter=0,
        weighting=WEIGHTING, robust=ROBUST, savemodel=SAVEMODEL, pbcor=False,
    )

    dirty_img = dirty_name + '.image'
    full_stats = imstat(dirty_img)
    print(f"  Full image: peak={full_stats['max'][0]*1e3:.3f} mJy, "
          f"RMS={full_stats['rms'][0]*1e3:.3f} mJy, "
          f"SNR(wrong)={full_stats['max'][0]/full_stats['rms'][0]:.1f}")

    print(f"  Off-source corner noise:")
    corner_rms = measure_off_source_noise(dirty_img)

    print(f"\n  *** CORRECT SNR: {full_stats['max'][0]/corner_rms:.1f} ***")
    print(f"  *** Peak / corner_RMS: {full_stats['max'][0]*1e3:.3f} / {corner_rms*1e3:.3f} mJy ***")

    exportfits(imagename=dirty_img,
               fitsimage=os.path.join(output_dir, 'dirty_image.fits'), overwrite=True)
    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = dirty_name + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)

    # Step 2: Test cleaning at several threshold levels
    print(f"\n--- Step 2: Test cleaning at various thresholds ---")
    peak_mjy = full_stats['max'][0] * 1e3
    corner_rms_mjy = corner_rms * 1e3

    # Thresholds to test (in mJy)
    tests = [
        ('3sig_corner', 3.0 * corner_rms_mjy, 10000),
        ('5sig_corner', 5.0 * corner_rms_mjy, 10000),
        ('10sig_corner', 10.0 * corner_rms_mjy, 10000),
        ('1pct_peak', 0.01 * peak_mjy, 10000),
        ('5pct_peak', 0.05 * peak_mjy, 10000),
        ('10pct_peak', 0.10 * peak_mjy, 10000),
        ('auto_multithresh', None, 10000),  # auto-multithresh with lower noisethreshold
    ]

    results = []
    for label, thr_mjy, niter in tests:
        imagename = os.path.join(output_dir, label)

        if label == 'auto_multithresh':
            # Test auto-multithresh with much lower noisethreshold
            for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
                p = imagename + suffix
                if os.path.isdir(p):
                    shutil.rmtree(p)
            print(f"\n  [auto_multithresh] niter={niter}, noisethreshold=3.0, sidelobethreshold=1.25")
            t0 = time.time()
            tclean(
                vis=input_ms, imagename=imagename,
                specmode='cube', nchan=1, start=chan, width=1,
                imsize=IMSIZE, cell=CELL, niter=niter,
                weighting=WEIGHTING, robust=ROBUST, savemodel=SAVEMODEL, pbcor=False,
                usemask='auto-multithresh',
                sidelobethreshold=1.25,
                noisethreshold=3.0,
                lownoisethreshold=1.0,
                minbeamfrac=0.1,
            )
            dt = time.time() - t0

            img_stats = imstat(imagename + '.image')
            model_stats = imstat(imagename + '.model')
            resid_stats = imstat(imagename + '.residual')

            ia = casatools.image()
            ia.open(imagename + '.model')
            model_data = ia.getchunk()
            ia.close()

            print(f"    Peak:         {img_stats['max'][0]*1e3:.3f} mJy")
            print(f"    Model flux:   {model_stats['sum'][0]*1e3:.3f} mJy")
            print(f"    Model pixels: {np.count_nonzero(model_data)}/{model_data.size}")
            print(f"    Resid RMS:    {resid_stats['rms'][0]*1e3:.3f} mJy")
            print(f"    Resid corner:")
            rc = measure_off_source_noise(imagename + '.residual')
            print(f"    Time:         {dt:.1f}s")

            exportfits(imagename=imagename + '.image',
                       fitsimage=imagename + '_image.fits', overwrite=True)
            exportfits(imagename=imagename + '.residual',
                       fitsimage=imagename + '_residual.fits', overwrite=True)

            for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
                p = imagename + suffix
                if os.path.isdir(p):
                    shutil.rmtree(p)

            results.append({
                'label': 'auto_multithresh',
                'threshold_mjy': 'auto',
                'model_flux_mJy': model_stats['sum'][0]*1e3,
                'n_model_pixels': int(np.count_nonzero(model_data)),
                'resid_rms_corner_mJy': rc*1e3,
                'time_s': dt,
            })
        else:
            r = run_tclean(input_ms, imagename, chan, niter, thr_mjy, label)
            results.append(r)

    # Summary
    print(f"\n\n{'='*80}")
    print(f"=== SUMMARY ===")
    print(f"  Peak flux: {peak_mjy:.3f} mJy")
    print(f"  Full-image RMS (WRONG): {full_stats['rms'][0]*1e3:.3f} mJy")
    print(f"  Corner RMS (CORRECT):   {corner_rms_mjy:.3f} mJy")
    print(f"  Correct SNR:            {peak_mjy/corner_rms_mjy:.1f}")
    print(f"{'='*80}")
    print(f"{'Label':<20} {'Thr mJy':>10} {'ModelFlux':>12} {'ModelPix':>10} "
          f"{'ResidRMS_c':>12} {'Time':>8}")
    print(f"{'-'*75}")
    for r in results:
        thr = f"{r['threshold_mjy']:.3f}" if isinstance(r['threshold_mjy'], float) else r['threshold_mjy']
        print(f"{r['label']:<20} {thr:>10} {r.get('model_flux_mJy',0):.3f}"
              f"{'':>3} {r.get('n_model_pixels',0):>10} "
              f"{r.get('resid_rms_corner_mJy',0):.3f}{'':>5} "
              f"{r.get('time_s',0):.1f}")

    print(f"\nEnd: {time.strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == '__main__':
    main()
