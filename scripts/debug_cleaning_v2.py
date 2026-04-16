#!/usr/bin/env python3
"""
Debug cleaning v2: test multiscale clean with proper masking.

Previous attempts used default Hogbom (point-source) deconvolver with no mask,
which is wrong for extended emission like the Galactic Center.

This script tests:
1. Hogbom + no mask (baseline, what we've been doing)
2. Multiscale + no mask
3. Hogbom + auto-multithresh mask
4. Multiscale + auto-multithresh mask
5. Multiscale + continuum-derived mask

For each, reports: peak, RMS (corner + full), dynamic range, model flux, time.

Usage:
    python debug_cleaning_v2.py [subband_idx] [channel]
    Default: subband 30, channel 100
"""

import os
import sys
import time
import shutil
import numpy as np
import casatools
from casatasks import tclean, exportfits, imsubimage, imstat, immath

PROJECT_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project'
)
SUBBANDS_DIR = os.path.join(PROJECT_DIR, 'subbands')
DEBUG_DIR = os.path.join(PROJECT_DIR, 'test_cleaning', 'debug_v2')

IMSIZE = [512, 512]
CELL = '2arcsec'
WEIGHTING = 'briggs'
ROBUST = 0.5
SAVEMODEL = 'none'

# Multiscale scales: 0 = point source, then ~2x, 5x, 15x beam
# MeerKAT L-band beam ~7 arcsec = 3.5 pixels at 2"/pix
MULTISCALE_SCALES = [0, 7, 21, 63]  # pixels


def measure_noise(imagename):
    """Measure RMS in four 64x64 corners and full image."""
    nx, ny = IMSIZE
    box_size = 64
    corners = [
        '0,0,%d,%d' % (box_size-1, box_size-1),
        '%d,0,%d,%d' % (nx-box_size, nx-1, box_size-1),
        '0,%d,%d,%d' % (ny-box_size, box_size-1, ny-1),
        '%d,%d,%d,%d' % (nx-box_size, ny-box_size, nx-1, ny-1),
    ]
    rms_vals = []
    for box in corners:
        s = imstat(imagename, box=box)
        rms_vals.append(s['rms'][0])
    corner_rms = np.median(rms_vals)
    full_stats = imstat(imagename)
    return corner_rms, full_stats['rms'][0], full_stats['max'][0], full_stats['min'][0]


def run_test(input_ms, output_dir, chan, label, niter=10000, threshold_mjy=None,
             deconvolver='hogbom', usemask='', mask='',
             multiscale_scales=None,
             sidelobethreshold=2.0, noisethreshold=4.25,
             lownoisethreshold=1.5, minbeamfrac=0.3):
    """Run a single tclean test and report diagnostics."""

    imagename = os.path.join(output_dir, label)

    # Clean up
    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = imagename + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)

    print(f"\n{'='*70}", flush=True)
    print(f"  [{label}]", flush=True)
    print(f"  deconvolver={deconvolver}, usemask={usemask or 'none'}", flush=True)
    print(f"  niter={niter}, threshold={threshold_mjy} mJy", flush=True)
    if multiscale_scales:
        print(f"  scales={multiscale_scales} pixels", flush=True)

    kwargs = dict(
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
        deconvolver=deconvolver,
    )

    if threshold_mjy is not None:
        kwargs['threshold'] = '%fmJy' % threshold_mjy

    if deconvolver == 'multiscale' and multiscale_scales:
        kwargs['scales'] = multiscale_scales
        kwargs['smallscalebias'] = 0.6

    if usemask == 'auto-multithresh':
        kwargs['usemask'] = 'auto-multithresh'
        kwargs['sidelobethreshold'] = sidelobethreshold
        kwargs['noisethreshold'] = noisethreshold
        kwargs['lownoisethreshold'] = lownoisethreshold
        kwargs['minbeamfrac'] = minbeamfrac
    elif usemask == 'user' and mask:
        kwargs['usemask'] = 'user'
        kwargs['mask'] = mask

    t0 = time.time()
    try:
        tclean(**kwargs)
    except Exception as e:
        dt = time.time() - t0
        print(f"\n  *** tclean CRASHED: {e} ***", flush=True)
        print(f"  Time before crash: {dt:.1f}s", flush=True)
        # Cleanup
        for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
            p = imagename + suffix
            if os.path.isdir(p):
                shutil.rmtree(p)
        return {
            'label': label, 'img_peak_mJy': 0, 'img_corner_rms_mJy': 0,
            'img_full_rms_mJy': 0, 'res_corner_rms_mJy': 0, 'res_full_rms_mJy': 0,
            'model_flux_mJy': 0, 'n_model_pix': 0, 'DR_image': 0, 'DR_residual': 0,
            'time_s': dt,
        }
    dt = time.time() - t0

    # Diagnostics
    img = imagename + '.image'
    resid = imagename + '.residual'
    model = imagename + '.model'

    img_corner_rms, img_full_rms, img_peak, img_min = measure_noise(img)
    res_corner_rms, res_full_rms, res_peak, res_min = measure_noise(resid)

    model_stats = imstat(model)
    model_flux = model_stats['sum'][0]
    model_max = model_stats['max'][0]

    # Count model pixels
    ia = casatools.image()
    ia.open(model)
    model_data = ia.getchunk()
    ia.close()
    n_model_pix = np.count_nonzero(model_data)

    # Dynamic range
    dr_dirty_style = img_peak / img_corner_rms if img_corner_rms > 0 else 0
    dr_residual = res_peak / res_corner_rms if res_corner_rms > 0 else 0

    print(f"\n  Image:    peak={img_peak*1e3:.2f} mJy, corner_rms={img_corner_rms*1e3:.2f} mJy, "
          f"full_rms={img_full_rms*1e3:.2f} mJy", flush=True)
    print(f"  Residual: peak={res_peak*1e3:.2f} mJy, corner_rms={res_corner_rms*1e3:.2f} mJy, "
          f"full_rms={res_full_rms*1e3:.2f} mJy", flush=True)
    print(f"  Model:    flux={model_flux*1e3:.2f} mJy, max={model_max*1e3:.2f} mJy, "
          f"pixels={n_model_pix}/{model_data.size}", flush=True)
    print(f"  DR(image): {dr_dirty_style:.1f}   DR(residual): {dr_residual:.1f}", flush=True)
    print(f"  Time: {dt:.1f}s ({dt/60:.1f} min)", flush=True)

    # Export FITS
    for suffix, tag in [('.image', 'image'), ('.residual', 'residual'), ('.model', 'model')]:
        exportfits(imagename=imagename + suffix,
                   fitsimage=os.path.join(output_dir, '%s_%s.fits' % (label, tag)),
                   overwrite=True)

    # Check if mask was created
    mask_file = imagename + '.mask'
    if os.path.isdir(mask_file):
        mask_stats = imstat(mask_file)
        mask_pix = int(mask_stats['sum'][0])
        print(f"  Mask: {mask_pix}/{IMSIZE[0]*IMSIZE[1]} pixels "
              f"({100*mask_pix/(IMSIZE[0]*IMSIZE[1]):.1f}%)", flush=True)
        exportfits(imagename=mask_file,
                   fitsimage=os.path.join(output_dir, '%s_mask.fits' % label),
                   overwrite=True)

    # Cleanup CASA products
    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = imagename + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)

    return {
        'label': label,
        'img_peak_mJy': img_peak * 1e3,
        'img_corner_rms_mJy': img_corner_rms * 1e3,
        'img_full_rms_mJy': img_full_rms * 1e3,
        'res_corner_rms_mJy': res_corner_rms * 1e3,
        'res_full_rms_mJy': res_full_rms * 1e3,
        'model_flux_mJy': model_flux * 1e3,
        'n_model_pix': n_model_pix,
        'DR_image': dr_dirty_style,
        'DR_residual': dr_residual,
        'time_s': dt,
    }


def main():
    subband_idx = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    chan = int(sys.argv[2]) if len(sys.argv) > 2 else 100

    input_ms = os.path.join(SUBBANDS_DIR, 'subband_%03d.ms' % subband_idx)
    output_dir = os.path.join(DEBUG_DIR, 'sb%03d_ch%04d' % (subband_idx, chan))
    os.makedirs(output_dir, exist_ok=True)

    msmd = casatools.msmetadata()
    msmd.open(input_ms)
    nchan = msmd.nchan(0)
    freqs = msmd.chanfreqs(0)
    msmd.close()
    freq_mhz = freqs[chan] / 1e6

    print("=" * 70, flush=True)
    print(f"  Debug Cleaning v2: Multiscale + Masking Tests", flush=True)
    print(f"  Subband {subband_idx}, Channel {chan}/{nchan}, Freq {freq_mhz:.3f} MHz", flush=True)
    print(f"  MS: {input_ms}", flush=True)
    print(f"  Output: {output_dir}", flush=True)
    print(f"  Theoretical noise: ~0.24 mJy (MeerKAT 64 ant, 8h, 26 kHz)", flush=True)
    print(f"  Start: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    print("=" * 70, flush=True)

    results = []

    # 0. Dirty image (baseline)
    r = run_test(input_ms, output_dir, chan, 'dirty',
                 niter=0)
    results.append(r)

    # 1. Hogbom, no mask (what we've been doing)
    r = run_test(input_ms, output_dir, chan, 'hogbom_nomask',
                 niter=10000, threshold_mjy=0.5,
                 deconvolver='hogbom')
    results.append(r)

    # 2. Multiscale, no mask
    r = run_test(input_ms, output_dir, chan, 'multiscale_nomask',
                 niter=10000, threshold_mjy=0.5,
                 deconvolver='multiscale',
                 multiscale_scales=MULTISCALE_SCALES)
    results.append(r)

    # 3. Hogbom + auto-multithresh (tuned for per-channel)
    r = run_test(input_ms, output_dir, chan, 'hogbom_automask',
                 niter=10000, threshold_mjy=0.5,
                 deconvolver='hogbom',
                 usemask='auto-multithresh',
                 sidelobethreshold=1.25,
                 noisethreshold=3.0,
                 lownoisethreshold=1.0,
                 minbeamfrac=0.1)
    results.append(r)

    # 4. Multiscale + auto-multithresh (the recommended approach)
    r = run_test(input_ms, output_dir, chan, 'multiscale_automask',
                 niter=10000, threshold_mjy=0.5,
                 deconvolver='multiscale',
                 multiscale_scales=MULTISCALE_SCALES,
                 usemask='auto-multithresh',
                 sidelobethreshold=1.25,
                 noisethreshold=3.0,
                 lownoisethreshold=1.0,
                 minbeamfrac=0.1)
    results.append(r)

    # 5. Multiscale + auto-multithresh, more aggressive masking
    r = run_test(input_ms, output_dir, chan, 'multiscale_automask_aggr',
                 niter=10000, threshold_mjy=0.5,
                 deconvolver='multiscale',
                 multiscale_scales=MULTISCALE_SCALES,
                 usemask='auto-multithresh',
                 sidelobethreshold=1.0,
                 noisethreshold=2.0,
                 lownoisethreshold=0.5,
                 minbeamfrac=0.05)
    results.append(r)

    # 6. Multiscale + auto-multithresh, deep clean (niter=50000)
    r = run_test(input_ms, output_dir, chan, 'multiscale_automask_deep',
                 niter=50000, threshold_mjy=0.5,
                 deconvolver='multiscale',
                 multiscale_scales=MULTISCALE_SCALES,
                 usemask='auto-multithresh',
                 sidelobethreshold=1.25,
                 noisethreshold=3.0,
                 lownoisethreshold=1.0,
                 minbeamfrac=0.1)
    results.append(r)

    # Summary
    print("\n\n" + "=" * 100, flush=True)
    print("SUMMARY", flush=True)
    print("=" * 100, flush=True)
    print(f"Theoretical noise: ~0.24 mJy", flush=True)
    print(f"{'Label':<30} {'Peak':>8} {'C_RMS':>8} {'F_RMS':>8} {'R_cRMS':>8} "
          f"{'Model':>10} {'DR_img':>7} {'DR_res':>7} {'Time':>7}", flush=True)
    print("-" * 100, flush=True)
    for r in results:
        print(f"{r['label']:<30} {r['img_peak_mJy']:>8.2f} {r['img_corner_rms_mJy']:>8.2f} "
              f"{r['img_full_rms_mJy']:>8.2f} {r['res_corner_rms_mJy']:>8.2f} "
              f"{r['model_flux_mJy']:>10.2f} {r['DR_image']:>7.1f} {r['DR_residual']:>7.1f} "
              f"{r['time_s']:>7.1f}", flush=True)

    print(f"\nEnd: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)


if __name__ == '__main__':
    main()
