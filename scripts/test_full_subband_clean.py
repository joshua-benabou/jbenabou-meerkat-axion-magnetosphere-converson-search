#!/usr/bin/env python3
"""
Test cleaning on a full 383-channel subband with niter=100000.
This gives us the real wall-clock time for production planning.
Also compares dirty vs cleaned corner RMS across channels.
"""

import os
import sys
import time
import shutil
import glob
import numpy as np
import casatools
from casatasks import tclean, exportfits, imsubimage, imstat

SUBBANDS_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/subbands'
)
IMAGES_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/images'
)
TEST_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/test_cleaning/full_subband'
)

IMSIZE = [512, 512]
CELL = '2arcsec'
WEIGHTING = 'briggs'
ROBUST = 0.5
SAVEMODEL = 'none'


def measure_peak_from_dirty(subband_idx):
    dirty_dir = os.path.join(IMAGES_DIR, f'subband_{subband_idx:03d}')
    fits_files = sorted(glob.glob(os.path.join(dirty_dir, 'chan_*.fits')))
    if not fits_files:
        return None
    n = len(fits_files)
    indices = np.linspace(max(1, n // 20), n - 1, min(10, n)).astype(int)
    from astropy.io import fits as pyfits
    peaks = []
    for idx in indices:
        with pyfits.open(fits_files[idx]) as hdul:
            data = hdul[0].data
            while data.ndim > 2:
                data = data[0]
            peaks.append(np.nanmax(data))
    return np.median(peaks)


def main():
    subband_idx = int(sys.argv[1]) if len(sys.argv) > 1 else 30

    input_ms = os.path.join(SUBBANDS_DIR, f'subband_{subband_idx:03d}.ms')
    output_dir = os.path.join(TEST_DIR, f'subband_{subband_idx:03d}')
    os.makedirs(output_dir, exist_ok=True)

    msmd = casatools.msmetadata()
    msmd.open(input_ms)
    nchan = msmd.nchan(0)
    freqs = msmd.chanfreqs(0)
    msmd.close()

    peak_jy = measure_peak_from_dirty(subband_idx)
    threshold_mjy = 0.10 * peak_jy * 1e3
    niter = 100000

    print(f"=== Full subband cleaning test ===", flush=True)
    print(f"  Subband: {subband_idx}, Channels: {nchan}", flush=True)
    print(f"  Freq: {freqs[0]/1e6:.3f}-{freqs[-1]/1e6:.3f} MHz", flush=True)
    print(f"  Peak from dirty: {peak_jy*1e3:.3f} mJy", flush=True)
    print(f"  Threshold: {threshold_mjy:.3f} mJy (10% of peak)", flush=True)
    print(f"  niter: {niter}", flush=True)
    print(f"  Start: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)

    imagename = os.path.join(output_dir, f'subband_{subband_idx:03d}_cube')

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
        niter=niter,
        weighting=WEIGHTING,
        robust=ROBUST,
        savemodel=SAVEMODEL,
        pbcor=False,
        threshold=f'{threshold_mjy}mJy',
    )
    t_tclean = time.time() - t0
    print(f"  tclean done: {t_tclean:.1f}s ({t_tclean/60:.1f} min)", flush=True)

    # Check model stats
    model_stats = imstat(imagename + '.model')
    print(f"  Model flux: {model_stats['sum'][0]*1e3:.1f} mJy", flush=True)
    print(f"  Model max:  {model_stats['max'][0]*1e3:.3f} mJy", flush=True)

    # Export a few test channels and compare to dirty
    test_chans = [10, 50, 100, 150, 200, 300]
    print(f"\n  === Dirty vs Cleaned comparison ===", flush=True)
    print(f"  {'Chan':>6} {'D_crms':>8} {'C_crms':>8} {'improv%':>8}", flush=True)

    from astropy.io import fits as pyfits

    for ch in test_chans:
        if ch >= nchan:
            continue
        freq_mhz = freqs[ch] / 1e6

        # Export cleaned channel
        chan_image = os.path.join(output_dir, f'_tmp_chan_{ch:04d}.image')
        if os.path.isdir(chan_image):
            shutil.rmtree(chan_image)
        imsubimage(imagename=imagename + '.image', outfile=chan_image, chans=str(ch))
        clean_fits = os.path.join(output_dir, f'clean_chan_{ch:04d}.fits')
        exportfits(imagename=chan_image, fitsimage=clean_fits, overwrite=True)
        shutil.rmtree(chan_image)

        # Read cleaned
        with pyfits.open(clean_fits) as h:
            cd = h[0].data
            while cd.ndim > 2:
                cd = cd[0]

        # Read dirty
        dirty_files = sorted(glob.glob(os.path.join(
            IMAGES_DIR, f'subband_{subband_idx:03d}', f'chan_{ch:04d}_*.fits')))
        with pyfits.open(dirty_files[0]) as h:
            dd = h[0].data
            while dd.ndim > 2:
                dd = dd[0]

        # Corner RMS (all 4 corners, 64x64 boxes)
        def corner_rms(d):
            corners = [d[:64,:64], d[:64,-64:], d[-64:,:64], d[-64:,-64:]]
            return np.median([np.sqrt(np.mean(c**2)) for c in corners]) * 1e3

        d_crms = corner_rms(dd)
        c_crms = corner_rms(cd)
        improv = (1 - c_crms / d_crms) * 100
        print(f"  {ch:>6} {d_crms:>8.2f} {c_crms:>8.2f} {improv:>7.1f}%", flush=True)

    # Cleanup CASA products (keep FITS)
    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = imagename + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)

    t_total = time.time() - t0
    print(f"\n  Total time: {t_total:.1f}s ({t_total/60:.1f} min)", flush=True)
    print(f"  End: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)

    # Budget estimate
    hours_per_subband = t_total / 3600
    total_hours = hours_per_subband * 86
    print(f"\n  === Budget estimate ===", flush=True)
    print(f"  Wall time per subband: {hours_per_subband:.2f} hours", flush=True)
    print(f"  Total for 86 subbands: {total_hours:.1f} node-hours", flush=True)
    print(f"  Core-hours (56 cores): {total_hours * 56:.0f}", flush=True)


if __name__ == '__main__':
    main()
