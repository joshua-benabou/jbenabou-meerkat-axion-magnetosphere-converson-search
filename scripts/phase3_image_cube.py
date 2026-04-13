#!/usr/bin/env python3
"""
Phase 3: Image all channels in a subband using tclean cube mode.

One tclean call per subband produces a 383-channel image cube.
Then exports each channel to individual FITS files.

Usage:
    Called by phase3_cube_submit.sh via SLURM array job.
    SLURM_ARRAY_TASK_ID selects which subband to process.
"""

import os
import sys
import time
import shutil
import casatools
from casatasks import tclean, exportfits, imsubimage

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
NITER = 0
WEIGHTING = 'briggs'
ROBUST = 0.5
SAVEMODEL = 'none'


def main():
    subband_idx = int(os.environ.get('SLURM_ARRAY_TASK_ID', sys.argv[1] if len(sys.argv) > 1 else '0'))

    input_ms = os.path.join(SUBBANDS_DIR, f'subband_{subband_idx:03d}.ms')
    output_dir = os.path.join(IMAGES_DIR, f'subband_{subband_idx:03d}')

    if not os.path.exists(input_ms):
        print(f"Subband {subband_idx}: {input_ms} does not exist, skipping.", flush=True)
        return

    os.makedirs(output_dir, exist_ok=True)

    # Get channel info
    msmd = casatools.msmetadata()
    msmd.open(input_ms)
    nchan = msmd.nchan(0)
    freqs = msmd.chanfreqs(0)
    msmd.close()

    # Check if already complete
    existing_fits = [f for f in os.listdir(output_dir) if f.endswith('.fits')]
    if len(existing_fits) >= nchan:
        print(f"Subband {subband_idx}: already has {len(existing_fits)} FITS files "
              f"(expected {nchan}). Skipping.", flush=True)
        return

    print(f"=== Phase 3: Imaging subband {subband_idx} ===", flush=True)
    print(f"  Input: {input_ms}", flush=True)
    print(f"  Channels: {nchan}", flush=True)
    print(f"  Freq range: {freqs[0]/1e6:.3f}-{freqs[-1]/1e6:.3f} MHz", flush=True)
    print(f"  Existing FITS: {len(existing_fits)}", flush=True)
    print(f"  Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)

    # Step 1: tclean cube mode
    imagename = os.path.join(output_dir, f'subband_{subband_idx:03d}_cube')

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
    )
    t_tclean = time.time() - t0
    print(f"  tclean cube done: {t_tclean:.1f}s ({t_tclean/60:.1f} min)", flush=True)

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
        # Extract single channel with imsubimage, then export to FITS
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
            print(f"    Exported {n_exported}/{nchan} "
                  f"(chan {chan})", flush=True)

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
        print(f"  WARNING: Only {n_fits}/{nchan} FITS files — keeping CASA products "
              f"for recovery", flush=True)

    t_total = time.time() - t0
    print(f"  Total: {t_total:.1f}s ({t_total/60:.1f} min), {n_fits} FITS files", flush=True)
    print(f"  End time: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)


if __name__ == '__main__':
    main()
