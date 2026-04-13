#!/usr/bin/env python3
"""
Recovery script: export FITS from existing tclean .image cubes.
Handles subbands where tclean succeeded but exportfits failed due to
the 'channel' keyword bug.

Usage:
    Called by phase3_export_recovery_submit.sh via SLURM array job.
    SLURM_ARRAY_TASK_ID selects which subband to process.
"""

import os
import sys
import time
import shutil
import casatools
from casatasks import exportfits, imsubimage

IMAGES_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/images'
)
SUBBANDS_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/subbands'
)


def main():
    subband_idx = int(os.environ.get('SLURM_ARRAY_TASK_ID',
                                      sys.argv[1] if len(sys.argv) > 1 else '0'))

    output_dir = os.path.join(IMAGES_DIR, f'subband_{subband_idx:03d}')
    cube_image = os.path.join(output_dir, f'subband_{subband_idx:03d}_cube.image')

    if not os.path.isdir(cube_image):
        print(f"Subband {subband_idx}: no .image cube at {cube_image}, skipping.", flush=True)
        return

    # Get frequencies from the subband MS
    input_ms = os.path.join(SUBBANDS_DIR, f'subband_{subband_idx:03d}.ms')
    msmd = casatools.msmetadata()
    msmd.open(input_ms)
    nchan = msmd.nchan(0)
    freqs = msmd.chanfreqs(0)
    msmd.close()

    existing_fits = [f for f in os.listdir(output_dir) if f.endswith('.fits')]
    if len(existing_fits) >= nchan:
        print(f"Subband {subband_idx}: already has {len(existing_fits)} FITS, skipping.", flush=True)
        return

    print(f"=== Export recovery: subband {subband_idx} ===", flush=True)
    print(f"  Cube: {cube_image}", flush=True)
    print(f"  Channels: {nchan}", flush=True)
    print(f"  Existing FITS: {len(existing_fits)}", flush=True)

    t0 = time.time()
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

    t_export = time.time() - t0
    print(f"  Export done: {n_exported} new files in {t_export:.1f}s ({t_export/60:.1f} min)",
          flush=True)

    # Verify and clean up CASA products
    n_fits = len([f for f in os.listdir(output_dir) if f.endswith('.fits')])
    if n_fits >= nchan:
        print(f"  All {n_fits} FITS present, cleaning up CASA products", flush=True)
        imagename_base = os.path.join(output_dir, f'subband_{subband_idx:03d}_cube')
        for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
            p = imagename_base + suffix
            if os.path.isdir(p):
                shutil.rmtree(p)
    else:
        print(f"  WARNING: Only {n_fits}/{nchan} FITS — keeping CASA products", flush=True)


if __name__ == '__main__':
    main()
