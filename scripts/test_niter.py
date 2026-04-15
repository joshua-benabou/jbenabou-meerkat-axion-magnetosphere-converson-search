#!/usr/bin/env python3
"""
Quick test: cube-mode tclean with high niter to let threshold control stopping.
Tests niter=1000000 on a small subset (20 channels) of subband 30.
"""

import os
import sys
import time
import shutil
import numpy as np
import casatools
from casatasks import tclean, imstat

SUBBANDS_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/subbands'
)
TEST_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/test_cleaning/niter_test'
)

IMSIZE = [512, 512]
CELL = '2arcsec'

def run_test(input_ms, niter, nchan_test, start_chan, threshold_mjy, label):
    output_dir = os.path.join(TEST_DIR, label)
    os.makedirs(output_dir, exist_ok=True)
    imagename = os.path.join(output_dir, 'cube')

    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = imagename + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)

    print(f"\n=== {label}: niter={niter}, nchan={nchan_test}, threshold={threshold_mjy} mJy ===", flush=True)
    t0 = time.time()
    tclean(
        vis=input_ms,
        imagename=imagename,
        specmode='cube',
        nchan=nchan_test,
        start=start_chan,
        width=1,
        imsize=IMSIZE,
        cell=CELL,
        niter=niter,
        weighting='briggs',
        robust=0.5,
        savemodel='none',
        pbcor=False,
        threshold=f'{threshold_mjy}mJy',
    )
    dt = time.time() - t0
    print(f"  Time: {dt:.1f}s ({dt/60:.1f} min)", flush=True)

    # Check model
    model_stats = imstat(imagename + '.model')
    resid_stats = imstat(imagename + '.residual')
    print(f"  Model flux: {model_stats['sum'][0]*1e3:.1f} mJy", flush=True)
    print(f"  Model max:  {model_stats['max'][0]*1e3:.3f} mJy", flush=True)
    print(f"  Resid RMS:  {resid_stats['rms'][0]*1e3:.3f} mJy", flush=True)

    # Cleanup
    for suffix in ['.image', '.model', '.pb', '.psf', '.residual', '.sumwt', '.mask']:
        p = imagename + suffix
        if os.path.isdir(p):
            shutil.rmtree(p)

    return dt


def main():
    input_ms = os.path.join(SUBBANDS_DIR, 'subband_030.ms')
    os.makedirs(TEST_DIR, exist_ok=True)

    msmd = casatools.msmetadata()
    msmd.open(input_ms)
    freqs = msmd.chanfreqs(0)
    msmd.close()

    nchan_test = 20  # small subset for speed
    start_chan = 90   # near channel 100
    threshold_mjy = 27.237

    # Test different niter values
    for niter in [10000, 100000, 1000000]:
        run_test(input_ms, niter, nchan_test, start_chan, threshold_mjy,
                 f'niter_{niter}')

    print("\nDone.", flush=True)


if __name__ == '__main__':
    main()
