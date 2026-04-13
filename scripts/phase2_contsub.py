#!/usr/bin/env python3
"""
Phase 2: Continuum subtraction via uvcontsub on each subband.

For each subband MS, fits a polynomial to the visibilities across channels
and subtracts it, leaving only the spectral line signal. This removes the
bright Galactic Center continuum emission.

Usage:
    Called by phase2_submit.sh via SLURM array job.
    SLURM_ARRAY_TASK_ID selects which subband to process.
"""

import os
import sys
import time
from casatasks import uvcontsub

# === Configuration ===
SUBBANDS_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/subbands'
)

# uvcontsub parameters
FITORDER = 1  # Linear fit to continuum; increase if needed


def main():
    subband_idx = int(os.environ.get('SLURM_ARRAY_TASK_ID', sys.argv[1] if len(sys.argv) > 1 else '0'))

    input_ms = os.path.join(SUBBANDS_DIR, f'subband_{subband_idx:03d}.ms')
    # uvcontsub writes output to <vis>.contsub by default
    output_ms = input_ms + '.contsub'

    if not os.path.exists(input_ms):
        print(f"Subband {subband_idx}: input {input_ms} does not exist, skipping.")
        return

    if os.path.exists(output_ms):
        print(f"Subband {subband_idx}: output {output_ms} already exists, skipping.")
        return

    print(f"=== Phase 2: Continuum subtraction for subband {subband_idx} ===")
    print(f"  Input: {input_ms}")
    print(f"  Output: {output_ms}")
    print(f"  Fit order: {FITORDER}")
    print(f"  Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    t0 = time.time()

    uvcontsub(
        vis=input_ms,
        fitorder=FITORDER,
    )

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s ({elapsed/60:.1f} min)")
    print(f"  End time: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    if os.path.exists(output_ms):
        print(f"  Output verified: {output_ms} exists")
    else:
        print(f"  WARNING: Expected output {output_ms} not found!")


if __name__ == '__main__':
    main()
