#!/usr/bin/env python3
"""
Validate Phase 3 FITS output: check completeness, file integrity, WCS headers,
and frequency monotonicity across all subbands.

Usage:
    python3 validate_fits.py [--images-dir DIR] [--subbands-dir DIR] [--fix-list FILE]

Produces a text report to stdout. Optionally writes a file listing subbands that
need re-imaging (for easy re-submission).
"""

import argparse
import os
import sys
import glob
import re
from collections import defaultdict

import numpy as np

try:
    from astropy.io import fits as pyfits
    HAS_ASTROPY = True
except ImportError:
    HAS_ASTROPY = False


BASE = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project'
)
DEFAULT_IMAGES = os.path.join(BASE, 'images')
DEFAULT_SUBBANDS = os.path.join(BASE, 'subbands')

EXPECTED_NCHAN = 383
EXPECTED_IMSIZE = 512
EXPECTED_CELL_DEG = 2.0 / 3600  # 2 arcsec in degrees

# Regex for FITS filenames produced by phase3_image_cube.py
FITS_RE = re.compile(r'^chan_(\d{4})_([0-9.]+)MHz\.fits$')


def discover_subbands(subbands_dir):
    """Return sorted list of subband indices that have a .ms on disk."""
    indices = []
    for name in os.listdir(subbands_dir):
        m = re.match(r'^subband_(\d{3})\.ms$', name)
        if m:
            indices.append(int(m.group(1)))
    return sorted(indices)


def validate_subband(images_dir, idx, check_headers=True):
    """Validate one subband's FITS output. Returns a dict of findings."""
    result = {
        'idx': idx,
        'status': 'ok',
        'n_fits': 0,
        'n_expected': EXPECTED_NCHAN,
        'missing_chans': [],
        'zero_size': [],
        'bad_header': [],
        'freqs_mhz': [],
        'warnings': [],
    }

    subdir = os.path.join(images_dir, f'subband_{idx:03d}')
    if not os.path.isdir(subdir):
        result['status'] = 'missing_dir'
        return result

    # Gather FITS files
    fits_files = {}
    for fname in os.listdir(subdir):
        m = FITS_RE.match(fname)
        if m:
            chan = int(m.group(1))
            freq = float(m.group(2))
            fits_files[chan] = (fname, freq)

    result['n_fits'] = len(fits_files)

    if result['n_fits'] == 0:
        # Check for lingering CASA products (tclean still running or failed)
        casa_products = [f for f in os.listdir(subdir)
                         if any(f.endswith(s) for s in
                                ('.image', '.psf', '.pb', '.residual', '.model', '.sumwt'))]
        if casa_products:
            result['status'] = 'in_progress'
            result['warnings'].append(
                f'No FITS but CASA products present: {", ".join(casa_products)}')
        else:
            result['status'] = 'empty'
        return result

    # Check channel completeness
    expected_chans = set(range(EXPECTED_NCHAN))
    present_chans = set(fits_files.keys())
    result['missing_chans'] = sorted(expected_chans - present_chans)

    # Check for zero-size files
    for chan, (fname, freq) in fits_files.items():
        fpath = os.path.join(subdir, fname)
        if os.path.getsize(fpath) == 0:
            result['zero_size'].append(chan)

    # Collect frequencies from filenames
    freq_list = []
    for chan in sorted(fits_files.keys()):
        freq_list.append((chan, fits_files[chan][1]))
    result['freqs_mhz'] = freq_list

    # Check frequency monotonicity from filenames
    if len(freq_list) > 1:
        freqs = [f for _, f in freq_list]
        diffs = np.diff(freqs)
        if not np.all(diffs > 0):
            result['warnings'].append('Filename frequencies are not strictly increasing')

    # Spot-check headers on first, middle, last channel
    if check_headers and HAS_ASTROPY and len(fits_files) > 0:
        chans_to_check = sorted(fits_files.keys())
        spot = [chans_to_check[0], chans_to_check[len(chans_to_check) // 2],
                chans_to_check[-1]]
        for chan in spot:
            fname, freq = fits_files[chan]
            fpath = os.path.join(subdir, fname)
            try:
                with pyfits.open(fpath) as hdul:
                    hdr = hdul[0].header
                    # Check image dimensions
                    naxis1 = hdr.get('NAXIS1', 0)
                    naxis2 = hdr.get('NAXIS2', 0)
                    if naxis1 != EXPECTED_IMSIZE or naxis2 != EXPECTED_IMSIZE:
                        result['bad_header'].append(
                            f'chan {chan}: NAXIS1={naxis1}, NAXIS2={naxis2} '
                            f'(expected {EXPECTED_IMSIZE})')
                    # Check pixel scale
                    cdelt1 = abs(hdr.get('CDELT1', 0))
                    cdelt2 = abs(hdr.get('CDELT2', 0))
                    if cdelt1 > 0 and abs(cdelt1 - EXPECTED_CELL_DEG) / EXPECTED_CELL_DEG > 0.01:
                        result['bad_header'].append(
                            f'chan {chan}: CDELT1={cdelt1*3600:.2f}arcsec '
                            f'(expected {EXPECTED_CELL_DEG*3600:.1f})')
                    # Check data isn't all NaN/zero
                    data = hdul[0].data
                    if data is not None:
                        finite = data[np.isfinite(data)]
                        if finite.size == 0:
                            result['bad_header'].append(f'chan {chan}: all NaN/Inf')
                        elif np.all(finite == 0):
                            result['bad_header'].append(f'chan {chan}: all zeros')
            except Exception as e:
                result['bad_header'].append(f'chan {chan}: read error: {e}')

    # Determine overall status
    if result['missing_chans'] or result['zero_size'] or result['bad_header']:
        result['status'] = 'incomplete' if result['missing_chans'] else 'issues'
    elif result['n_fits'] < EXPECTED_NCHAN:
        result['status'] = 'incomplete'

    return result


def main():
    parser = argparse.ArgumentParser(description='Validate Phase 3 FITS output')
    parser.add_argument('--images-dir', default=DEFAULT_IMAGES)
    parser.add_argument('--subbands-dir', default=DEFAULT_SUBBANDS)
    parser.add_argument('--fix-list', default=None,
                        help='Write file listing subband indices that need re-imaging')
    parser.add_argument('--no-headers', action='store_true',
                        help='Skip FITS header checks (faster)')
    args = parser.parse_args()

    subband_indices = discover_subbands(args.subbands_dir)
    n_subbands = len(subband_indices)
    print(f"=== FITS Validation Report ===")
    print(f"Subbands on disk: {n_subbands}")
    print(f"Expected channels per subband: {EXPECTED_NCHAN}")
    print(f"Expected total FITS files: {n_subbands * EXPECTED_NCHAN}")
    print()

    results = []
    total_fits = 0
    status_counts = defaultdict(int)
    needs_rerun = []

    for idx in subband_indices:
        r = validate_subband(args.images_dir, idx, check_headers=not args.no_headers)
        results.append(r)
        total_fits += r['n_fits']
        status_counts[r['status']] += 1

        if r['status'] != 'ok':
            needs_rerun.append(idx)

    # Summary
    print(f"--- Summary ---")
    print(f"Total FITS files found: {total_fits} / {n_subbands * EXPECTED_NCHAN}")
    print(f"Subbands complete (ok):   {status_counts['ok']}")
    print(f"Subbands in progress:     {status_counts['in_progress']}")
    print(f"Subbands incomplete:      {status_counts['incomplete']}")
    print(f"Subbands with issues:     {status_counts['issues']}")
    print(f"Subbands missing dir:     {status_counts['missing_dir']}")
    print(f"Subbands empty:           {status_counts['empty']}")
    print()

    # Frequency coverage
    all_freqs = []
    for r in results:
        all_freqs.extend([f for _, f in r['freqs_mhz']])
    if all_freqs:
        all_freqs.sort()
        print(f"--- Frequency Coverage ---")
        print(f"Min: {all_freqs[0]:.3f} MHz")
        print(f"Max: {all_freqs[-1]:.3f} MHz")
        print(f"Total unique channels with FITS: {len(all_freqs)}")
        # Check for gaps between subbands
        subband_ranges = []
        for r in results:
            if r['freqs_mhz']:
                freqs = [f for _, f in r['freqs_mhz']]
                subband_ranges.append((r['idx'], min(freqs), max(freqs)))
        subband_ranges.sort(key=lambda x: x[1])
        gaps = []
        for i in range(len(subband_ranges) - 1):
            gap = subband_ranges[i + 1][1] - subband_ranges[i][2]
            if gap > 0.1:  # >100 kHz gap
                gaps.append((subband_ranges[i][0], subband_ranges[i + 1][0], gap))
        if gaps:
            print(f"Frequency gaps detected ({len(gaps)}):")
            for sb_a, sb_b, gap_mhz in gaps:
                print(f"  Between subband {sb_a:03d} and {sb_b:03d}: {gap_mhz:.3f} MHz")
        else:
            print("No frequency gaps detected")
        print()

    # Per-subband detail for problematic ones
    problem_results = [r for r in results if r['status'] != 'ok']
    if problem_results:
        print(f"--- Problem Subbands ({len(problem_results)}) ---")
        for r in problem_results:
            print(f"\nSubband {r['idx']:03d}: {r['status']} "
                  f"({r['n_fits']}/{r['n_expected']} FITS)")
            if r['missing_chans']:
                if len(r['missing_chans']) <= 10:
                    print(f"  Missing channels: {r['missing_chans']}")
                else:
                    print(f"  Missing {len(r['missing_chans'])} channels: "
                          f"{r['missing_chans'][:5]}...{r['missing_chans'][-5:]}")
            if r['zero_size']:
                print(f"  Zero-size files: channels {r['zero_size']}")
            if r['bad_header']:
                for msg in r['bad_header']:
                    print(f"  Header issue: {msg}")
            if r['warnings']:
                for w in r['warnings']:
                    print(f"  Warning: {w}")
    else:
        print("All subbands passed validation!")

    # Write fix list
    if args.fix_list and needs_rerun:
        with open(args.fix_list, 'w') as f:
            for idx in needs_rerun:
                f.write(f"{idx}\n")
        print(f"\nWrote {len(needs_rerun)} subband indices to {args.fix_list}")

    # Exit code: 0 if all ok, 1 if any problems
    return 0 if not needs_rerun else 1


if __name__ == '__main__':
    sys.exit(main())
