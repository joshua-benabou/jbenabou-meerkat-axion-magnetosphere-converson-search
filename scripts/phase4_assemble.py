#!/usr/bin/env python3
"""
Phase 4: Assemble FITS images into a master catalog and validate.

Scans all subband image directories, builds a CSV catalog mapping
filename -> frequency (MHz), and validates FITS headers. Optionally
creates a single FITS cube.

Usage:
    python3 phase4_assemble.py [--cube]
"""

import os
import sys
import glob
import csv
from astropy.io import fits
import numpy as np

IMAGES_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/images'
)
CATALOG_PATH = os.path.join(IMAGES_DIR, 'channel_catalog.csv')


def build_catalog(img_type='dirty'):
    """Scan all subband directories and build frequency-sorted catalog."""
    pattern = os.path.join(IMAGES_DIR, 'subband_*', img_type, '*.fits')
    fits_files = sorted(glob.glob(pattern))

    print(f"Found {len(fits_files)} {img_type} FITS files")

    catalog = []
    bad_files = []

    for fpath in fits_files:
        try:
            hdr = fits.getheader(fpath)
            freq_hz = hdr.get('CRVAL3') or hdr.get('CRVAL4') or hdr.get('RESTFRQ')
            if freq_hz is None:
                bad_files.append((fpath, 'no frequency in header'))
                continue

            freq_mhz = float(freq_hz) / 1e6
            bmaj = hdr.get('BMAJ', np.nan)
            bmin = hdr.get('BMIN', np.nan)
            bpa = hdr.get('BPA', np.nan)

            catalog.append({
                'filename': os.path.relpath(fpath, IMAGES_DIR),
                'freq_mhz': freq_mhz,
                'bmaj_deg': bmaj,
                'bmin_deg': bmin,
                'bpa_deg': bpa,
                'naxis1': hdr.get('NAXIS1'),
                'naxis2': hdr.get('NAXIS2'),
            })
        except Exception as e:
            bad_files.append((fpath, str(e)))

    # Sort by frequency
    catalog.sort(key=lambda x: x['freq_mhz'])

    # Write CSV
    catalog_path = os.path.join(IMAGES_DIR, f'channel_catalog_{img_type}.csv')
    with open(catalog_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=catalog[0].keys())
        writer.writeheader()
        writer.writerows(catalog)

    print(f"Catalog written: {catalog_path}")
    print(f"  Channels: {len(catalog)}")
    if catalog:
        print(f"  Freq range: {catalog[0]['freq_mhz']:.3f} - {catalog[-1]['freq_mhz']:.3f} MHz")

        # Check for gaps
        freqs = [c['freq_mhz'] for c in catalog]
        diffs = np.diff(freqs)
        median_step = np.median(diffs)
        gaps = np.where(diffs > 2 * median_step)[0]
        if len(gaps) > 0:
            print(f"  WARNING: {len(gaps)} frequency gaps detected (>{2*median_step:.3f} MHz):")
            for g in gaps[:10]:
                print(f"    Gap at {freqs[g]:.3f} - {freqs[g+1]:.3f} MHz "
                      f"(missing ~{diffs[g]/median_step:.0f} channels)")

    if bad_files:
        print(f"  WARNING: {len(bad_files)} bad files:")
        for fpath, reason in bad_files[:10]:
            print(f"    {fpath}: {reason}")

    return catalog


def main():
    make_cube = '--cube' in sys.argv

    for img_type in ['dirty', 'cleaned']:
        pattern = os.path.join(IMAGES_DIR, 'subband_*', img_type, '*.fits')
        if glob.glob(pattern):
            print(f"\n=== Building catalog for {img_type} images ===")
            catalog = build_catalog(img_type)
        else:
            print(f"\nNo {img_type} images found, skipping.")

    if make_cube:
        print("\n=== FITS cube creation not yet implemented ===")
        print("For the axion search, individual channel FITS files are preferred.")


if __name__ == '__main__':
    main()
