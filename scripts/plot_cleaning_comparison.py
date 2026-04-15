#!/usr/bin/env python3
"""
Side-by-side comparison of dirty vs cleaned images from the debug cleaning test.
Shows dirty, best cleaned (10% peak threshold), residual, and model.
"""

import os
os.environ["PATH"] = '/global/software/sl-7.x86_64/modules/tools/texlive/2016/bin/x86_64-linux/:' + os.environ["PATH"]

import sys
sys.path.insert(0, '/global/scratch/projects/pc_heptheory/jbenabou')
import plotting_defaults

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits

DATA_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/test_cleaning/debug/sb030_ch0100'
)
OUT_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/plots'
)


def load_fits(name):
    with fits.open(f'{DATA_DIR}/{name}') as hdul:
        data = hdul[0].data
    while data.ndim > 2:
        data = data[0]
    return data * 1e3  # Jy -> mJy


def main():
    dirty = load_fits('dirty_image.fits')
    cleaned = load_fits('10pct_peak_image.fits')
    residual = load_fits('10pct_peak_residual.fits')
    model = load_fits('10pct_peak_model.fits')

    vmin = np.percentile(dirty, 1)
    vmax = np.percentile(dirty, 99.5)

    fig, axes = plt.subplots(2, 2, figsize=(14, 13))

    # Dirty
    im0 = axes[0, 0].imshow(dirty, origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
    axes[0, 0].set_title(r'Dirty Image ($n_{\rm iter}=0$)', fontsize=18)
    plt.colorbar(im0, ax=axes[0, 0], label=r'mJy/beam', shrink=0.85)

    # Cleaned
    im1 = axes[0, 1].imshow(cleaned, origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
    axes[0, 1].set_title(r'Cleaned ($10\%$ peak threshold)', fontsize=18)
    plt.colorbar(im1, ax=axes[0, 1], label=r'mJy/beam', shrink=0.85)

    # Residual
    resid_vmax = np.percentile(np.abs(residual), 99)
    im2 = axes[1, 0].imshow(residual, origin='lower', cmap='RdBu_r',
                             vmin=-resid_vmax, vmax=resid_vmax)
    axes[1, 0].set_title(r'Residual', fontsize=18)
    plt.colorbar(im2, ax=axes[1, 0], label=r'mJy/beam', shrink=0.85)

    # Model
    model_nonzero = model[model > 0]
    model_vmax = np.percentile(model_nonzero, 99) if len(model_nonzero) > 0 else 1
    im3 = axes[1, 1].imshow(model, origin='lower', cmap='inferno',
                             vmin=0, vmax=model_vmax)
    axes[1, 1].set_title(r'Model (%d pixels)' % int(np.count_nonzero(model)), fontsize=18)
    plt.colorbar(im3, ax=axes[1, 1], label=r'mJy/beam', shrink=0.85)

    for ax in axes.flat:
        ax.set_xlabel(r'Pixel')
        ax.set_ylabel(r'Pixel')

    # Stats annotation
    dirty_rms_corner = 40.06
    clean_rms_corner = 26.44
    peak_dirty = np.max(dirty)
    peak_clean = np.max(cleaned)
    model_flux = np.sum(model)

    stats_text = (
        r'Subband 30, Channel 100 (1158.8 MHz), $\Delta\nu = 26$ kHz'
        '\n'
        r'Dirty: peak $= %.1f$ mJy, corner RMS $= %.1f$ mJy' % (peak_dirty, dirty_rms_corner)
        + '\n'
        + r'Cleaned: peak $= %.1f$ mJy, corner RMS $= %.1f$ mJy' % (peak_clean, clean_rms_corner)
        + '\n'
        + r'Model flux: $%.0f$ mJy in %d pixels' % (model_flux, int(np.count_nonzero(model)))
        + '\n'
        + r'RMS reduction: $%.1f \to %.1f$ mJy ($%d\%%$)' % (
            dirty_rms_corner, clean_rms_corner,
            int((1 - clean_rms_corner / dirty_rms_corner) * 100))
    )
    fig.text(0.5, 0.01, stats_text, ha='center', fontsize=14,
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

    plt.suptitle(r'\texttt{tclean} Debug: Dirty vs Cleaned (single 26 kHz channel, GC field)',
                 fontsize=18, y=0.98)
    plt.tight_layout(rect=[0, 0.10, 1, 0.96])

    os.makedirs(OUT_DIR, exist_ok=True)
    outpath = f'{OUT_DIR}/cleaning_comparison_sb030_ch100.png'
    plt.savefig(outpath, dpi=150)
    print(f'Saved: {outpath}')
    plt.close()

    # === Multi-strategy comparison ===
    fig2, axes2 = plt.subplots(2, 3, figsize=(20, 12))
    strategies = [
        ('dirty_image.fits', r'Dirty ($n_{\rm iter}=0$)'),
        ('3sig_corner_image.fits', r'$3\sigma$ corner (120 mJy)'),
        ('5sig_corner_image.fits', r'$5\sigma$ corner (200 mJy)'),
        ('10pct_peak_image.fits', r'$10\%$ peak (27 mJy)'),
        ('1pct_peak_image.fits', r'$1\%$ peak (2.7 mJy)'),
        ('auto_multithresh_image.fits', r'auto-multithresh'),
    ]

    for ax, (fname, title) in zip(axes2.flat, strategies):
        data = load_fits(fname)
        im = ax.imshow(data, origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
        ax.set_title(title, fontsize=16)
        ax.set_xlabel(r'Pixel')
        ax.set_ylabel(r'Pixel')
        plt.colorbar(im, ax=ax, label=r'mJy/beam', shrink=0.8)

    plt.suptitle(r'Cleaning Strategy Comparison (subband 30, chan 100, 1158.8 MHz)',
                 fontsize=18, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    outpath2 = f'{OUT_DIR}/cleaning_strategies_comparison.png'
    plt.savefig(outpath2, dpi=150)
    print(f'Saved: {outpath2}')
    plt.close()


if __name__ == '__main__':
    main()
