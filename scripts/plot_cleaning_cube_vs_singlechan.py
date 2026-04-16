#!/usr/bin/env python3
"""
Compare cube-mode cleaned vs single-channel cleaned images.

Both for subband 30, channel 100 (1158.8 MHz).
- Cube-mode: from production run (images/subband_030/cleaned/)
- Single-channel: from debug test (test_cleaning/debug/sb030_ch0100/)
- Dirty: from dirty imaging (images/subband_030/)

All data loaded from FITS files. No hardcoded values.
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
import glob

IMG_DIR = '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/images'
DEBUG_DIR = '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/test_cleaning/debug/sb030_ch0100'
PLOT_DIR = '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/plots'


def load_fits(path):
    assert os.path.exists(path), f"FITS not found: {path}"
    with fits.open(path) as h:
        d = h[0].data
    while d.ndim > 2:
        d = d[0]
    return d * 1e3  # mJy


def corner_rms(d):
    corners = [d[:64, :64], d[:64, -64:], d[-64:, :64], d[-64:, -64:]]
    return np.median([np.sqrt(np.mean(c**2)) for c in corners])


# Load data
dirty_path = sorted(glob.glob(IMG_DIR + '/subband_030/chan_0100_*.fits'))[0]
cube_path = sorted(glob.glob(IMG_DIR + '/subband_030/cleaned/chan_0100_*.fits'))[0]
debug_path = DEBUG_DIR + '/10pct_peak_image.fits'

dirty = load_fits(dirty_path)
cube_cleaned = load_fits(cube_path)
debug_cleaned = load_fits(debug_path)

# Compute stats
d_peak, d_crms, d_frms = np.max(dirty), corner_rms(dirty), np.std(dirty)
c_peak, c_crms, c_frms = np.max(cube_cleaned), corner_rms(cube_cleaned), np.std(cube_cleaned)
s_peak, s_crms, s_frms = np.max(debug_cleaned), corner_rms(debug_cleaned), np.std(debug_cleaned)

print(f"Dirty:           peak={d_peak:.1f}, corner_rms={d_crms:.1f}, full_rms={d_frms:.1f} mJy")
print(f"Cube cleaned:    peak={c_peak:.1f}, corner_rms={c_crms:.1f}, full_rms={c_frms:.1f} mJy")
print(f"Single-ch clean: peak={s_peak:.1f}, corner_rms={s_crms:.1f}, full_rms={s_frms:.1f} mJy")

assert d_peak > c_peak, f"Cube cleaned peak ({c_peak}) >= dirty peak ({d_peak})"
assert d_peak > s_peak, f"Single-chan cleaned peak ({s_peak}) >= dirty peak ({d_peak})"

# Plot
vmin = np.percentile(dirty, 1)
vmax = np.percentile(dirty, 99.5)

fig, axes = plt.subplots(2, 3, figsize=(20, 13))

im0 = axes[0, 0].imshow(dirty, origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
axes[0, 0].set_title(r'Dirty', fontsize=15)
plt.colorbar(im0, ax=axes[0, 0], label=r'mJy/beam', shrink=0.8)

im1 = axes[0, 1].imshow(cube_cleaned, origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
axes[0, 1].set_title(r'Cube-mode cleaned (niter=100k)', fontsize=15)
plt.colorbar(im1, ax=axes[0, 1], label=r'mJy/beam', shrink=0.8)

im2 = axes[0, 2].imshow(debug_cleaned, origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
axes[0, 2].set_title(r'Single-chan cleaned (niter=10k)', fontsize=15)
plt.colorbar(im2, ax=axes[0, 2], label=r'mJy/beam', shrink=0.8)

diff_cube = dirty - cube_cleaned
diff_debug = dirty - debug_cleaned
diff_lim = max(np.percentile(np.abs(diff_cube), 99), np.percentile(np.abs(diff_debug), 99))

im3 = axes[1, 0].imshow(diff_cube, origin='lower', cmap='RdBu_r', vmin=-diff_lim, vmax=diff_lim)
axes[1, 0].set_title(r'Dirty $-$ Cube-cleaned', fontsize=15)
plt.colorbar(im3, ax=axes[1, 0], label=r'mJy/beam', shrink=0.8)

im4 = axes[1, 1].imshow(diff_debug, origin='lower', cmap='RdBu_r', vmin=-diff_lim, vmax=diff_lim)
axes[1, 1].set_title(r'Dirty $-$ Single-chan cleaned', fontsize=15)
plt.colorbar(im4, ax=axes[1, 1], label=r'mJy/beam', shrink=0.8)

row = np.unravel_index(np.argmax(dirty), dirty.shape)[0]
axes[1, 2].plot(dirty[row, :], color='maroon', lw=1, label='Dirty')
axes[1, 2].plot(cube_cleaned[row, :], color='cornflowerblue', lw=1, label='Cube cleaned')
axes[1, 2].plot(debug_cleaned[row, :], color='forestgreen', lw=1, label='Single-chan cleaned')
axes[1, 2].set_xlabel(r'Pixel')
axes[1, 2].set_ylabel(r'Flux [mJy/beam]')
axes[1, 2].set_title(r'Slice through peak (row %d)' % row, fontsize=15)
axes[1, 2].legend(fontsize=13)

for ax in axes.flat[:5]:
    ax.set_xlabel(r'Pixel')
    ax.set_ylabel(r'Pixel')

stats = (
    r'Subband 30, Channel 100 (1158.8 MHz)' + '\n'
    + r'Dirty: peak=%.1f, corner RMS=%.1f, full RMS=%.1f mJy' % (d_peak, d_crms, d_frms) + '\n'
    + r'Cube cleaned: peak=%.1f, corner RMS=%.1f, full RMS=%.1f mJy' % (c_peak, c_crms, c_frms) + '\n'
    + r'Single-chan cleaned: peak=%.1f, corner RMS=%.1f, full RMS=%.1f mJy' % (s_peak, s_crms, s_frms)
)
fig.text(0.5, 0.01, stats, ha='center', fontsize=13,
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

plt.suptitle(r'Cleaning Comparison: Cube-mode vs Single-channel', fontsize=17, y=0.98)
plt.tight_layout(rect=[0, 0.06, 1, 0.96])
plt.savefig(PLOT_DIR + '/cleaning_cube_vs_singlechan.png', dpi=150)
plt.close()
print("Saved: cleaning_cube_vs_singlechan.png")
