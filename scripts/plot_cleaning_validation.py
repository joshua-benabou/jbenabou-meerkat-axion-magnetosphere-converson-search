#!/usr/bin/env python3
"""
Cleaning validation plots:
1. Side-by-side dirty vs cleaned for several subbands
2. Dirty - Cleaned difference images
3. Model flux and RMS improvement vs frequency
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
PLOT_DIR = '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/plots'
os.makedirs(PLOT_DIR, exist_ok=True)


def load_fits(path):
    with fits.open(path) as h:
        d = h[0].data
    while d.ndim > 2:
        d = d[0]
    return d * 1e3  # mJy


def corner_rms(d):
    corners = [d[:64,:64], d[:64,-64:], d[-64:,:64], d[-64:,-64:]]
    return np.median([np.sqrt(np.mean(c**2)) for c in corners])


# === Plot 1: Side-by-side dirty vs cleaned for 6 subbands ===
sample_sbs = [10, 20, 30, 50, 70, 80]
fig, axes = plt.subplots(len(sample_sbs), 3, figsize=(18, 4 * len(sample_sbs)))

for row, sb in enumerate(sample_sbs):
    dirty_files = sorted(glob.glob(IMG_DIR + '/subband_%03d/chan_0100_*.fits' % sb))
    clean_files = sorted(glob.glob(IMG_DIR + '/subband_%03d/cleaned/chan_0100_*.fits' % sb))
    if not dirty_files or not clean_files:
        continue

    dirty = load_fits(dirty_files[0])
    cleaned = load_fits(clean_files[0])
    diff = dirty - cleaned

    freq = os.path.basename(dirty_files[0]).split('_')[2].replace('MHz.fits', '')
    vmin = np.percentile(dirty, 1)
    vmax = np.percentile(dirty, 99.5)

    im0 = axes[row, 0].imshow(dirty, origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
    axes[row, 0].set_title(r'Dirty (SB %d, %s MHz)' % (sb, freq), fontsize=13)
    plt.colorbar(im0, ax=axes[row, 0], shrink=0.8)

    im1 = axes[row, 1].imshow(cleaned, origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
    axes[row, 1].set_title(r'Cleaned', fontsize=13)
    plt.colorbar(im1, ax=axes[row, 1], shrink=0.8)

    diff_lim = np.percentile(np.abs(diff), 99)
    im2 = axes[row, 2].imshow(diff, origin='lower', cmap='RdBu_r', vmin=-diff_lim, vmax=diff_lim)
    axes[row, 2].set_title(r'Dirty $-$ Cleaned (model)', fontsize=13)
    plt.colorbar(im2, ax=axes[row, 2], shrink=0.8)

    for ax in axes[row, :]:
        ax.set_xlabel(r'Pixel')
        ax.set_ylabel(r'Pixel')

plt.suptitle(r'Cleaning Validation: Dirty vs Cleaned (channel 100 of each subband)', fontsize=16, y=1.0)
plt.tight_layout()
plt.savefig(PLOT_DIR + '/cleaning_validation_grid.png', dpi=120)
plt.close()
print('Saved: cleaning_validation_grid.png')


# === Plot 2: Model flux and RMS vs subband ===
subbands = []
model_flux = []
dirty_peak = []
clean_peak = []
dirty_crms = []
clean_crms = []
freq_centers = []

for sb in range(86):
    log_file = '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/scripts/logs/phase3_clean_%d.out' % sb
    if not os.path.exists(log_file):
        continue
    with open(log_file) as f:
        txt = f.read()
    if 'Exit code: 0' not in txt:
        continue

    # Parse model flux from log
    for line in txt.split('\n'):
        if 'Model flux (total)' in line:
            mf = float(line.split(':')[1].strip().split()[0])
        if 'Freq range' in line:
            parts = line.split(':')[1].strip().split('-')
            f0 = float(parts[0])
            f1 = float(parts[1].replace(' MHz', ''))
            fc = (f0 + f1) / 2

    # Measure from FITS
    dirty_files = sorted(glob.glob(IMG_DIR + '/subband_%03d/chan_0100_*.fits' % sb))
    clean_files = sorted(glob.glob(IMG_DIR + '/subband_%03d/cleaned/chan_0100_*.fits' % sb))
    if not dirty_files or not clean_files:
        continue

    dd = load_fits(dirty_files[0])
    cd = load_fits(clean_files[0])

    subbands.append(sb)
    model_flux.append(mf / 1e3)  # Jy
    freq_centers.append(fc)
    dirty_peak.append(np.nanmax(dd))
    clean_peak.append(np.nanmax(cd))
    dirty_crms.append(corner_rms(dd))
    clean_crms.append(corner_rms(cd))

subbands = np.array(subbands)
model_flux = np.array(model_flux)
freq_centers = np.array(freq_centers)
dirty_peak = np.array(dirty_peak)
clean_peak = np.array(clean_peak)
dirty_crms = np.array(dirty_crms)
clean_crms = np.array(clean_crms)

fig, axes = plt.subplots(2, 2, figsize=(18, 12))

# Model flux vs frequency
axes[0, 0].plot(freq_centers, model_flux, 'o-', color='cornflowerblue', markersize=4)
axes[0, 0].set_xlabel(r'Frequency [MHz]')
axes[0, 0].set_ylabel(r'Total Model Flux [Jy]')
axes[0, 0].set_title(r'Model Flux per Subband')

# Peak flux: dirty vs cleaned
axes[0, 1].plot(freq_centers, dirty_peak, 'o-', color='maroon', markersize=4, label=r'Dirty')
axes[0, 1].plot(freq_centers, clean_peak, 's-', color='forestgreen', markersize=4, label=r'Cleaned')
axes[0, 1].set_xlabel(r'Frequency [MHz]')
axes[0, 1].set_ylabel(r'Peak Flux [mJy/beam]')
axes[0, 1].set_title(r'Peak Flux: Dirty vs Cleaned (chan 100)')
axes[0, 1].legend()

# Corner RMS: dirty vs cleaned
axes[1, 0].plot(freq_centers, dirty_crms, 'o-', color='maroon', markersize=4, label=r'Dirty')
axes[1, 0].plot(freq_centers, clean_crms, 's-', color='forestgreen', markersize=4, label=r'Cleaned')
axes[1, 0].set_xlabel(r'Frequency [MHz]')
axes[1, 0].set_ylabel(r'Corner RMS [mJy/beam]')
axes[1, 0].set_title(r'Corner RMS: Dirty vs Cleaned (chan 100)')
axes[1, 0].legend()

# Fractional improvement
improv = (1 - clean_crms / dirty_crms) * 100
axes[1, 1].bar(freq_centers, improv, width=8, color='cornflowerblue', alpha=0.8)
axes[1, 1].axhline(0, color='k', lw=0.5)
axes[1, 1].set_xlabel(r'Frequency [MHz]')
axes[1, 1].set_ylabel(r'Corner RMS Improvement [\%]')
axes[1, 1].set_title(r'RMS Improvement from Cleaning')

plt.suptitle(r'Cleaning Validation Summary (all subbands, niter=100000, threshold=10\% peak)',
             fontsize=16, y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(PLOT_DIR + '/cleaning_validation_summary.png', dpi=150)
plt.close()
print('Saved: cleaning_validation_summary.png')
print(f'\nMedian RMS improvement: {np.median(improv):.1f}%')
print(f'Mean model flux: {np.mean(model_flux):.1f} Jy')
