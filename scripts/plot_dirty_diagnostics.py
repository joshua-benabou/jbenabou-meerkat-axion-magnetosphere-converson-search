#!/usr/bin/env python3
"""
Generate all dirty image diagnostic plots from the FITS files on disk.

Produces:
  - dirty_image_example.png: example GC image at ~1159 MHz
  - rms_peak_vs_freq.png: RMS and peak flux vs frequency across all subbands
  - beam_vs_freq.png: synthesized beam size vs frequency
  - cross_subband.png: continuity check at subband boundaries
  - waterfall.png: central pixel intensity vs frequency
  - noise_histogram.png: histogram of per-channel RMS values
  - sgra_zoom.png: zoomed view of Sgr A* region

All data loaded from FITS files in images/subband_XXX/.
No hardcoded values — everything computed at runtime.
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
from astropy.wcs import WCS
import glob

IMG_DIR = '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/images'
PLOT_DIR = '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/plots'
os.makedirs(PLOT_DIR, exist_ok=True)


def load_fits(path):
    """Load FITS, squeeze to 2D, return data in mJy and header."""
    with fits.open(path) as h:
        d = h[0].data.copy()
        hdr = h[0].header.copy()
    while d.ndim > 2:
        d = d[0]
    return d * 1e3, hdr  # mJy


def freq_from_filename(fname):
    """Extract frequency in MHz from filename like chan_0100_1158.766MHz.fits."""
    base = os.path.basename(fname)
    parts = base.split('_')
    freq_str = parts[2].replace('MHz.fits', '')
    return float(freq_str)


def corner_rms_mJy(data_mJy):
    """Measure median RMS in four 64x64 corners."""
    corners = [data_mJy[:64, :64], data_mJy[:64, -64:],
               data_mJy[-64:, :64], data_mJy[-64:, -64:]]
    return np.median([np.sqrt(np.mean(c**2)) for c in corners])


# ============================================================
# Collect per-channel stats across all subbands
# ============================================================
print("Scanning all subbands...", flush=True)
all_freqs = []
all_peaks = []
all_rms = []
all_corner_rms = []
all_bmaj = []
all_bmin = []

for sb in range(86):
    sb_dir = os.path.join(IMG_DIR, 'subband_%03d' % sb)
    fits_files = sorted(glob.glob(sb_dir + '/chan_*.fits'))
    if not fits_files:
        continue

    # Sample every 10th channel, skip first/last 5 (edge artifacts)
    safe_start = min(5, len(fits_files) - 1)
    safe_end = max(0, len(fits_files) - 5)
    sample_indices = list(range(safe_start, safe_end, 10))
    if not sample_indices:
        sample_indices = [len(fits_files) // 2]
    for idx in sample_indices:
        f = fits_files[idx]
        freq = freq_from_filename(f)
        data, hdr = load_fits(f)

        all_freqs.append(freq)
        all_peaks.append(np.nanmax(data))
        all_rms.append(np.std(data))
        all_corner_rms.append(corner_rms_mJy(data))

        if 'BMAJ' in hdr:
            all_bmaj.append(hdr['BMAJ'] * 3600)  # deg to arcsec
            all_bmin.append(hdr['BMIN'] * 3600)
        else:
            all_bmaj.append(np.nan)
            all_bmin.append(np.nan)

    if (sb + 1) % 10 == 0:
        print(f"  Scanned {sb+1}/86 subbands", flush=True)

all_freqs = np.array(all_freqs)
all_peaks = np.array(all_peaks)
all_rms = np.array(all_rms)
all_corner_rms = np.array(all_corner_rms)
all_bmaj = np.array(all_bmaj)
all_bmin = np.array(all_bmin)

sort_idx = np.argsort(all_freqs)
all_freqs = all_freqs[sort_idx]
all_peaks = all_peaks[sort_idx]
all_rms = all_rms[sort_idx]
all_corner_rms = all_corner_rms[sort_idx]
all_bmaj = all_bmaj[sort_idx]
all_bmin = all_bmin[sort_idx]

# Remove extreme outliers (RFI-corrupted channels)
# Flag channels where peak > 10x median peak or corner RMS > 10x median
median_peak = np.median(all_peaks)
median_crms = np.median(all_corner_rms)
good = (all_peaks < 10 * median_peak) & (all_corner_rms < 10 * median_crms) & (all_peaks > 0)

print(f"Total channels sampled: {len(all_freqs)} ({np.sum(~good)} outliers removed)", flush=True)
print(f"Freq range: {all_freqs.min():.3f} - {all_freqs.max():.3f} MHz", flush=True)
print(f"Peak range (good): {all_peaks[good].min():.2f} - {all_peaks[good].max():.2f} mJy", flush=True)
print(f"Corner RMS range (good): {all_corner_rms[good].min():.2f} - {all_corner_rms[good].max():.2f} mJy", flush=True)

# Known RFI bands for shading
RFI_BANDS = [
    (925, 960, 'GSM-900'),
    (1164, 1189, 'GPS/Galileo'),
    (1190, 1215, 'GPS/Galileo'),
    (1525, 1559, 'Inmarsat'),
    (1559, 1610, 'GPS L1'),
    (1610, 1626, 'Iridium'),
    (1670, 1712, 'Band edge'),
    (856, 880, 'Band edge'),
]

# ============================================================
# 1. dirty_image_example.png
# ============================================================
print("\nGenerating dirty_image_example.png...", flush=True)
example_files = sorted(glob.glob(IMG_DIR + '/subband_030/chan_0100_*.fits'))
assert len(example_files) > 0, "No FITS found for subband 030 chan 100"
example_data, example_hdr = load_fits(example_files[0])
example_freq = freq_from_filename(example_files[0])

fig, ax = plt.subplots(1, 1, figsize=(10, 9))
vmin = np.percentile(example_data, 1)
vmax = np.percentile(example_data, 99.5)
im = ax.imshow(example_data, origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
plt.colorbar(im, ax=ax, label=r'mJy/beam')
ax.set_xlabel(r'Pixel')
ax.set_ylabel(r'Pixel')
ax.set_title(r'Dirty Image: Subband 30, Channel 100 (%.1f MHz, $\Delta\nu = 26$ kHz)' % example_freq,
             fontsize=16)

peak = np.nanmax(example_data)
crms = corner_rms_mJy(example_data)
ax.text(0.02, 0.98, r'Peak = %.1f mJy, Corner RMS = %.1f mJy, DR = %.1f' % (peak, crms, peak/crms),
        transform=ax.transAxes, va='top', fontsize=14,
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig(PLOT_DIR + '/dirty_image_example.png', dpi=150)
plt.close()
print(f"  Peak={peak:.1f} mJy, Corner RMS={crms:.1f} mJy, DR={peak/crms:.1f}")
print("  Saved: dirty_image_example.png")

# ============================================================
# 2. rms_peak_vs_freq.png
# ============================================================
print("\nGenerating rms_peak_vs_freq.png...", flush=True)
fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True)

# Shade known RFI bands
for ax in axes:
    for flo, fhi, rfi_label in RFI_BANDS:
        ax.axvspan(flo, fhi, alpha=0.15, color='red')

axes[0].plot(all_freqs[good], all_peaks[good], 'o-', color='cornflowerblue',
             markersize=3, lw=0.5, alpha=0.7)
axes[0].set_ylabel(r'Peak Flux [mJy/beam]')
axes[0].set_title(r'Peak Flux vs Frequency (dirty images, mid-channel samples)')

axes[1].plot(all_freqs[good], all_corner_rms[good], 'o-', color='forestgreen',
             markersize=3, lw=0.5, alpha=0.7, label=r'Corner RMS')
axes[1].plot(all_freqs[good], all_rms[good], 'o-', color='cornflowerblue',
             markersize=3, lw=0.5, alpha=0.3, label=r'Full-image RMS')
axes[1].axhline(0.24, color='maroon', ls='--', lw=1.5, label=r'Theoretical (0.24 mJy)')
axes[1].set_xlabel(r'Frequency [MHz]')
axes[1].set_ylabel(r'RMS [mJy/beam]')
axes[1].set_title(r'RMS Noise vs Frequency (red shading = known RFI bands)')
axes[1].legend(fontsize=14)

plt.tight_layout()
plt.savefig(PLOT_DIR + '/rms_peak_vs_freq.png', dpi=150)
plt.close()
print(f"  Median corner RMS: {np.nanmedian(all_corner_rms[good]):.2f} mJy")
print(f"  Median full RMS: {np.nanmedian(all_rms[good]):.2f} mJy")
print("  Saved: rms_peak_vs_freq.png")

# ============================================================
# 3. beam_vs_freq.png
# ============================================================
print("\nGenerating beam_vs_freq.png...", flush=True)
valid_beam = ~np.isnan(all_bmaj)
assert np.any(valid_beam), "No valid beam info in any FITS header"

fig, ax = plt.subplots(1, 1, figsize=(16, 6))
ax.plot(all_freqs[valid_beam], all_bmaj[valid_beam], '.', color='cornflowerblue',
        markersize=2, alpha=0.7, label=r'BMAJ')
ax.plot(all_freqs[valid_beam], all_bmin[valid_beam], '.', color='maroon',
        markersize=2, alpha=0.7, label=r'BMIN')
ax.set_xlabel(r'Frequency [MHz]')
ax.set_ylabel(r'Beam Size [arcsec]')
ax.set_title(r'Synthesized Beam vs Frequency')
ax.legend(fontsize=14)

plt.tight_layout()
plt.savefig(PLOT_DIR + '/beam_vs_freq.png', dpi=150)
plt.close()
print(f"  BMAJ range: {np.nanmin(all_bmaj[valid_beam]):.2f} - {np.nanmax(all_bmaj[valid_beam]):.2f} arcsec")
print(f"  BMIN range: {np.nanmin(all_bmin[valid_beam]):.2f} - {np.nanmax(all_bmin[valid_beam]):.2f} arcsec")
print("  Saved: beam_vs_freq.png")

# ============================================================
# 4. cross_subband.png
# ============================================================
print("\nGenerating cross_subband.png...", flush=True)
fig, axes = plt.subplots(2, 1, figsize=(16, 10))

# Pick a few subband boundaries and compare last/first channels
boundary_sbs = [10, 20, 30, 40, 50]
colors = ['cornflowerblue', 'forestgreen', 'maroon', 'goldenrod', 'firebrick']

for i, sb in enumerate(boundary_sbs):
    # Last channel of subband sb
    files_a = sorted(glob.glob(IMG_DIR + '/subband_%03d/chan_*.fits' % sb))
    # First channel of subband sb+1
    files_b = sorted(glob.glob(IMG_DIR + '/subband_%03d/chan_*.fits' % (sb + 1)))
    if not files_a or not files_b:
        continue

    data_a, _ = load_fits(files_a[-1])
    data_b, _ = load_fits(files_b[0])
    freq_a = freq_from_filename(files_a[-1])
    freq_b = freq_from_filename(files_b[0])

    # Central row slice
    row = data_a.shape[0] // 2
    axes[0].plot(data_a[row, :], color=colors[i], lw=0.8, alpha=0.8,
                 label=r'SB %d last (%.1f MHz)' % (sb, freq_a))
    axes[0].plot(data_b[row, :], '--', color=colors[i], lw=0.8, alpha=0.8,
                 label=r'SB %d first (%.1f MHz)' % (sb + 1, freq_b))

    # Peak difference
    diff = np.abs(np.nanmax(data_a) - np.nanmax(data_b))
    axes[1].bar(freq_a, diff, width=5, color=colors[i], alpha=0.7)

axes[0].set_xlabel(r'Pixel')
axes[0].set_ylabel(r'Flux [mJy/beam]')
axes[0].set_title(r'Cross-Subband Continuity: Central Row Slices at Boundaries')
axes[0].legend(fontsize=10, ncol=2)

axes[1].set_xlabel(r'Boundary Frequency [MHz]')
axes[1].set_ylabel(r'$|\Delta$ Peak$|$ [mJy]')
axes[1].set_title(r'Peak Flux Difference at Subband Boundaries')

plt.tight_layout()
plt.savefig(PLOT_DIR + '/cross_subband.png', dpi=150)
plt.close()
print("  Saved: cross_subband.png")

# ============================================================
# 5. waterfall.png
# ============================================================
print("\nGenerating waterfall.png...", flush=True)
center_pix = (256, 256)
center_vals = []
center_freqs = []

for sb in range(86):
    sb_dir = os.path.join(IMG_DIR, 'subband_%03d' % sb)
    fits_files = sorted(glob.glob(sb_dir + '/chan_*.fits'))
    safe_s = min(5, len(fits_files) - 1)
    safe_e = max(0, len(fits_files) - 5)
    for idx in range(safe_s, safe_e, 5):  # every 5th channel, skip edges
        f = fits_files[idx]
        data, _ = load_fits(f)
        center_vals.append(data[center_pix[0], center_pix[1]])
        center_freqs.append(freq_from_filename(f))

center_freqs = np.array(center_freqs)
center_vals = np.array(center_vals)
sort_idx = np.argsort(center_freqs)
center_freqs = center_freqs[sort_idx]
center_vals = center_vals[sort_idx]

# Clip extreme RFI outliers for plotting
med_cv = np.median(center_vals)
clip_lo = np.percentile(center_vals, 1)
clip_hi = np.percentile(center_vals, 99)

fig, ax = plt.subplots(1, 1, figsize=(16, 6))
for flo, fhi, rfi_label in RFI_BANDS:
    ax.axvspan(flo, fhi, alpha=0.15, color='red')
ax.plot(center_freqs, np.clip(center_vals, clip_lo, clip_hi), '-',
        color='cornflowerblue', lw=0.5, alpha=0.8)
ax.set_ylim(clip_lo - 10, clip_hi + 10)
ax.set_xlabel(r'Frequency [MHz]')
ax.set_ylabel(r'Flux at Central Pixel [mJy/beam]')
ax.set_title(r'Waterfall: Central Pixel (256, 256) vs Frequency (clipped at 1-99 pct)')

plt.tight_layout()
plt.savefig(PLOT_DIR + '/waterfall.png', dpi=150)
plt.close()
print(f"  Central pixel flux range: {center_vals.min():.2f} - {center_vals.max():.2f} mJy")
print("  Saved: waterfall.png")

# ============================================================
# 6. noise_histogram.png
# ============================================================
print("\nGenerating noise_histogram.png...", flush=True)
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

axes[0].hist(all_corner_rms, bins=50, color='cornflowerblue', alpha=0.8)
axes[0].axvline(np.median(all_corner_rms), color='maroon', ls='--', lw=2,
                label=r'Median = %.1f mJy' % np.median(all_corner_rms))
axes[0].axvline(0.24, color='forestgreen', ls='--', lw=1.5,
                label=r'Theoretical = 0.24 mJy')
axes[0].set_xlabel(r'Corner RMS [mJy/beam]')
axes[0].set_ylabel(r'Count')
axes[0].set_title(r'Corner RMS Distribution')
axes[0].legend(fontsize=13)

axes[1].hist(all_rms, bins=50, color='goldenrod', alpha=0.8)
axes[1].axvline(np.median(all_rms), color='maroon', ls='--', lw=2,
                label=r'Median = %.1f mJy' % np.median(all_rms))
axes[1].set_xlabel(r'Full-Image RMS [mJy/beam]')
axes[1].set_ylabel(r'Count')
axes[1].set_title(r'Full-Image RMS Distribution')
axes[1].legend(fontsize=13)

plt.tight_layout()
plt.savefig(PLOT_DIR + '/noise_histogram.png', dpi=150)
plt.close()
print("  Saved: noise_histogram.png")

# ============================================================
# 7. sgra_zoom.png
# ============================================================
print("\nGenerating sgra_zoom.png...", flush=True)
# Sgr A* should be near the image center (pointing center)
# Zoom to central 128x128 pixels
data_zoom = example_data[192:320, 192:320]
fig, ax = plt.subplots(1, 1, figsize=(10, 9))
vmin_z = np.percentile(data_zoom, 1)
vmax_z = np.percentile(data_zoom, 99.5)
im = ax.imshow(data_zoom, origin='lower', cmap='inferno', vmin=vmin_z, vmax=vmax_z,
               extent=[192, 320, 192, 320])
plt.colorbar(im, ax=ax, label=r'mJy/beam')
ax.set_xlabel(r'Pixel')
ax.set_ylabel(r'Pixel')
ax.set_title(r'Sgr A* Region (central 128$\times$128 pixels, %.1f MHz)' % example_freq, fontsize=16)

plt.tight_layout()
plt.savefig(PLOT_DIR + '/sgra_zoom.png', dpi=150)
plt.close()
print("  Saved: sgra_zoom.png")

print("\nAll plots generated.", flush=True)
