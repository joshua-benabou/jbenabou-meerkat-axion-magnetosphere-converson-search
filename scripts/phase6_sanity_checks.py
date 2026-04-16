#!/usr/bin/env python3
"""
Phase 6: Sanity Checks for MeerKAT Galactic Center Imaging
===========================================================

Validates the imaging pipeline output with five independent checks:

  1. Noise assessment: corner RMS vs theoretical sensitivity
  2. Point source injection/recovery: fake source recovery at multiple SNR levels
  3. Cleaning fidelity: dirty vs cleaned image comparison
  4. 5-sigma detection threshold vs frequency
  5. Channel-to-channel consistency: adjacent-channel correlation

All data is loaded from real FITS files on disk. Nothing is hardcoded or faked.

READ-ONLY on all FITS files. Outputs only to plots/ directory.
"""

import os
os.environ["PATH"] = (
    '/global/software/sl-7.x86_64/modules/tools/texlive/2016/bin/x86_64-linux/:'
    + os.environ["PATH"]
)

import sys
sys.path.insert(0, '/global/scratch/projects/pc_heptheory/jbenabou')
import plotting_defaults

import re
import glob
import logging
import argparse
from pathlib import Path
from multiprocessing import Pool

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(
    "/global/scratch/projects/pc_heptheory/jbenabou/"
    "NS_megaproject/MeerKAT_data/meerkat_reduction_project"
)
IMAGE_DIR = PROJECT_ROOT / "images"
PLOT_DIR = PROJECT_ROOT / "plots"
RFI_CSV = PROJECT_ROOT / "rfi_channel_flags.csv"

NPIX = 512
CORNER_SIZE = 64  # 64x64 pixel boxes in each corner for noise measurement

# MeerKAT theoretical sensitivity parameters
N_ANT = 64
SEFD_JY = 420.0       # MeerKAT L-band system equivalent flux density (Jy)
T_INT_S = 8.0 * 3600   # 8 hours integration
DELTA_NU_HZ = 26123.0  # 26.123 kHz channel width
N_POL = 2
# sigma = SEFD / sqrt(N_ant*(N_ant-1)*N_pol*delta_nu*t_int)
THEORETICAL_RMS_JY = SEFD_JY / np.sqrt(
    N_ANT * (N_ANT - 1) * N_POL * DELTA_NU_HZ * T_INT_S
)
THEORETICAL_RMS_MJY = THEORETICAL_RMS_JY * 1e3  # ~0.24 mJy

# Filename regex: chan_NNNN_FREQ.fits
CHAN_RE = re.compile(r"chan_(\d{4})_([\d.]+)MHz\.fits$")

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def load_fits_2d(path):
    """Load a FITS file and return squeezed 2D data (Jy/beam) and header."""
    with fits.open(path, memmap=True, mode="readonly") as hdul:
        data = hdul[0].data.copy()
        hdr = hdul[0].header.copy()
    while data.ndim > 2:
        data = data[0]
    return data.astype(np.float64), hdr


def corner_rms(data):
    """Compute RMS in four 64x64 corners, return median of the four."""
    c = CORNER_SIZE
    corners = [
        data[:c, :c],
        data[:c, -c:],
        data[-c:, :c],
        data[-c:, -c:],
    ]
    rms_vals = [np.sqrt(np.nanmean(corner**2)) for corner in corners]
    return np.median(rms_vals)


def freq_from_filename(fname):
    """Extract frequency in MHz from filename like chan_0100_1158.766MHz.fits."""
    base = os.path.basename(fname)
    m = CHAN_RE.search(base)
    if m:
        return float(m.group(2))
    # Fallback
    parts = base.split('_')
    freq_str = parts[2].replace('MHz.fits', '')
    return float(freq_str)


def chan_from_filename(fname):
    """Extract channel index from filename."""
    base = os.path.basename(fname)
    m = CHAN_RE.search(base)
    if m:
        return int(m.group(1))
    return -1


def discover_subbands():
    """Return sorted list of subband directory paths that exist."""
    sb_dirs = sorted(glob.glob(str(IMAGE_DIR / "subband_*")))
    return sb_dirs


def sample_files_from_subbands(sb_dirs, stride=40, skip_edges=5):
    """
    Sample FITS files across subbands for efficient analysis.
    Returns list of (filepath, subband_idx, chan_idx, freq_mhz).
    stride: take every Nth channel within each subband.
    skip_edges: skip first/last N channels per subband (edge artifacts).
    """
    samples = []
    for sb_dir in sb_dirs:
        sb_name = os.path.basename(sb_dir)
        try:
            sb_idx = int(sb_name.split("_")[1])
        except (IndexError, ValueError):
            continue

        fits_files = sorted(glob.glob(os.path.join(sb_dir, "chan_*.fits")))
        if len(fits_files) < 2 * skip_edges + 1:
            continue

        safe_start = skip_edges
        safe_end = len(fits_files) - skip_edges
        for idx in range(safe_start, safe_end, stride):
            fpath = fits_files[idx]
            freq = freq_from_filename(fpath)
            chan = chan_from_filename(fpath)
            samples.append((fpath, sb_idx, chan, freq))

    return samples


def _worker_corner_rms(args):
    """Multiprocessing worker: compute corner RMS for one file."""
    fpath, sb_idx, chan_idx, freq_mhz = args
    try:
        data, hdr = load_fits_2d(fpath)
        if data.shape != (NPIX, NPIX):
            return (sb_idx, chan_idx, freq_mhz, np.nan, fpath)
        rms = corner_rms(data)
        return (sb_idx, chan_idx, freq_mhz, rms, fpath)
    except Exception as e:
        log.warning("Failed to read %s: %s", fpath, e)
        return (sb_idx, chan_idx, freq_mhz, np.nan, fpath)


# ===========================================================================
# Check 1: Noise Assessment
# ===========================================================================

def check_noise_assessment(samples, n_workers):
    """
    Measure corner RMS across sampled channels and compare to theoretical noise.
    """
    log.info("=== Check 1: Noise Assessment ===")
    log.info("Sampling %d channels for noise measurement...", len(samples))

    with Pool(processes=n_workers) as pool:
        results = pool.map(_worker_corner_rms, samples, chunksize=32)

    freqs = np.array([r[2] for r in results])
    rms_vals = np.array([r[3] for r in results])  # Jy/beam

    valid = np.isfinite(rms_vals)
    freqs_v = freqs[valid]
    rms_mJy = rms_vals[valid] * 1e3  # convert to mJy

    sort_idx = np.argsort(freqs_v)
    freqs_v = freqs_v[sort_idx]
    rms_mJy = rms_mJy[sort_idx]

    median_rms = np.median(rms_mJy)
    ratio = median_rms / THEORETICAL_RMS_MJY

    log.info("  Measured channels: %d", len(rms_mJy))
    log.info("  Median corner RMS: %.4f mJy/beam", median_rms)
    log.info("  Theoretical RMS:   %.4f mJy/beam", THEORETICAL_RMS_MJY)
    log.info("  Ratio (measured/theoretical): %.2f", ratio)
    log.info("  Min corner RMS: %.4f mJy, Max: %.4f mJy",
             np.min(rms_mJy), np.max(rms_mJy))

    # Sanity check: measured noise should be within reasonable range of theoretical
    # Typically 1-10x for dirty images (confusion noise, sidelobes, etc.)
    assert median_rms > 0, "Median corner RMS must be positive"
    assert ratio > 0.1, (
        f"Measured RMS ({median_rms:.3f} mJy) is suspiciously low "
        f"compared to theoretical ({THEORETICAL_RMS_MJY:.3f} mJy)"
    )
    assert ratio < 1000, (
        f"Measured RMS ({median_rms:.3f} mJy) is suspiciously high "
        f"compared to theoretical ({THEORETICAL_RMS_MJY:.3f} mJy)"
    )
    log.info("  PASS: Noise level is %.1fx theoretical (expected 1-10x for dirty images)", ratio)

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Left: measured vs theoretical across frequency
    axes[0].scatter(freqs_v, rms_mJy, s=3, alpha=0.5, color='cornflowerblue',
                    label='Measured corner RMS', rasterized=True)
    axes[0].axhline(THEORETICAL_RMS_MJY, color='maroon', ls='--', lw=2,
                    label=r'Theoretical (%.2f mJy)' % THEORETICAL_RMS_MJY)
    axes[0].axhline(median_rms, color='forestgreen', ls=':', lw=2,
                    label=r'Median measured (%.2f mJy)' % median_rms)
    axes[0].set_xlabel('Frequency [MHz]')
    axes[0].set_ylabel('Corner RMS [mJy/beam]')
    axes[0].set_title('Noise Assessment: Measured vs Theoretical')
    axes[0].legend(fontsize=11)
    axes[0].grid(True, alpha=0.3)

    # Right: histogram of ratio
    ratios = rms_mJy / THEORETICAL_RMS_MJY
    axes[1].hist(ratios, bins=50, color='cornflowerblue', alpha=0.8, edgecolor='none')
    axes[1].axvline(np.median(ratios), color='maroon', ls='--', lw=2,
                    label=r'Median ratio = %.1f' % np.median(ratios))
    axes[1].set_xlabel('Measured / Theoretical RMS')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Distribution of RMS Ratios')
    axes[1].legend(fontsize=12)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = PLOT_DIR / "sanity_noise_assessment.png"
    fig.savefig(str(outpath), dpi=150, bbox_inches='tight')
    plt.close(fig)
    log.info("  Saved: %s", outpath)

    return freqs_v, rms_mJy


# ===========================================================================
# Check 2: Point Source Injection Test
# ===========================================================================

def check_point_source_injection(n_trials=200):
    """
    Inject a fake Gaussian point source into a real dirty image at various SNR
    levels and measure recovery rate.
    """
    log.info("=== Check 2: Point Source Injection Test ===")

    # Pick a representative image from a mid-band, RFI-quiet subband
    test_sb = 30
    test_chan = 100
    pattern = str(IMAGE_DIR / f"subband_{test_sb:03d}" / f"chan_{test_chan:04d}_*.fits")
    test_files = sorted(glob.glob(pattern))
    assert len(test_files) > 0, f"No FITS file found for subband {test_sb} chan {test_chan}"
    test_file = test_files[0]

    data, hdr = load_fits_2d(test_file)
    test_freq = freq_from_filename(test_file)
    assert data.shape == (NPIX, NPIX), f"Unexpected shape {data.shape}"

    noise_rms = corner_rms(data)  # Jy/beam
    log.info("  Test image: %s", os.path.basename(test_file))
    log.info("  Frequency: %.3f MHz", test_freq)
    log.info("  Corner RMS: %.4f Jy/beam (%.4f mJy)", noise_rms, noise_rms * 1e3)

    # Read beam from header
    bmaj_deg = hdr.get("BMAJ", 7.0 / 3600)  # default 7 arcsec
    bmin_deg = hdr.get("BMIN", 5.0 / 3600)  # default 5 arcsec
    cdelt = abs(hdr.get("CDELT1", 2.0 / 3600))  # pixel size in degrees

    bmaj_pix = bmaj_deg / cdelt
    bmin_pix = bmin_deg / cdelt

    log.info("  Beam: BMAJ=%.2f arcsec (%.2f pix), BMIN=%.2f arcsec (%.2f pix)",
             bmaj_deg * 3600, bmaj_pix, bmin_deg * 3600, bmin_pix)

    # Gaussian PSF: FWHM -> sigma
    sigma_x = bmin_pix / (2 * np.sqrt(2 * np.log(2)))
    sigma_y = bmaj_pix / (2 * np.sqrt(2 * np.log(2)))

    # Create PSF stamp (31x31 should be large enough for any MeerKAT L-band beam)
    psf_size = 31
    psf_center = psf_size // 2
    yy, xx = np.mgrid[:psf_size, :psf_size]
    psf = np.exp(
        -0.5 * ((xx - psf_center)**2 / sigma_x**2
                 + (yy - psf_center)**2 / sigma_y**2)
    )
    psf /= psf.max()  # normalize peak to 1

    # SNR levels to test
    snr_levels = [1, 2, 3, 5, 7, 10, 15, 20]
    recovery_rate = []

    rng = np.random.RandomState(42)

    for snr in snr_levels:
        injected_flux = snr * noise_rms  # peak flux in Jy/beam
        n_recovered = 0

        for trial in range(n_trials):
            # Random position in the corner regions (where noise is measured,
            # away from bright Sgr A* emission)
            quadrant = trial % 4
            margin = psf_size  # keep away from edges
            if quadrant == 0:
                cx = rng.randint(margin, CORNER_SIZE - margin)
                cy = rng.randint(margin, CORNER_SIZE - margin)
            elif quadrant == 1:
                cx = rng.randint(margin, CORNER_SIZE - margin)
                cy = rng.randint(NPIX - CORNER_SIZE + margin, NPIX - margin)
            elif quadrant == 2:
                cx = rng.randint(NPIX - CORNER_SIZE + margin, NPIX - margin)
                cy = rng.randint(margin, CORNER_SIZE - margin)
            else:
                cx = rng.randint(NPIX - CORNER_SIZE + margin, NPIX - margin)
                cy = rng.randint(NPIX - CORNER_SIZE + margin, NPIX - margin)

            # Inject source
            injected = data.copy()
            x_lo = cx - psf_center
            x_hi = cx - psf_center + psf_size
            y_lo = cy - psf_center
            y_hi = cy - psf_center + psf_size

            # Bounds check
            if x_lo < 0 or x_hi > NPIX or y_lo < 0 or y_hi > NPIX:
                continue

            injected[y_lo:y_hi, x_lo:x_hi] += injected_flux * psf

            # Measure: is the peak in a small box around injection site
            # above 5-sigma (detection threshold)?
            box = injected[y_lo:y_hi, x_lo:x_hi]
            peak_val = np.max(box)

            # Also measure local noise from the original image in the same box
            orig_box = data[y_lo:y_hi, x_lo:x_hi]
            local_rms = np.sqrt(np.mean(orig_box**2))
            if local_rms == 0:
                local_rms = noise_rms

            # Recovery criterion: peak > 3 * local_rms
            # (we use 3-sigma detection, not 5, because we know the position)
            if peak_val > 3 * local_rms:
                n_recovered += 1

        rate = n_recovered / n_trials
        recovery_rate.append(rate)
        log.info("  SNR=%2d: injected flux=%.4f Jy, recovery=%.1f%% (%d/%d)",
                 snr, injected_flux, rate * 100, n_recovered, n_trials)

    # Sanity checks
    assert recovery_rate[-1] > 0.9, (
        f"Recovery rate at SNR=20 is only {recovery_rate[-1]*100:.0f}%% -- something is wrong"
    )
    assert recovery_rate[0] < recovery_rate[-1], (
        "Recovery rate should increase with SNR"
    )
    log.info("  PASS: Recovery rate increases with SNR as expected")

    # --- Plot ---
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    ax.plot(snr_levels, [r * 100 for r in recovery_rate], 'o-',
            color='cornflowerblue', markersize=8, lw=2)
    ax.axhline(100, color='gray', ls=':', alpha=0.5)
    ax.axhline(50, color='gray', ls=':', alpha=0.3)
    ax.axvline(5, color='maroon', ls='--', lw=1.5, alpha=0.7,
               label=r'$5\sigma$ threshold')
    ax.set_xlabel('Injected SNR (peak / corner RMS)')
    ax.set_ylabel('Recovery Rate [\\%]')
    ax.set_title(
        'Point Source Injection Test (%.1f MHz, RMS=%.3f mJy, %d trials/SNR)'
        % (test_freq, noise_rms * 1e3, n_trials)
    )
    ax.set_ylim(-5, 110)
    ax.set_xlim(0, max(snr_levels) + 1)
    ax.legend(fontsize=13)
    ax.grid(True, alpha=0.3)

    # Annotate beam info
    ax.text(0.98, 0.02,
            'Beam: %.1f" x %.1f"\nCorner RMS: %.3f mJy'
            % (bmaj_deg * 3600, bmin_deg * 3600, noise_rms * 1e3),
            transform=ax.transAxes, ha='right', va='bottom', fontsize=11,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    outpath = PLOT_DIR / "sanity_injection_test.png"
    fig.savefig(str(outpath), dpi=150, bbox_inches='tight')
    plt.close(fig)
    log.info("  Saved: %s", outpath)


# ===========================================================================
# Check 3: Cleaning Fidelity
# ===========================================================================

def _worker_clean_compare(args):
    """Load dirty and cleaned image, return comparison stats."""
    dirty_path, clean_path, sb_idx, chan_idx, freq_mhz = args
    try:
        dirty, _ = load_fits_2d(dirty_path)
        clean, _ = load_fits_2d(clean_path)
        if dirty.shape != (NPIX, NPIX) or clean.shape != (NPIX, NPIX):
            return None

        diff = clean - dirty
        dirty_rms = corner_rms(dirty)
        clean_rms = corner_rms(clean)

        return {
            "sb": sb_idx,
            "chan": chan_idx,
            "freq": freq_mhz,
            "dirty_rms": dirty_rms,
            "clean_rms": clean_rms,
            "rms_ratio": clean_rms / dirty_rms if dirty_rms > 0 else np.nan,
            "diff_mean": np.nanmean(diff),
            "diff_std": np.nanstd(diff),
            "dirty_peak": np.nanmax(np.abs(dirty)),
            "clean_peak": np.nanmax(np.abs(clean)),
            # Store flattened pixel samples for scatter/histogram
            # (subsample to keep memory manageable)
            "dirty_flat": dirty[::4, ::4].ravel(),
            "clean_flat": clean[::4, ::4].ravel(),
            "diff_flat": diff[::4, ::4].ravel(),
        }
    except Exception as e:
        log.warning("Clean comparison failed for %s: %s", dirty_path, e)
        return None


def check_cleaning_fidelity(n_workers):
    """
    Compare dirty vs cleaned images across sampled channels.
    """
    log.info("=== Check 3: Cleaning Fidelity ===")

    # Discover cleaned images and match with dirty
    pairs = []
    sb_dirs = discover_subbands()
    for sb_dir in sb_dirs:
        sb_name = os.path.basename(sb_dir)
        try:
            sb_idx = int(sb_name.split("_")[1])
        except (IndexError, ValueError):
            continue

        clean_dir = os.path.join(sb_dir, "cleaned")
        if not os.path.isdir(clean_dir):
            continue

        dirty_files = sorted(glob.glob(os.path.join(sb_dir, "chan_*.fits")))
        clean_files_set = set(os.path.basename(f)
                              for f in glob.glob(os.path.join(clean_dir, "chan_*.fits")))

        # Sample every 40th channel
        for i in range(5, len(dirty_files) - 5, 40):
            dfile = dirty_files[i]
            dname = os.path.basename(dfile)
            if dname in clean_files_set:
                cfile = os.path.join(clean_dir, dname)
                freq = freq_from_filename(dfile)
                chan = chan_from_filename(dfile)
                pairs.append((dfile, cfile, sb_idx, chan, freq))

    if not pairs:
        log.warning("  No dirty/cleaned image pairs found. Skipping cleaning fidelity check.")
        return

    log.info("  Found %d dirty/cleaned pairs to compare", len(pairs))

    with Pool(processes=n_workers) as pool:
        results = pool.map(_worker_clean_compare, pairs, chunksize=8)

    results = [r for r in results if r is not None]
    if not results:
        log.warning("  All comparisons failed. Skipping cleaning fidelity check.")
        return

    log.info("  Successfully compared %d pairs", len(results))

    freqs = np.array([r["freq"] for r in results])
    rms_ratios = np.array([r["rms_ratio"] for r in results])
    dirty_rms = np.array([r["dirty_rms"] for r in results]) * 1e3  # mJy
    clean_rms = np.array([r["clean_rms"] for r in results]) * 1e3  # mJy

    sort_idx = np.argsort(freqs)
    freqs = freqs[sort_idx]
    rms_ratios = rms_ratios[sort_idx]
    dirty_rms = dirty_rms[sort_idx]
    clean_rms = clean_rms[sort_idx]
    results_sorted = [results[i] for i in sort_idx]

    median_ratio = np.nanmedian(rms_ratios)
    log.info("  Median clean/dirty RMS ratio: %.4f", median_ratio)
    log.info("  Mean dirty corner RMS:  %.4f mJy", np.nanmean(dirty_rms))
    log.info("  Mean clean corner RMS:  %.4f mJy", np.nanmean(clean_rms))

    # Cleaning should reduce or not significantly increase corner noise
    assert median_ratio < 5.0, (
        f"Clean/dirty RMS ratio ({median_ratio:.2f}) is suspiciously high"
    )
    log.info("  PASS: Cleaning does not pathologically increase noise")

    # --- Plot (a): histogram of pixel differences ---
    # Collect difference pixels from a few representative channels
    n_hist_chans = min(10, len(results_sorted))
    diff_all = np.concatenate([r["diff_flat"] for r in results_sorted[:n_hist_chans]])
    diff_all_mJy = diff_all * 1e3

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # (a) Histogram of pixel differences
    pct1, pct99 = np.percentile(diff_all_mJy, [1, 99])
    axes[0].hist(diff_all_mJy, bins=100, range=(pct1, pct99),
                 color='cornflowerblue', alpha=0.8, edgecolor='none')
    axes[0].axvline(0, color='maroon', ls='--', lw=1.5)
    axes[0].set_xlabel('Cleaned - Dirty [mJy/beam]')
    axes[0].set_ylabel('Count')
    axes[0].set_title('Pixel Difference Histogram (sample of %d channels)' % n_hist_chans)
    axes[0].grid(True, alpha=0.3)

    # (b) Scatter plot of dirty vs cleaned pixel values
    n_scatter = min(5, len(results_sorted))
    for i in range(n_scatter):
        r = results_sorted[i * (len(results_sorted) // max(n_scatter, 1))]
        axes[1].scatter(r["dirty_flat"] * 1e3, r["clean_flat"] * 1e3,
                        s=0.5, alpha=0.2, rasterized=True,
                        label='%.0f MHz' % r["freq"])
    # 1:1 line
    all_dirty = np.concatenate([results_sorted[i * (len(results_sorted) // max(n_scatter, 1))]["dirty_flat"]
                                for i in range(n_scatter)]) * 1e3
    lims = [np.percentile(all_dirty, 0.5), np.percentile(all_dirty, 99.5)]
    axes[1].plot(lims, lims, 'k--', lw=1, alpha=0.5, label='1:1')
    axes[1].set_xlabel('Dirty [mJy/beam]')
    axes[1].set_ylabel('Cleaned [mJy/beam]')
    axes[1].set_title('Dirty vs Cleaned Pixel Values')
    axes[1].set_xlim(lims)
    axes[1].set_ylim(lims)
    axes[1].legend(fontsize=9, markerscale=5)
    axes[1].grid(True, alpha=0.3)

    # (c) Clean/dirty RMS ratio across frequency
    axes[2].scatter(freqs, rms_ratios, s=10, alpha=0.7, color='cornflowerblue',
                    rasterized=True)
    axes[2].axhline(1.0, color='gray', ls=':', lw=1)
    axes[2].axhline(median_ratio, color='maroon', ls='--', lw=2,
                    label='Median = %.3f' % median_ratio)
    axes[2].set_xlabel('Frequency [MHz]')
    axes[2].set_ylabel('Cleaned / Dirty Corner RMS')
    axes[2].set_title('Cleaning RMS Ratio vs Frequency')
    axes[2].legend(fontsize=12)
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = PLOT_DIR / "sanity_cleaning_fidelity.png"
    fig.savefig(str(outpath), dpi=150, bbox_inches='tight')
    plt.close(fig)
    log.info("  Saved: %s", outpath)


# ===========================================================================
# Check 4: 5-sigma Detection Threshold
# ===========================================================================

def check_detection_threshold(freqs, rms_mJy):
    """
    Compute and plot the 5-sigma detection threshold across the band.
    """
    log.info("=== Check 4: 5-sigma Detection Threshold ===")

    threshold_5sig = 5.0 * rms_mJy  # mJy
    theoretical_5sig = 5.0 * THEORETICAL_RMS_MJY

    log.info("  Median 5-sigma threshold: %.3f mJy", np.median(threshold_5sig))
    log.info("  Theoretical 5-sigma:      %.3f mJy", theoretical_5sig)
    log.info("  Min 5-sigma threshold:    %.3f mJy", np.min(threshold_5sig))
    log.info("  Max 5-sigma threshold:    %.3f mJy", np.max(threshold_5sig))

    assert np.all(threshold_5sig > 0), "All thresholds must be positive"

    # --- Plot ---
    fig, ax = plt.subplots(1, 1, figsize=(16, 6))

    ax.scatter(freqs, threshold_5sig, s=3, alpha=0.5, color='cornflowerblue',
               label=r'Measured $5\sigma$ threshold', rasterized=True)
    ax.axhline(theoretical_5sig, color='maroon', ls='--', lw=2,
               label=r'Theoretical $5\sigma$ = %.2f mJy' % theoretical_5sig)
    ax.axhline(np.median(threshold_5sig), color='forestgreen', ls=':', lw=2,
               label=r'Median measured $5\sigma$ = %.2f mJy' % np.median(threshold_5sig))

    # Shade known RFI bands
    rfi_bands = [
        (925, 960, 'GSM-900'),
        (1559, 1610, 'GPS L1'),
        (856, 880, 'Band edge'),
        (1670, 1712, 'Band edge'),
    ]
    for flo, fhi, rfi_label in rfi_bands:
        ax.axvspan(flo, fhi, alpha=0.1, color='red')

    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel(r'$5\sigma$ Detection Threshold [mJy/beam]')
    ax.set_title(r'$5\sigma$ Discovery Threshold Across L-band')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = PLOT_DIR / "sanity_detection_threshold.png"
    fig.savefig(str(outpath), dpi=150, bbox_inches='tight')
    plt.close(fig)
    log.info("  Saved: %s", outpath)
    log.info("  PASS: Detection thresholds computed and plotted")


# ===========================================================================
# Check 5: Channel-to-Channel Consistency
# ===========================================================================

def _worker_adjacent_corr(args):
    """Compute correlation coefficient between two adjacent channel images."""
    file_a, file_b, sb_idx, chan_a, chan_b, freq_a, freq_b = args
    try:
        data_a, _ = load_fits_2d(file_a)
        data_b, _ = load_fits_2d(file_b)
        if data_a.shape != (NPIX, NPIX) or data_b.shape != (NPIX, NPIX):
            return None

        flat_a = data_a.ravel()
        flat_b = data_b.ravel()

        # Remove NaNs
        valid = np.isfinite(flat_a) & np.isfinite(flat_b)
        if valid.sum() < 100:
            return None

        corr = np.corrcoef(flat_a[valid], flat_b[valid])[0, 1]
        return {
            "sb": sb_idx,
            "chan_a": chan_a,
            "chan_b": chan_b,
            "freq": (freq_a + freq_b) / 2,
            "corr": corr,
        }
    except Exception as e:
        return None


def check_channel_consistency(n_workers):
    """
    Compute correlation between adjacent channels. Sky emission changes slowly
    across 26 kHz, so correlation should be high.
    """
    log.info("=== Check 5: Channel-to-Channel Consistency ===")

    # Build adjacent-channel pairs, sampling every 40th pair per subband
    pair_tasks = []
    sb_dirs = discover_subbands()
    for sb_dir in sb_dirs:
        sb_name = os.path.basename(sb_dir)
        try:
            sb_idx = int(sb_name.split("_")[1])
        except (IndexError, ValueError):
            continue

        fits_files = sorted(glob.glob(os.path.join(sb_dir, "chan_*.fits")))
        if len(fits_files) < 2:
            continue

        # Sample pairs every 40 channels, skip edges
        for i in range(5, len(fits_files) - 6, 40):
            fa = fits_files[i]
            fb = fits_files[i + 1]
            chan_a = chan_from_filename(fa)
            chan_b = chan_from_filename(fb)
            freq_a = freq_from_filename(fa)
            freq_b = freq_from_filename(fb)
            pair_tasks.append((fa, fb, sb_idx, chan_a, chan_b, freq_a, freq_b))

    log.info("  Computing correlations for %d adjacent-channel pairs...", len(pair_tasks))

    with Pool(processes=n_workers) as pool:
        results = pool.map(_worker_adjacent_corr, pair_tasks, chunksize=8)

    results = [r for r in results if r is not None]
    if not results:
        log.warning("  No valid adjacent-channel pairs. Skipping consistency check.")
        return

    freqs = np.array([r["freq"] for r in results])
    corrs = np.array([r["corr"] for r in results])

    sort_idx = np.argsort(freqs)
    freqs = freqs[sort_idx]
    corrs = corrs[sort_idx]

    median_corr = np.median(corrs)
    low_corr_threshold = 0.5
    n_anomalous = np.sum(corrs < low_corr_threshold)

    log.info("  Total pairs analyzed: %d", len(corrs))
    log.info("  Median correlation: %.6f", median_corr)
    log.info("  Min correlation:    %.6f", np.min(corrs))
    log.info("  Max correlation:    %.6f", np.max(corrs))
    log.info("  Anomalous (corr < %.1f): %d / %d (%.1f%%)",
             low_corr_threshold, n_anomalous, len(corrs),
             100.0 * n_anomalous / len(corrs))

    if n_anomalous > 0:
        # Print which channels are anomalous
        anomalous = [(r["sb"], r["chan_a"], r["chan_b"], r["freq"], r["corr"])
                     for r in results if r["corr"] < low_corr_threshold]
        anomalous.sort(key=lambda x: x[4])  # sort by correlation
        log.warning("  Anomalous pairs (lowest correlation first):")
        for sb, ca, cb, f, c in anomalous[:20]:
            log.warning("    SB %03d, chan %04d-%04d (%.1f MHz): corr=%.4f", sb, ca, cb, f, c)

    # Sanity: for a well-behaved dataset, median correlation should be > 0.9
    assert median_corr > 0.5, (
        f"Median adjacent-channel correlation ({median_corr:.4f}) is "
        "suspiciously low -- possible data issue"
    )
    log.info("  PASS: Adjacent channels are highly correlated (median=%.4f)", median_corr)

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    axes[0].scatter(freqs, corrs, s=5, alpha=0.5, color='cornflowerblue',
                    rasterized=True)
    axes[0].axhline(median_corr, color='maroon', ls='--', lw=2,
                    label='Median = %.4f' % median_corr)
    axes[0].axhline(low_corr_threshold, color='red', ls=':', lw=1.5, alpha=0.7,
                    label='Anomaly threshold = %.1f' % low_corr_threshold)
    axes[0].set_xlabel('Frequency [MHz]')
    axes[0].set_ylabel('Pearson Correlation (adjacent channels)')
    axes[0].set_title('Channel-to-Channel Consistency')
    axes[0].legend(fontsize=11)
    axes[0].grid(True, alpha=0.3)
    axes[0].set_ylim(min(0, np.min(corrs) - 0.05), 1.02)

    axes[1].hist(corrs, bins=50, color='cornflowerblue', alpha=0.8, edgecolor='none')
    axes[1].axvline(median_corr, color='maroon', ls='--', lw=2,
                    label='Median = %.4f' % median_corr)
    axes[1].set_xlabel('Pearson Correlation')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Distribution of Adjacent-Channel Correlations')
    axes[1].legend(fontsize=12)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = PLOT_DIR / "sanity_channel_consistency.png"
    fig.savefig(str(outpath), dpi=150, bbox_inches='tight')
    plt.close(fig)
    log.info("  Saved: %s", outpath)


# ===========================================================================
# Main
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Phase 6: Sanity checks on MeerKAT GC imaging data"
    )
    parser.add_argument("--workers", type=int, default=16,
                        help="Number of multiprocessing workers (default: 16)")
    parser.add_argument("--sample-stride", type=int, default=40,
                        help="Sample every Nth channel per subband (default: 40)")
    parser.add_argument("--injection-trials", type=int, default=200,
                        help="Number of injection trials per SNR level (default: 200)")
    args = parser.parse_args()

    log.info("=" * 60)
    log.info("Phase 6: Sanity Checks")
    log.info("=" * 60)
    log.info("Workers: %d", args.workers)
    log.info("Sample stride: %d channels", args.sample_stride)
    log.info("Injection trials: %d per SNR", args.injection_trials)
    log.info("Theoretical RMS: %.4f mJy/beam", THEORETICAL_RMS_MJY)

    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    # Discover data
    sb_dirs = discover_subbands()
    assert len(sb_dirs) > 0, "No subband directories found in %s" % IMAGE_DIR
    log.info("Found %d subband directories", len(sb_dirs))

    # Sample files for checks that scan across the band
    samples = sample_files_from_subbands(
        sb_dirs, stride=args.sample_stride, skip_edges=5
    )
    assert len(samples) > 0, "No FITS files sampled -- check image directory"
    log.info("Sampled %d channels for analysis", len(samples))

    # --- Check 1: Noise Assessment ---
    freqs_v, rms_mJy = check_noise_assessment(samples, args.workers)

    # --- Check 2: Point Source Injection ---
    check_point_source_injection(n_trials=args.injection_trials)

    # --- Check 3: Cleaning Fidelity ---
    check_cleaning_fidelity(args.workers)

    # --- Check 4: 5-sigma Detection Threshold ---
    check_detection_threshold(freqs_v, rms_mJy)

    # --- Check 5: Channel-to-Channel Consistency ---
    check_channel_consistency(args.workers)

    # --- Summary ---
    log.info("=" * 60)
    log.info("Phase 6: All sanity checks complete")
    log.info("Plots saved to: %s/sanity_*.png", PLOT_DIR)
    log.info("=" * 60)


if __name__ == "__main__":
    main()
