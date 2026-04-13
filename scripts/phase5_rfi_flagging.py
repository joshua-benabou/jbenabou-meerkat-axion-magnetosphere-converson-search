#!/usr/bin/env python3
"""
Phase 5: RFI Channel Flagging for MeerKAT Galactic Center Data
===============================================================

Scans all dirty-image FITS files across 86 subbands, computes per-channel
statistics (RMS in an annular region avoiding Sgr A*, peak, median, kurtosis),
and flags channels contaminated by RFI or instrumental artifacts.

Flagging criteria:
  1. Statistical outlier: RMS > 3x the local sliding-window median (window=21 channels)
  2. Known RFI bands: GSM-900 (925-960 MHz), GPS L1 (1575.42 +/- 2 MHz)
  3. Band-edge artifacts: first and last channel of each subband

Outputs:
  - rfi_channel_flags.csv  (per-channel statistics + flag info)
  - plots/rfi_overview.png (RMS vs frequency + per-subband flag fraction)

THIS SCRIPT IS READ-ONLY ON ALL FITS FILES. It only reads images and writes
new output files (CSV + PNG).
"""

import os
import re
import sys
import glob
import logging
import argparse
from multiprocessing import Pool
from pathlib import Path

import csv
import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scipy.stats import kurtosis as scipy_kurtosis
except ImportError:
    def scipy_kurtosis(a, fisher=True):
        """Fallback kurtosis if scipy unavailable."""
        m = np.mean(a)
        s = np.std(a)
        if s == 0:
            return 0.0
        k = np.mean(((a - m) / s) ** 4)
        return k - 3.0 if fisher else k

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path("/global/scratch/projects/pc_heptheory/jbenabou/"
                    "NS_megaproject/MeerKAT_data/meerkat_reduction_project")
IMAGE_DIR = PROJECT_ROOT / "images"
OUTPUT_CSV = PROJECT_ROOT / "rfi_channel_flags.csv"
OUTPUT_PLOT = PROJECT_ROOT / "plots" / "rfi_overview.png"

# Image geometry (512x512, 2 arcsec/pixel)
NPIX = 512
CENTER = NPIX // 2  # 256

# Annular region for RMS: avoid central source (Sgr A*) and image edges.
# Inner radius: 80 px = 160 arcsec (safely outside Sgr A* + sidelobes)
# Outer radius: 220 px = 440 arcsec (avoids edge effects at 256 px)
ANNULUS_INNER = 80
ANNULUS_OUTER = 220

# RFI flagging parameters
RMS_OUTLIER_FACTOR = 3.0
SLIDING_WINDOW = 21  # channels, must be odd

# Known RFI bands (MHz)
RFI_BANDS = [
    ("GSM-900",  925.0, 960.0),
    ("GPS_L1", 1573.42, 1577.42),
]

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
# Pre-compute the annular mask (shared across workers via fork)
# ---------------------------------------------------------------------------

def _make_annular_mask(npix, r_inner, r_outer):
    """Boolean mask: True inside the annulus."""
    y, x = np.ogrid[:npix, :npix]
    r2 = (x - npix // 2) ** 2 + (y - npix // 2) ** 2
    return (r2 >= r_inner ** 2) & (r2 <= r_outer ** 2)

ANNULAR_MASK = _make_annular_mask(NPIX, ANNULUS_INNER, ANNULUS_OUTER)

# ---------------------------------------------------------------------------
# Worker function
# ---------------------------------------------------------------------------

def compute_channel_stats(args):
    """
    Read one FITS file and return per-channel statistics.

    Parameters
    ----------
    args : tuple
        (subband_idx: int, channel_idx: int, freq_mhz: float, filepath: str)

    Returns
    -------
    dict with keys: subband, channel, frequency_mhz, rms, peak, median, kurtosis
    """
    subband_idx, channel_idx, freq_mhz, filepath = args

    result = {
        "subband": subband_idx,
        "channel": channel_idx,
        "frequency_mhz": freq_mhz,
        "rms": np.nan,
        "peak": np.nan,
        "median": np.nan,
        "kurtosis": np.nan,
    }

    try:
        with fits.open(filepath, memmap=True, mode="readonly") as hdul:
            data = hdul[0].data
            # Handle extra dimensions: could be (1,1,ny,nx) or (ny,nx)
            data = np.squeeze(data).astype(np.float64)

            if data.shape != (NPIX, NPIX):
                log.warning("Unexpected shape %s for %s", data.shape, filepath)
                return result

            annulus_data = data[ANNULAR_MASK]
            # Drop NaNs
            annulus_data = annulus_data[np.isfinite(annulus_data)]

            if annulus_data.size == 0:
                return result

            result["rms"] = float(np.std(annulus_data))
            result["peak"] = float(np.nanmax(np.abs(data)))
            result["median"] = float(np.nanmedian(data))
            result["kurtosis"] = float(scipy_kurtosis(annulus_data, fisher=True))

    except Exception as e:
        log.error("Failed to read %s: %s", filepath, e)

    return result

# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

def discover_fits_files(image_dir):
    """
    Walk the image directory and return a sorted list of
    (subband_idx, channel_idx, freq_mhz, filepath) tuples.
    """
    tasks = []
    subband_dirs = sorted(glob.glob(str(image_dir / "subband_*")))

    for sb_dir in subband_dirs:
        sb_name = os.path.basename(sb_dir)
        # Extract subband index
        try:
            sb_idx = int(sb_name.split("_")[1])
        except (IndexError, ValueError):
            continue

        fits_files = sorted(glob.glob(os.path.join(sb_dir, "chan_*.fits")))
        for fpath in fits_files:
            fname = os.path.basename(fpath)
            m = CHAN_RE.match(fname)
            if m is None:
                continue
            chan_idx = int(m.group(1))
            freq_mhz = float(m.group(2))
            tasks.append((sb_idx, chan_idx, freq_mhz, fpath))

    log.info("Discovered %d FITS files across %d subbands",
             len(tasks), len(subband_dirs))
    return tasks

# ---------------------------------------------------------------------------
# Flagging logic
# ---------------------------------------------------------------------------

def flag_channels(records):
    """
    Add rfi_flag and flag_reason to each record dict.
    records: list of dicts with keys subband, channel, frequency_mhz, rms, ...
    """
    for r in records:
        r["rfi_flag"] = False
        r["flag_reason"] = ""

    # Group by subband
    by_subband = {}
    for i, r in enumerate(records):
        by_subband.setdefault(r["subband"], []).append(i)

    # --- 1. Band-edge channels (first and last per subband) ----------------
    for sb_idx, indices in by_subband.items():
        channels = [records[i]["channel"] for i in indices]
        ch_min, ch_max = min(channels), max(channels)
        for i in indices:
            if records[i]["channel"] in (ch_min, ch_max):
                records[i]["rfi_flag"] = True
                records[i]["flag_reason"] = _append_reason(
                    records[i]["flag_reason"], "band_edge")

    # --- 2. Known RFI bands ------------------------------------------------
    for name, flo, fhi in RFI_BANDS:
        for r in records:
            if flo <= r["frequency_mhz"] <= fhi:
                r["rfi_flag"] = True
                r["flag_reason"] = _append_reason(r["flag_reason"], name)

    # --- 3. Statistical outliers per subband (sliding-window median) -------
    half_w = SLIDING_WINDOW // 2
    for sb_idx, indices in by_subband.items():
        # Sort indices by channel
        indices_sorted = sorted(indices, key=lambda i: records[i]["channel"])
        rms_arr = np.array([records[i]["rms"] for i in indices_sorted])

        n = len(rms_arr)
        local_median = np.full(n, np.nan)
        for j in range(n):
            lo = max(0, j - half_w)
            hi = min(n, j + half_w + 1)
            window = rms_arr[lo:hi]
            window = window[np.isfinite(window)]
            if window.size > 0:
                local_median[j] = np.median(window)

        for j in range(n):
            if (np.isfinite(rms_arr[j]) and np.isfinite(local_median[j])
                    and rms_arr[j] > RMS_OUTLIER_FACTOR * local_median[j]):
                idx = indices_sorted[j]
                records[idx]["rfi_flag"] = True
                records[idx]["flag_reason"] = _append_reason(
                    records[idx]["flag_reason"], "rms_outlier")

    return records


def _append_reason(existing, new):
    """Append a flag reason string, semicolon-separated."""
    if existing:
        return existing + "; " + new
    return new

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def make_overview_plot(records, output_path):
    """
    Two-panel plot:
      Top:    RMS vs frequency, flagged channels in red
      Bottom: fraction of flagged channels per subband
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    freqs = np.array([r["frequency_mhz"] for r in records])
    rms = np.array([r["rms"] for r in records])
    flags = np.array([r["rfi_flag"] for r in records])
    subbands = np.array([r["subband"] for r in records])

    good = ~flags
    bad = flags

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 8),
                                    gridspec_kw={"height_ratios": [3, 1]},
                                    sharex=False)

    # --- Top panel: RMS vs frequency ---
    ax1.scatter(freqs[good], rms[good],
                s=0.3, alpha=0.4, color="steelblue", label="Good", rasterized=True)
    ax1.scatter(freqs[bad], rms[bad],
                s=1.0, alpha=0.6, color="red", label="Flagged", rasterized=True)

    for name, flo, fhi in RFI_BANDS:
        ax1.axvspan(flo, fhi, color="salmon", alpha=0.15, label=name)

    ax1.set_ylabel("Annular RMS (Jy/beam)")
    ax1.set_title("RFI Flagging Overview: RMS vs Frequency")
    ax1.set_yscale("log")
    handles, labels = ax1.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax1.legend(by_label.values(), by_label.keys(), loc="upper right", fontsize=8)
    ax1.grid(True, alpha=0.3)

    # --- Bottom panel: flagged fraction per subband ---
    sb_unique = sorted(set(subbands))
    frac_list = []
    freq_min_list = []
    for sb in sb_unique:
        mask = subbands == sb
        n_total = mask.sum()
        n_flag = (flags & mask).sum()
        frac_list.append(n_flag / n_total if n_total > 0 else 0)
        freq_min_list.append(freqs[mask].min())

    ax2.bar(sb_unique, frac_list, color="indianred", edgecolor="none")
    ax2.set_xlabel("Subband index")
    ax2.set_ylabel("Flagged fraction")
    ax2.set_title("Fraction of Flagged Channels per Subband")
    ax2.set_ylim(0, min(1.05, max(frac_list) * 1.3 + 0.05) if frac_list else 1.0)
    ax2.grid(True, alpha=0.3, axis="y")

    ax2b = ax2.twiny()
    ax2b.set_xlim(ax2.get_xlim())
    n_ticks = 10
    tick_positions = np.linspace(0, 85, n_ticks, dtype=int)
    sb_freq_map = dict(zip(sb_unique, freq_min_list))
    tick_labels = [f"{sb_freq_map[t]:.0f}" if t in sb_freq_map else "" for t in tick_positions]
    ax2b.set_xticks(tick_positions)
    ax2b.set_xticklabels(tick_labels, fontsize=7)
    ax2b.set_xlabel("Approx. frequency (MHz)", fontsize=8)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved overview plot to %s", output_path)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Phase 5: RFI channel flagging")
    parser.add_argument("--workers", type=int, default=16,
                        help="Number of multiprocessing workers (default: 16)")
    parser.add_argument("--image-dir", type=str, default=str(IMAGE_DIR),
                        help="Path to images/ directory")
    parser.add_argument("--output-csv", type=str, default=str(OUTPUT_CSV),
                        help="Output CSV path")
    parser.add_argument("--output-plot", type=str, default=str(OUTPUT_PLOT),
                        help="Output plot path")
    args = parser.parse_args()

    image_dir = Path(args.image_dir)
    output_csv = Path(args.output_csv)
    output_plot = Path(args.output_plot)

    log.info("Phase 5: RFI channel flagging")
    log.info("Image directory: %s", image_dir)
    log.info("Workers: %d", args.workers)

    # --- Discover files ---
    tasks = discover_fits_files(image_dir)
    if not tasks:
        log.error("No FITS files found in %s", image_dir)
        sys.exit(1)

    # --- Compute statistics in parallel ---
    log.info("Computing per-channel statistics with %d workers ...", args.workers)
    with Pool(processes=args.workers) as pool:
        results = pool.map(compute_channel_stats, tasks, chunksize=64)

    log.info("Statistics computed for %d channels", len(results))

    # Sort by subband, channel
    results.sort(key=lambda r: (r["subband"], r["channel"]))

    # --- Flag channels ---
    log.info("Applying flagging criteria ...")
    results = flag_channels(results)

    n_flagged = sum(1 for r in results if r["rfi_flag"])
    n_total = len(results)
    log.info("Flagged %d / %d channels (%.1f%%)",
             n_flagged, n_total, 100.0 * n_flagged / n_total)

    # --- Write CSV ---
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["subband", "channel", "frequency_mhz", "rms", "peak",
                  "median", "kurtosis", "rfi_flag", "flag_reason"]
    with open(str(output_csv), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    log.info("Saved channel flags to %s", output_csv)

    # --- Summary per flag reason ---
    if n_flagged > 0:
        reason_counts = {}
        for r in results:
            if r["rfi_flag"] and r["flag_reason"]:
                for reason in r["flag_reason"].split("; "):
                    reason_counts[reason] = reason_counts.get(reason, 0) + 1
        log.info("Flag breakdown:")
        for reason, count in sorted(reason_counts.items(), key=lambda x: -x[1]):
            log.info("  %-20s %5d channels", reason, count)

    # --- Plot ---
    log.info("Generating overview plot ...")
    make_overview_plot(results, output_plot)

    log.info("Phase 5 complete.")


if __name__ == "__main__":
    main()
