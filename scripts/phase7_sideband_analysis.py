#!/usr/bin/env python3
"""
Phase 7: Sideband Analysis for Axion Search

Connects the signal templates (from phase7_signal_templates.py) to the
actual MeerKAT FITS channel images and performs sideband background
subtraction to search for excess flux consistent with axion-photon
conversion.

Strategy:
    For each channel with predicted axion signal, we:
    1. Load the FITS image for that channel
    2. Load images from neighboring "sideband" channels (excluding a
       guard band around the signal channel)
    3. Build a pixel-by-pixel background model from the sideband images
    4. Subtract the background and measure the residual flux
    5. Compare to the expected signal from the template

    The sideband naturally removes continuum emission and any smooth
    spectral structure, isolating frequency-dependent features like
    a narrow axion line.

This script uses functions from phase7_axion_search.py for the
image loading / sideband machinery, and from phase7_signal_templates.py
for the signal template loading.

Usage:
    python phase7_sideband_analysis.py --mass 4.13e-6 --coupling 1e-11
    python phase7_sideband_analysis.py --mass 4.13e-6 --use_dirty --n_sideband 30
    python phase7_sideband_analysis.py --mass 4.13e-6 --top_channels 20

Output:
    - Results CSV: results/phase7_sideband_analysis/<mass>_results.csv
    - Diagnostic plots: plots/phase7_sideband_*.png
"""

import os
os.environ["PATH"] = (
    '/global/software/sl-7.x86_64/modules/tools/texlive/2016/bin/x86_64-linux/:'
    + os.environ["PATH"]
)

import sys
sys.path.insert(0, '/global/scratch/projects/pc_heptheory/jbenabou')
import plotting_defaults

# Add scripts dir to path so we can import sibling modules
SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPTS_DIR)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import csv
import glob
import warnings
from typing import List, Tuple, Optional, Dict

# Import from our other Phase 7 modules
from phase7_axion_search import (
    IMAGES_DIR, N_SUBBANDS, CHANS_PER_SUBBAND, TOTAL_CHANS,
    CHAN_WIDTH_KHZ, FIRST_FREQ_MHZ, IMSIZE, PIXEL_SCALE_ARCSEC,
    freq_to_global_channel, global_channel_to_subband,
    channel_to_freq_mhz, find_channel_fits,
    load_fits_image, save_fits_image,
    get_sideband_channels, build_background_model, sideband_subtract,
    measure_local_rms, measure_peak_flux,
    load_rfi_flags,
)
from phase7_signal_templates import (
    SignalTemplate, build_template, build_averaged_template,
    load_combined_flux, find_combined_flux_file,
    AVAILABLE_MASSES, POP_DIRS, N_REALIZATIONS,
    H_EV_S, G_REF, CHAN_WIDTH_HZ,
    channel_to_freq_hz as template_channel_to_freq_hz,
)


# =============================================================================
# Paths
# =============================================================================

PROJECT_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project'
)
RESULTS_DIR = os.path.join(PROJECT_DIR, 'results', 'phase7_sideband_analysis')
PLOT_DIR = os.path.join(PROJECT_DIR, 'plots')


# =============================================================================
# Channel-level analysis
# =============================================================================

def analyze_channel(global_chan: int, template: SignalTemplate,
                    rfi_flags: dict,
                    n_sideband: int = 50, n_guard: int = 5,
                    use_cleaned: bool = True,
                    bg_method: str = 'median') -> Optional[dict]:
    """
    Perform sideband analysis on a single channel.

    Parameters
    ----------
    global_chan : int
        Global channel index (0-32767).
    template : SignalTemplate
        Signal template (used to get expected flux for this channel).
    rfi_flags : dict
        RFI flag map.
    n_sideband : int
        Number of sideband channels on each side.
    n_guard : int
        Guard band half-width in channels.
    use_cleaned : bool
        Use cleaned images (True) or dirty (False).
    bg_method : str
        'median' or 'mean'.

    Returns
    -------
    dict or None
        Results dictionary, or None if channel could not be analyzed.
    """
    subband_idx, local_chan = global_channel_to_subband(global_chan)
    freq_mhz = channel_to_freq_mhz(subband_idx, local_chan)

    # Find expected signal in this channel
    chan_mask = template.channel_indices == global_chan
    if np.any(chan_mask):
        expected_flux_jy = float(template.flux_density_jy[chan_mask][0])
    else:
        expected_flux_jy = 0.0

    # Load signal channel image
    fits_path = find_channel_fits(subband_idx, local_chan, cleaned=use_cleaned)
    if fits_path is None:
        return None

    try:
        signal_image, header = load_fits_image(fits_path)
    except Exception as e:
        print(f"  WARNING: Failed to load {fits_path}: {e}")
        return None

    # Build sideband channel list
    sideband_channels = get_sideband_channels(
        subband_idx, local_chan, n_sideband, n_guard, rfi_flags
    )

    if len(sideband_channels) < 5:
        print(f"  WARNING: Only {len(sideband_channels)} sideband channels for "
              f"chan {global_chan}. Skipping.")
        return None

    # Build background model
    background, n_loaded = build_background_model(
        sideband_channels, cleaned=use_cleaned, method=bg_method
    )
    if background is None:
        return None

    # Subtract
    residual = sideband_subtract(signal_image, background)

    # Measure statistics at image center (Sgr A* / GC position)
    # The GC is at the image center for our pointing
    x_center = IMSIZE // 2
    y_center = IMSIZE // 2

    # Measure peak flux in a central aperture
    aperture_pix = 20  # ~40 arcsec diameter
    peak_flux = measure_peak_flux(residual, x_center, y_center,
                                  aperture_radius_pix=aperture_pix)
    local_rms = measure_local_rms(residual, x_center, y_center,
                                  exclude_radius_pix=aperture_pix + 5)

    # Also measure the total flux in the central aperture
    ny, nx = residual.shape
    yy, xx = np.ogrid[:ny, :nx]
    r2 = (xx - x_center) ** 2 + (yy - y_center) ** 2
    aperture_mask = r2 <= aperture_pix ** 2
    aperture_pixels = residual[aperture_mask]
    aperture_pixels = aperture_pixels[np.isfinite(aperture_pixels)]
    total_aperture_flux = float(np.sum(aperture_pixels)) if len(aperture_pixels) > 0 else np.nan
    mean_aperture_flux = float(np.mean(aperture_pixels)) if len(aperture_pixels) > 0 else np.nan

    # Significance
    significance = peak_flux / local_rms if (np.isfinite(local_rms) and local_rms > 0) else 0.0

    # Global image statistics
    valid = residual[np.isfinite(residual)]
    global_rms = float(np.std(valid)) if len(valid) > 0 else np.nan
    global_mean = float(np.mean(valid)) if len(valid) > 0 else np.nan

    return {
        'global_channel': global_chan,
        'subband': subband_idx,
        'local_channel': local_chan,
        'freq_mhz': freq_mhz,
        'expected_flux_jy': expected_flux_jy,
        'expected_flux_ujy': expected_flux_jy * 1e6,
        'peak_residual_jy': peak_flux,
        'mean_aperture_jy': mean_aperture_flux,
        'total_aperture_jy': total_aperture_flux,
        'local_rms_jy': local_rms,
        'global_rms_jy': global_rms,
        'significance_sigma': significance,
        'n_sideband_loaded': n_loaded,
        'n_sideband_requested': len(sideband_channels),
        'fits_path': fits_path,
    }


# =============================================================================
# Full analysis pipeline
# =============================================================================

def run_sideband_analysis(template: SignalTemplate,
                          rfi_flags: dict,
                          n_sideband: int = 50,
                          n_guard: int = 5,
                          use_cleaned: bool = True,
                          bg_method: str = 'median',
                          top_n: Optional[int] = None,
                          ) -> List[dict]:
    """
    Run sideband analysis for all channels with predicted signal.

    Parameters
    ----------
    template : SignalTemplate
        The signal template.
    rfi_flags : dict
        RFI flag map.
    n_sideband, n_guard : int
        Sideband parameters.
    use_cleaned : bool
        Use cleaned or dirty images.
    bg_method : str
        Background estimation method.
    top_n : int, optional
        If set, only analyze the top N channels by expected flux.

    Returns
    -------
    list of dict
        Results for each analyzed channel.
    """
    print("\n" + "=" * 70)
    print("SIDEBAND ANALYSIS")
    print("=" * 70)
    print(f"  Mass:         {template.mass_ev:.4e} eV")
    print(f"  Coupling:     {template.coupling_gev:.2e} GeV^-1")
    print(f"  Population:   {template.pop_model} (Pop {template.pop_idx})")
    print(f"  n_sideband:   {n_sideband}")
    print(f"  n_guard:      {n_guard}")
    print(f"  bg_method:    {bg_method}")
    print(f"  Image type:   {'cleaned' if use_cleaned else 'dirty'}")
    print(f"  Total active channels: {len(template.channel_indices)}")

    # Select channels to analyze
    channels_to_analyze = template.channel_indices.copy()
    fluxes = template.flux_density_jy.copy()

    if top_n is not None and top_n < len(channels_to_analyze):
        # Sort by flux and take top N
        sort_idx = np.argsort(fluxes)[::-1]
        channels_to_analyze = channels_to_analyze[sort_idx[:top_n]]
        fluxes = fluxes[sort_idx[:top_n]]
        print(f"  Analyzing top {top_n} channels by expected flux")

    print(f"  Channels to analyze: {len(channels_to_analyze)}")

    results = []
    n_success = 0
    n_fail = 0

    for i, global_chan in enumerate(channels_to_analyze):
        if (i + 1) % 50 == 0 or i == 0:
            print(f"\n  Processing channel {i+1}/{len(channels_to_analyze)} "
                  f"(global chan {global_chan})...", flush=True)

        result = analyze_channel(
            global_chan, template, rfi_flags,
            n_sideband=n_sideband, n_guard=n_guard,
            use_cleaned=use_cleaned, bg_method=bg_method,
        )

        if result is not None:
            results.append(result)
            n_success += 1
        else:
            n_fail += 1

    print(f"\n  Analysis complete: {n_success} channels analyzed, {n_fail} failed/skipped")
    return results


def save_results(results: List[dict], template: SignalTemplate,
                 outdir: str = RESULTS_DIR) -> str:
    """Save results to CSV."""
    os.makedirs(outdir, exist_ok=True)
    mass_label = f"{template.mass_ev:.3e}eV"
    pop_label = f"{template.pop_model}_pop{template.pop_idx}"
    fname = f"sideband_results_{mass_label}_{pop_label}.csv"
    outpath = os.path.join(outdir, fname)

    if not results:
        print("  No results to save.")
        return ""

    fieldnames = list(results[0].keys())
    with open(outpath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    print(f"  Results saved to: {outpath}")
    return outpath


# =============================================================================
# Diagnostic plots
# =============================================================================

def plot_sideband_diagnostics(results: List[dict], template: SignalTemplate,
                              outdir: str = PLOT_DIR) -> List[str]:
    """
    Generate diagnostic plots from the sideband analysis results.

    Returns list of output paths.
    """
    os.makedirs(outdir, exist_ok=True)
    if not results:
        print("  No results to plot.")
        return []

    outpaths = []
    mass_label = f"{template.mass_ev*1e6:.2f}"
    pop_label = f"{template.pop_model}"
    pidx = template.pop_idx if template.pop_idx >= 0 else 'avg'

    freqs = np.array([r['freq_mhz'] for r in results])
    expected = np.array([r['expected_flux_ujy'] for r in results])
    peak_res = np.array([r['peak_residual_jy'] for r in results]) * 1e6  # to uJy
    rms = np.array([r['local_rms_jy'] for r in results]) * 1e6
    significance = np.array([r['significance_sigma'] for r in results])
    mean_ap = np.array([r['mean_aperture_jy'] for r in results]) * 1e6

    # --- Figure 1: Expected vs Observed ---
    fig, axes = plt.subplots(2, 2, figsize=(18, 12))

    # Panel 1: Expected signal template
    ax = axes[0, 0]
    ax.plot(freqs, expected, 'o', ms=3, color='cornflowerblue', alpha=0.6)
    ax.set_xlabel(r'Frequency [MHz]')
    ax.set_ylabel(r'Expected Flux Density [$\mu$Jy]')
    ax.set_title(r'Signal Template')
    ax.set_yscale('log')

    # Panel 2: Measured peak residual
    ax = axes[0, 1]
    ax.plot(freqs, peak_res, 'o', ms=3, color='forestgreen', alpha=0.6)
    ax.axhline(0, color='gray', ls='--', alpha=0.5)
    ax.set_xlabel(r'Frequency [MHz]')
    ax.set_ylabel(r'Peak Residual [$\mu$Jy/beam]')
    ax.set_title(r'Sideband-Subtracted Peak Residual')

    # Panel 3: Local RMS
    ax = axes[1, 0]
    ax.plot(freqs, rms, 'o', ms=3, color='maroon', alpha=0.6)
    ax.set_xlabel(r'Frequency [MHz]')
    ax.set_ylabel(r'Local RMS [$\mu$Jy/beam]')
    ax.set_title(r'Noise Level per Channel')

    # Panel 4: Significance
    ax = axes[1, 1]
    ax.plot(freqs, significance, 'o', ms=3, color='goldenrod', alpha=0.6)
    ax.axhline(5.0, color='maroon', ls='--', alpha=0.5, label=r'$5\sigma$')
    ax.axhline(3.0, color='goldenrod', ls='--', alpha=0.5, label=r'$3\sigma$')
    ax.set_xlabel(r'Frequency [MHz]')
    ax.set_ylabel(r'Significance [$\sigma$]')
    ax.set_title(r'Detection Significance')
    ax.legend(fontsize=12)

    plt.suptitle(
        r'Phase 7 Sideband Analysis: $m_a = %s\;\mu$eV, %s Pop %s' %
        (mass_label, pop_label.capitalize(), pidx),
        fontsize=16
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    outpath = os.path.join(
        outdir,
        f'phase7_sideband_diagnostics_ma{template.mass_ev:.3e}eV_{pop_label}_pop{pidx}.png'
    )
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Plot saved: {outpath}")
    outpaths.append(outpath)

    # --- Figure 2: Sensitivity comparison ---
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Panel 1: Expected signal vs noise level
    ax = axes[0]
    ax.plot(freqs, expected, 'o', ms=3, color='cornflowerblue', alpha=0.6,
            label=r'Expected signal ($g = %.0e$ GeV$^{-1}$)' % template.coupling_gev)
    ax.plot(freqs, rms, 'o', ms=3, color='maroon', alpha=0.4,
            label=r'Local RMS (noise)')
    ax.set_xlabel(r'Frequency [MHz]')
    ax.set_ylabel(r'Flux Density [$\mu$Jy]')
    ax.set_title(r'Signal vs Noise')
    ax.set_yscale('log')
    ax.legend(fontsize=12)

    # Panel 2: Required coupling for 5-sigma detection
    # Need g such that expected_flux(g) = 5 * rms
    # Since flux ~ g^4: g_detect = g_ref * (5*rms / expected)^(1/4)
    ax = axes[1]
    valid_expected = expected > 0
    if np.any(valid_expected):
        g_detect = template.coupling_gev * (5.0 * rms[valid_expected] /
                                             expected[valid_expected]) ** 0.25
        ax.plot(freqs[valid_expected], g_detect * 1e12, 'o', ms=3,
                color='forestgreen', alpha=0.6)
        ax.set_xlabel(r'Frequency [MHz]')
        ax.set_ylabel(r'$g_{a\gamma\gamma}$ [$10^{-12}$ GeV$^{-1}$]')
        ax.set_title(r'Estimated $5\sigma$ Sensitivity per Channel')
        ax.set_yscale('log')

        # Print best sensitivity
        best_idx = np.argmin(g_detect)
        best_g = g_detect[best_idx]
        best_freq = freqs[valid_expected][best_idx]
        print(f"\n  Best single-channel sensitivity:")
        print(f"    g_5sigma = {best_g:.3e} GeV^-1 at {best_freq:.3f} MHz")

    plt.tight_layout()
    outpath = os.path.join(
        outdir,
        f'phase7_sideband_sensitivity_ma{template.mass_ev:.3e}eV_{pop_label}_pop{pidx}.png'
    )
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Plot saved: {outpath}")
    outpaths.append(outpath)

    # --- Figure 3: Significance histogram ---
    fig, ax = plt.subplots(figsize=(8, 6))
    finite_sig = significance[np.isfinite(significance)]
    if len(finite_sig) > 0:
        ax.hist(finite_sig, bins=50, color='cornflowerblue', alpha=0.8,
                edgecolor='navy', density=True)
        # Overlay expected Gaussian
        x_gauss = np.linspace(-5, 5, 200)
        ax.plot(x_gauss, np.exp(-x_gauss**2 / 2) / np.sqrt(2 * np.pi),
                color='maroon', lw=2, ls='--', label=r'Gaussian ($\mu=0, \sigma=1$)')
        ax.set_xlabel(r'Significance [$\sigma$]')
        ax.set_ylabel(r'Probability Density')
        ax.set_title(r'Distribution of Sideband Residual Significance')
        ax.legend(fontsize=12)

        # Print summary
        print(f"\n  Significance distribution:")
        print(f"    Mean:   {np.mean(finite_sig):.3f}")
        print(f"    Std:    {np.std(finite_sig):.3f}")
        print(f"    Max:    {np.max(finite_sig):.3f}")
        print(f"    >3sig:  {np.sum(finite_sig > 3)}")
        print(f"    >5sig:  {np.sum(finite_sig > 5)}")

    plt.tight_layout()
    outpath = os.path.join(
        outdir,
        f'phase7_sideband_significance_hist_ma{template.mass_ev:.3e}eV_{pop_label}_pop{pidx}.png'
    )
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Plot saved: {outpath}")
    outpaths.append(outpath)

    return outpaths


# =============================================================================
# Upper limit estimation
# =============================================================================

def estimate_upper_limit(results: List[dict], template: SignalTemplate,
                         confidence: float = 0.95) -> dict:
    """
    Estimate a simple frequentist upper limit on the coupling constant.

    For each channel, the upper limit on flux is roughly:
        S_upper = peak_residual + z_alpha * local_rms

    where z_alpha is the one-sided quantile (e.g., 1.645 for 95%).

    The coupling upper limit from channel i is:
        g_upper_i = g_ref * (S_upper_i / S_expected_i)^(1/4)

    The combined limit is the minimum over all channels (most constraining).

    Parameters
    ----------
    results : list of dict
    template : SignalTemplate
    confidence : float
        Confidence level (default 0.95).

    Returns
    -------
    dict with upper limit info.
    """
    from scipy import stats
    z_alpha = stats.norm.ppf(confidence)

    g_upper_per_chan = []
    freq_per_chan = []

    for r in results:
        expected = r['expected_flux_jy']
        if expected <= 0 or not np.isfinite(expected):
            continue
        peak = r['peak_residual_jy']
        rms = r['local_rms_jy']
        if not np.isfinite(peak) or not np.isfinite(rms) or rms <= 0:
            continue

        s_upper = max(peak, 0) + z_alpha * rms
        g_upper = template.coupling_gev * (s_upper / expected) ** 0.25
        g_upper_per_chan.append(g_upper)
        freq_per_chan.append(r['freq_mhz'])

    if not g_upper_per_chan:
        return {'g_upper': np.nan, 'freq_best_mhz': np.nan, 'n_channels': 0}

    g_upper_arr = np.array(g_upper_per_chan)
    freq_arr = np.array(freq_per_chan)

    best_idx = np.argmin(g_upper_arr)

    return {
        'g_upper': g_upper_arr[best_idx],
        'freq_best_mhz': freq_arr[best_idx],
        'n_channels': len(g_upper_arr),
        'g_median': np.median(g_upper_arr),
        'confidence': confidence,
    }


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Phase 7: Sideband analysis for axion search"
    )
    parser.add_argument('--mass', type=float, required=True,
                        help='Axion mass in eV (e.g., 4.13e-6)')
    parser.add_argument('--coupling', type=float, default=1e-12,
                        help='Coupling constant in GeV^-1 (default: 1e-12)')
    parser.add_argument('--pop_model', choices=['young', 'old'], default='young',
                        help='Population model (default: young)')
    parser.add_argument('--pop_idx', type=int, default=0,
                        help='Population realization index (default: 0)')
    parser.add_argument('--all_realizations', action='store_true',
                        help='Average over all realizations for the template')
    parser.add_argument('--n_sideband', type=int, default=50,
                        help='Sideband channels on each side (default: 50)')
    parser.add_argument('--n_guard', type=int, default=5,
                        help='Guard band channels on each side (default: 5)')
    parser.add_argument('--bg_method', choices=['median', 'mean'], default='median',
                        help='Background estimation method (default: median)')
    parser.add_argument('--use_dirty', action='store_true',
                        help='Use dirty images instead of cleaned')
    parser.add_argument('--top_channels', type=int, default=None,
                        help='Only analyze top N channels by expected flux')
    parser.add_argument('--no_plot', action='store_true',
                        help='Skip diagnostic plots')
    args = parser.parse_args()

    print("=" * 70)
    print("Phase 7: Sideband Analysis for Axion Search")
    print("=" * 70)

    # --- Validate mass ---
    mass_ev = args.mass
    closest_mass = min(AVAILABLE_MASSES.keys(), key=lambda m: abs(m - mass_ev))
    if abs(closest_mass - mass_ev) / closest_mass > 0.01:
        print(f"WARNING: Requested mass {mass_ev:.3e} not in pre-computed set.")
        print(f"  Available: {sorted(AVAILABLE_MASSES.keys())}")
        print(f"  Using closest: {closest_mass:.3e}")
    mass_ev = closest_mass

    # --- Build signal template ---
    print(f"\nBuilding signal template for mass = {mass_ev:.3e} eV...")
    if args.all_realizations:
        template, _ = build_averaged_template(
            mass_ev, args.pop_model, args.coupling
        )
    else:
        filepath = find_combined_flux_file(mass_ev, args.pop_model, args.pop_idx)
        data = load_combined_flux(filepath)
        template = build_template(data, mass_ev, args.coupling,
                                  args.pop_model, args.pop_idx)

    # --- Print template summary ---
    nu_0_mhz = mass_ev / H_EV_S / 1e6
    print(f"\n  Axion rest frequency: {nu_0_mhz:.3f} MHz")
    print(f"  Active channels:     {len(template.channel_indices)}")
    if len(template.flux_density_jy) > 0:
        print(f"  Peak expected flux:  {np.max(template.flux_density_jy)*1e6:.4f} uJy")
        print(f"  Total expected flux: {template.total_flux_jy_hz:.4e} Jy-Hz")

    # --- Check what images are available ---
    print(f"\n  Checking available FITS images...")
    use_cleaned = not args.use_dirty
    n_available = 0
    for ch in template.channel_indices:
        sb, lc = global_channel_to_subband(ch)
        if find_channel_fits(sb, lc, cleaned=use_cleaned) is not None:
            n_available += 1
    print(f"  Images found: {n_available} / {len(template.channel_indices)} "
          f"({'cleaned' if use_cleaned else 'dirty'})")

    if n_available == 0:
        print("\n  ERROR: No FITS images found for any signal channels.")
        print("  The imaging pipeline (Phases 2-3) must complete before "
              "running the sideband analysis.")
        print("  You can test the template generation with --no_plot or "
              "check image availability.")

        # Still save the template for reference
        from phase7_signal_templates import save_template, print_template_diagnostics
        print_template_diagnostics(template)
        save_template(template)

        sys.exit(1)

    # --- Load RFI flags ---
    rfi_flags = load_rfi_flags()

    # --- Run sideband analysis ---
    results = run_sideband_analysis(
        template, rfi_flags,
        n_sideband=args.n_sideband,
        n_guard=args.n_guard,
        use_cleaned=use_cleaned,
        bg_method=args.bg_method,
        top_n=args.top_channels,
    )

    # --- Save results ---
    save_results(results, template)

    # --- Print summary ---
    if results:
        print("\n" + "=" * 70)
        print("RESULTS SUMMARY")
        print("=" * 70)

        freqs = [r['freq_mhz'] for r in results]
        sigs = [r['significance_sigma'] for r in results]
        rms_vals = [r['local_rms_jy'] for r in results]

        print(f"  Channels analyzed:     {len(results)}")
        print(f"  Frequency range:       {min(freqs):.3f} -- {max(freqs):.3f} MHz")
        print(f"  Median local RMS:      {np.median(rms_vals)*1e6:.4f} uJy/beam")
        print(f"  Max significance:      {max(sigs):.2f} sigma")
        print(f"  Channels > 3 sigma:    {sum(1 for s in sigs if s > 3)}")
        print(f"  Channels > 5 sigma:    {sum(1 for s in sigs if s > 5)}")

        # Estimate upper limit
        try:
            ul = estimate_upper_limit(results, template)
            if np.isfinite(ul['g_upper']):
                print(f"\n  95% CL upper limit on coupling:")
                print(f"    g_upper = {ul['g_upper']:.3e} GeV^-1")
                print(f"    Best channel: {ul['freq_best_mhz']:.3f} MHz")
                print(f"    Median per-channel: {ul['g_median']:.3e} GeV^-1")
        except ImportError:
            print("\n  (scipy not available -- skipping upper limit estimate)")

    # --- Plots ---
    if not args.no_plot and results:
        plot_sideband_diagnostics(results, template)

    print(f"\nPhase 7 sideband analysis complete.")


if __name__ == '__main__':
    main()
