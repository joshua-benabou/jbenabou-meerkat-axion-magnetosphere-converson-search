#!/usr/bin/env python3
"""
Phase 7: Signal Template Generation from Sam Witte's NS Population Models

Loads pre-computed Combined_Flux.dat files from Sam's ray-tracing code,
bins photon energies into MeerKAT channel frequencies, and produces
per-channel expected flux density templates for a given axion mass and
coupling constant.

Combined_Flux.dat columns:
    [0] NS_index
    [1] flux_weight (Jy-Hz) at g_agamma = 1e-12 GeV^-1
    [2] photon_energy (eV)
    [3] conversion_probability at g = 1e-12 GeV^-1
    [4] x (kpc)
    [5] y (kpc)
    [6] z (kpc)

Key physics:
    - Photon frequency: nu = E_photon / h
    - Flux density: S_nu = sum(flux_weight) / channel_width  [Jy]
    - Coupling scaling: flux ~ g^4, so S(g) = S(g_ref) * (g/g_ref)^4
    - Reference coupling: g_ref = 1e-12 GeV^-1

Usage:
    python phase7_signal_templates.py --mass 4.13e-6
    python phase7_signal_templates.py --mass 3.54e-6 --coupling 1e-11 --pop_model young --pop_idx 0
    python phase7_signal_templates.py --mass 4.13e-6 --all_realizations

Output:
    - Signal template CSV: results/phase7_signal_templates/<mass>_<pop>_template.csv
    - Diagnostic plots: plots/phase7_template_*.png
"""

import os
os.environ["PATH"] = (
    '/global/software/sl-7.x86_64/modules/tools/texlive/2016/bin/x86_64-linux/:'
    + os.environ["PATH"]
)

import sys
sys.path.insert(0, '/global/scratch/projects/pc_heptheory/jbenabou')
import plotting_defaults

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import glob
import csv
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict

# =============================================================================
# Constants and paths
# =============================================================================

PROJECT_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project'
)
SAM_DATA_DIR = os.path.join(
    PROJECT_DIR,
    'reference/Real_Analysis/Real_Analysis/Output_Files'
)
TEMPLATE_DIR = os.path.join(PROJECT_DIR, 'results', 'phase7_signal_templates')
PLOT_DIR = os.path.join(PROJECT_DIR, 'plots')

# Physical constants
H_EV_S = 4.135667696e-15   # Planck constant in eV*s
G_REF = 1e-12               # Reference coupling in GeV^-1

# MeerKAT channel grid
FIRST_FREQ_HZ = 856.0e6     # First channel center frequency
CHAN_WIDTH_HZ = 26123.0      # Channel width in Hz
TOTAL_CHANS = 32768          # Total number of channels
LAST_FREQ_HZ = FIRST_FREQ_HZ + (TOTAL_CHANS - 1) * CHAN_WIDTH_HZ

# Population model directory names
POP_DIRS = {
    'young': 'PopYoung_WB_TauO_1.00e+07_B_12.79_sB_0.51_P_1.07_sP_14.56_',
    'old':   'PopOld_WB_TauO_1.00e+12_B_12.82_sB_0.56_P_0.19_sP_22.19_',
}

# Available masses (eV) and their directory labels
AVAILABLE_MASSES = {
    3.54e-6: 'Ma_3.540e-06',
    4.13e-6: 'Ma_4.130e-06',
    7.08e-6: 'Ma_7.080e-06',
}

# Number of population realizations available
N_REALIZATIONS = {
    'young': 10,  # Pop_0 through Pop_9
    'old':    4,  # Pop_0 through Pop_3
}


# =============================================================================
# Data loading
# =============================================================================

def find_combined_flux_file(mass_ev: float, pop_model: str, pop_idx: int) -> str:
    """
    Find the Combined_Flux.dat file for a given mass, population model, and
    realization index.

    Parameters
    ----------
    mass_ev : float
        Axion mass in eV.
    pop_model : str
        'young' or 'old'.
    pop_idx : int
        Population realization index.

    Returns
    -------
    str
        Full path to Combined_Flux.dat.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    """
    pop_dir = POP_DIRS[pop_model]
    mass_dir = AVAILABLE_MASSES[mass_ev]
    path = os.path.join(SAM_DATA_DIR, pop_dir, mass_dir,
                        f'Pop_{pop_idx}', 'Combined_Flux.dat')
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Combined_Flux.dat not found at {path}\n"
            f"Available masses for {pop_model}: check {os.path.join(SAM_DATA_DIR, pop_dir)}"
        )
    return path


def load_combined_flux(filepath: str) -> np.ndarray:
    """
    Load a Combined_Flux.dat file.

    Returns
    -------
    ndarray of shape (N, 7)
        Columns: NS_index, flux_weight(Jy-Hz), photon_energy(eV),
                 conversion_prob, x, y, z (kpc).
    """
    data = np.loadtxt(filepath)
    assert data.ndim == 2, f"Expected 2D array, got shape {data.shape}"
    assert data.shape[1] == 7, (
        f"Expected 7 columns [NS_idx, flux_weight, photon_energy, "
        f"conv_prob, x, y, z], got {data.shape[1]}"
    )
    n_photons = len(data)
    n_ns = len(np.unique(data[:, 0]))
    print(f"  Loaded {n_photons} ray-tracing samples from {n_ns} neutron stars")
    return data


def energy_to_freq_hz(energy_ev: float) -> float:
    """Convert photon energy in eV to frequency in Hz."""
    return energy_ev / H_EV_S


def freq_hz_to_channel(freq_hz: float) -> int:
    """Map frequency in Hz to nearest MeerKAT channel index (0-based)."""
    chan = round((freq_hz - FIRST_FREQ_HZ) / CHAN_WIDTH_HZ)
    return int(np.clip(chan, 0, TOTAL_CHANS - 1))


def channel_to_freq_hz(chan: int) -> float:
    """Map channel index to center frequency in Hz."""
    return FIRST_FREQ_HZ + chan * CHAN_WIDTH_HZ


# =============================================================================
# Template generation
# =============================================================================

@dataclass
class SignalTemplate:
    """Per-channel signal template for one axion mass + population realization."""
    mass_ev: float
    coupling_gev: float
    pop_model: str
    pop_idx: int
    # Arrays indexed by MeerKAT channel
    channel_indices: np.ndarray   # shape (N_active,) -- channels with nonzero signal
    freq_hz: np.ndarray           # shape (N_active,) -- center frequencies
    flux_density_jy: np.ndarray   # shape (N_active,) -- flux density in Jy
    n_photons_per_chan: np.ndarray # shape (N_active,) -- number of ray samples per channel
    # Summary
    total_flux_jy_hz: float       # total integrated flux weight
    n_ns_contributing: int        # number of NSs with nonzero flux


def build_template(data: np.ndarray, mass_ev: float,
                   coupling_gev: float = G_REF,
                   pop_model: str = 'young',
                   pop_idx: int = 0) -> SignalTemplate:
    """
    Bin photon energies into MeerKAT channels and compute per-channel flux density.

    The flux_weight column in Combined_Flux.dat is in Jy-Hz (integrated flux
    across the line for that ray). To get flux density in a channel of width
    delta_nu, we sum flux_weights of all rays landing in that channel and
    divide by delta_nu.

    Coupling scaling: flux ~ g^4, so
        S(g) = S(g_ref) * (g / g_ref)^4

    Parameters
    ----------
    data : ndarray
        Combined_Flux.dat contents (N, 7).
    mass_ev : float
        Axion mass in eV.
    coupling_gev : float
        Axion-photon coupling in GeV^-1.
    pop_model : str
        Population model name.
    pop_idx : int
        Population realization index.

    Returns
    -------
    SignalTemplate
    """
    # Convert photon energies to frequencies
    photon_energies_ev = data[:, 2]
    freqs_hz = photon_energies_ev / H_EV_S

    # Map to MeerKAT channels
    channels = np.array([freq_hz_to_channel(f) for f in freqs_hz])

    # Filter to rays within the MeerKAT band
    in_band = (freqs_hz >= FIRST_FREQ_HZ - 0.5 * CHAN_WIDTH_HZ) & \
              (freqs_hz <= LAST_FREQ_HZ + 0.5 * CHAN_WIDTH_HZ)
    n_out_of_band = np.sum(~in_band)
    if n_out_of_band > 0:
        print(f"  {n_out_of_band}/{len(freqs_hz)} rays fall outside MeerKAT band "
              f"({FIRST_FREQ_HZ/1e6:.1f}-{LAST_FREQ_HZ/1e6:.1f} MHz)")

    flux_weights = data[:, 1]  # Jy-Hz at g_ref

    # Coupling rescaling: flux ~ g^4
    g_scale = (coupling_gev / G_REF) ** 4

    # Bin flux weights into channels
    # Use np.bincount for speed
    flux_per_chan = np.bincount(channels, weights=flux_weights * g_scale,
                               minlength=TOTAL_CHANS)
    count_per_chan = np.bincount(channels, minlength=TOTAL_CHANS)

    # Convert integrated flux (Jy-Hz) to flux density (Jy) by dividing by channel width
    flux_density = flux_per_chan / CHAN_WIDTH_HZ

    # Extract only channels with nonzero signal
    active = np.where(flux_density > 0)[0]

    # Count unique NSs that contribute
    ns_in_band = data[in_band, 0] if np.any(in_band) else np.array([])
    n_ns = len(np.unique(ns_in_band)) if len(ns_in_band) > 0 else 0

    total_flux = np.sum(flux_weights) * g_scale

    template = SignalTemplate(
        mass_ev=mass_ev,
        coupling_gev=coupling_gev,
        pop_model=pop_model,
        pop_idx=pop_idx,
        channel_indices=active,
        freq_hz=np.array([channel_to_freq_hz(c) for c in active]),
        flux_density_jy=flux_density[active],
        n_photons_per_chan=count_per_chan[active].astype(int),
        total_flux_jy_hz=total_flux,
        n_ns_contributing=n_ns,
    )

    return template


def save_template(template: SignalTemplate, outdir: str = TEMPLATE_DIR) -> str:
    """
    Save template to CSV.

    Returns the output filepath.
    """
    os.makedirs(outdir, exist_ok=True)
    mass_label = f"{template.mass_ev:.3e}eV"
    g_label = f"g{template.coupling_gev:.2e}"
    fname = (f"template_{template.pop_model}_pop{template.pop_idx}_"
             f"ma{mass_label}_{g_label}.csv")
    outpath = os.path.join(outdir, fname)

    with open(outpath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['# Signal template for axion search'])
        writer.writerow([f'# mass_ev = {template.mass_ev}'])
        writer.writerow([f'# coupling_gev = {template.coupling_gev}'])
        writer.writerow([f'# pop_model = {template.pop_model}'])
        writer.writerow([f'# pop_idx = {template.pop_idx}'])
        writer.writerow([f'# total_flux_jy_hz = {template.total_flux_jy_hz}'])
        writer.writerow([f'# n_ns_contributing = {template.n_ns_contributing}'])
        writer.writerow([f'# n_active_channels = {len(template.channel_indices)}'])
        writer.writerow(['channel_index', 'freq_hz', 'freq_mhz',
                         'flux_density_jy', 'flux_density_ujy',
                         'n_ray_samples'])
        for i in range(len(template.channel_indices)):
            writer.writerow([
                template.channel_indices[i],
                f"{template.freq_hz[i]:.1f}",
                f"{template.freq_hz[i] / 1e6:.6f}",
                f"{template.flux_density_jy[i]:.6e}",
                f"{template.flux_density_jy[i] * 1e6:.6f}",
                template.n_photons_per_chan[i],
            ])

    print(f"  Template saved to: {outpath}")
    return outpath


# =============================================================================
# Averaging over population realizations
# =============================================================================

def build_averaged_template(mass_ev: float, pop_model: str,
                            coupling_gev: float = G_REF,
                            pop_indices: Optional[List[int]] = None
                            ) -> Tuple[SignalTemplate, List[SignalTemplate]]:
    """
    Build signal templates for all available realizations and compute the average.

    Parameters
    ----------
    mass_ev : float
        Axion mass in eV.
    pop_model : str
        'young' or 'old'.
    coupling_gev : float
        Coupling constant.
    pop_indices : list of int, optional
        Which realizations to include. If None, uses all available.

    Returns
    -------
    (averaged_template, list_of_individual_templates)
    """
    if pop_indices is None:
        n_max = N_REALIZATIONS[pop_model]
        pop_indices = list(range(n_max))

    templates = []
    all_flux_density = np.zeros(TOTAL_CHANS)
    all_count = np.zeros(TOTAL_CHANS, dtype=int)
    n_loaded = 0

    for idx in pop_indices:
        try:
            filepath = find_combined_flux_file(mass_ev, pop_model, idx)
        except FileNotFoundError as e:
            print(f"  Skipping Pop_{idx}: {e}")
            continue

        print(f"\n  Loading {pop_model} Pop_{idx}...")
        data = load_combined_flux(filepath)
        t = build_template(data, mass_ev, coupling_gev, pop_model, idx)
        templates.append(t)

        # Accumulate into full-band arrays
        for i, ch in enumerate(t.channel_indices):
            all_flux_density[ch] += t.flux_density_jy[i]
            all_count[ch] += t.n_photons_per_chan[i]
        n_loaded += 1

    assert n_loaded > 0, (
        f"No realizations loaded for mass={mass_ev}, pop={pop_model}"
    )

    # Average
    all_flux_density /= n_loaded
    active = np.where(all_flux_density > 0)[0]

    avg_template = SignalTemplate(
        mass_ev=mass_ev,
        coupling_gev=coupling_gev,
        pop_model=pop_model,
        pop_idx=-1,  # -1 indicates average
        channel_indices=active,
        freq_hz=np.array([channel_to_freq_hz(c) for c in active]),
        flux_density_jy=all_flux_density[active],
        n_photons_per_chan=(all_count[active] / n_loaded).astype(int),
        total_flux_jy_hz=np.mean([t.total_flux_jy_hz for t in templates]),
        n_ns_contributing=int(np.mean([t.n_ns_contributing for t in templates])),
    )

    print(f"\n  Average over {n_loaded} realizations:")
    print(f"    Active channels: {len(active)}")
    print(f"    Mean total flux: {avg_template.total_flux_jy_hz:.4e} Jy-Hz")

    return avg_template, templates


# =============================================================================
# Diagnostic plots
# =============================================================================

def plot_template_spectrum(template: SignalTemplate,
                           individual_templates: Optional[List[SignalTemplate]] = None,
                           outdir: str = PLOT_DIR) -> str:
    """
    Generate diagnostic plots for a signal template.

    Plot 1: Full spectrum (flux density vs frequency)
    Plot 2: Zoom on the core region
    Plot 3: Histogram of per-channel flux densities
    Plot 4: Cumulative flux vs channel

    Returns the output plot path.
    """
    os.makedirs(outdir, exist_ok=True)

    freq_mhz = template.freq_hz / 1e6
    flux_ujy = template.flux_density_jy * 1e6  # micro-Jy
    nu_0_mhz = template.mass_ev / H_EV_S / 1e6

    # Figure label
    mass_label_str = f"{template.mass_ev * 1e6:.2f}"
    pop_label = template.pop_model.capitalize()
    if template.pop_idx == -1:
        pop_label += " (avg)"
    else:
        pop_label += f" Pop {template.pop_idx}"

    fig, axes = plt.subplots(2, 2, figsize=(18, 12))

    # --- Plot 1: Full spectrum ---
    ax = axes[0, 0]
    if individual_templates is not None and len(individual_templates) > 1:
        for t in individual_templates:
            f_mhz_i = t.freq_hz / 1e6
            s_ujy_i = t.flux_density_jy * 1e6
            ax.plot(f_mhz_i, s_ujy_i, alpha=0.15, color='cornflowerblue', lw=0.5)
    ax.plot(freq_mhz, flux_ujy, color='navy', lw=1.0,
            label=f'{pop_label}')
    ax.axvline(nu_0_mhz, color='maroon', ls='--', alpha=0.6,
               label=r'$\nu_0 = m_a/h = %.1f$ MHz' % nu_0_mhz)
    ax.set_xlabel(r'Frequency [MHz]')
    ax.set_ylabel(r'Flux Density [$\mu$Jy] at $g = %.0e$ GeV$^{-1}$'
                  % template.coupling_gev)
    ax.set_title(r'Signal Template: $m_a = %s\;\mu$eV' % mass_label_str)
    ax.legend(fontsize=12)

    # --- Plot 2: Zoom on core ---
    ax = axes[0, 1]
    # Find the region containing 90% of the flux
    sort_idx = np.argsort(freq_mhz)
    sorted_freq = freq_mhz[sort_idx]
    sorted_flux = flux_ujy[sort_idx]
    cumflux = np.cumsum(sorted_flux)
    if cumflux[-1] > 0:
        cumflux_norm = cumflux / cumflux[-1]
        lo_idx = np.searchsorted(cumflux_norm, 0.01)
        hi_idx = np.searchsorted(cumflux_norm, 0.99)
        lo_freq = sorted_freq[max(0, lo_idx - 5)]
        hi_freq = sorted_freq[min(len(sorted_freq) - 1, hi_idx + 5)]
        # Add some padding
        padding = (hi_freq - lo_freq) * 0.1
        zoom_mask = (freq_mhz >= lo_freq - padding) & (freq_mhz <= hi_freq + padding)
    else:
        zoom_mask = np.ones(len(freq_mhz), dtype=bool)

    if individual_templates is not None and len(individual_templates) > 1:
        for t in individual_templates:
            f_mhz_i = t.freq_hz / 1e6
            s_ujy_i = t.flux_density_jy * 1e6
            zm = (f_mhz_i >= freq_mhz[zoom_mask].min()) & \
                 (f_mhz_i <= freq_mhz[zoom_mask].max())
            if np.any(zm):
                ax.plot(f_mhz_i[zm], s_ujy_i[zm], alpha=0.15,
                        color='cornflowerblue', lw=0.5)
    ax.plot(freq_mhz[zoom_mask], flux_ujy[zoom_mask], color='navy', lw=1.0)
    ax.axvline(nu_0_mhz, color='maroon', ls='--', alpha=0.6)
    ax.set_xlabel(r'Frequency [MHz]')
    ax.set_ylabel(r'Flux Density [$\mu$Jy]')
    ax.set_title(r'Zoom: Central 98\% of Flux')

    # --- Plot 3: Flux density histogram ---
    ax = axes[1, 0]
    nonzero = flux_ujy[flux_ujy > 0]
    if len(nonzero) > 0:
        ax.hist(np.log10(nonzero), bins=50, color='cornflowerblue', alpha=0.8,
                edgecolor='navy')
    ax.set_xlabel(r'$\log_{10}$(Flux Density / $\mu$Jy)')
    ax.set_ylabel(r'Number of Channels')
    ax.set_title(r'Distribution of Per-Channel Flux Densities')

    # --- Plot 4: Cumulative flux ---
    ax = axes[1, 1]
    sorted_flux_desc = np.sort(template.flux_density_jy)[::-1]
    cumulative = np.cumsum(sorted_flux_desc)
    if cumulative[-1] > 0:
        cumulative_pct = cumulative / cumulative[-1] * 100
    else:
        cumulative_pct = np.zeros_like(cumulative)
    ax.plot(np.arange(1, len(cumulative_pct) + 1), cumulative_pct,
            color='forestgreen', lw=2)
    ax.set_xlabel(r'Number of Channels (ranked by flux)')
    ax.set_ylabel(r'Cumulative Flux [\%]')
    ax.set_title(r'Flux Concentration Across Channels')
    ax.axhline(50, color='goldenrod', ls='--', alpha=0.5)
    ax.axhline(90, color='maroon', ls='--', alpha=0.5)
    if len(cumulative_pct) > 0 and cumulative_pct[-1] > 50:
        n50 = np.searchsorted(cumulative_pct, 50) + 1
        n90 = np.searchsorted(cumulative_pct, 90) + 1
        ax.text(n50 * 1.5, 53, f'{n50} channels = 50%', fontsize=12)
        ax.text(n90 * 1.1, 93, f'{n90} channels = 90%', fontsize=12)

    plt.suptitle(
        r'Phase 7 Signal Template: $m_a = %s\;\mu$eV, %s, $g = %.0e$ GeV$^{-1}$'
        % (mass_label_str, pop_label, template.coupling_gev),
        fontsize=16
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    pop_tag = template.pop_model
    pidx = template.pop_idx if template.pop_idx >= 0 else 'avg'
    outpath = os.path.join(
        outdir,
        f'phase7_template_ma{template.mass_ev:.3e}eV_{pop_tag}_pop{pidx}.png'
    )
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Plot saved: {outpath}")
    return outpath


def plot_meerkat_band_overlay(templates: Dict[str, SignalTemplate],
                              outdir: str = PLOT_DIR) -> str:
    """
    Overlay templates for multiple masses on the MeerKAT band.
    Shows where each axion mass's signal falls within 856-1712 MHz.
    """
    os.makedirs(outdir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(16, 6))

    colors = ['cornflowerblue', 'forestgreen', 'maroon', 'goldenrod']

    for i, (label, template) in enumerate(templates.items()):
        freq_mhz = template.freq_hz / 1e6
        flux_ujy = template.flux_density_jy * 1e6
        color = colors[i % len(colors)]
        ax.plot(freq_mhz, flux_ujy, color=color, lw=0.8, label=label)
        nu_0 = template.mass_ev / H_EV_S / 1e6
        ax.axvline(nu_0, color=color, ls='--', alpha=0.4)

    # Show MeerKAT band boundaries
    ax.axvline(856.0, color='gray', ls=':', alpha=0.5, label='MeerKAT L-band edges')
    ax.axvline(LAST_FREQ_HZ / 1e6, color='gray', ls=':', alpha=0.5)

    ax.set_xlabel(r'Frequency [MHz]')
    ax.set_ylabel(r'Flux Density [$\mu$Jy] at $g = 10^{-12}$ GeV$^{-1}$')
    ax.set_title(r'Signal Templates in MeerKAT L-band')
    ax.legend(fontsize=12)
    ax.set_xlim(850, LAST_FREQ_HZ / 1e6 + 10)

    plt.tight_layout()
    outpath = os.path.join(outdir, 'phase7_template_meerkat_band_overlay.png')
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Plot saved: {outpath}")
    return outpath


# =============================================================================
# Diagnostics
# =============================================================================

def print_template_diagnostics(template: SignalTemplate) -> None:
    """Print diagnostic summary of a signal template."""
    nu_0_mhz = template.mass_ev / H_EV_S / 1e6
    chan_width_mhz = CHAN_WIDTH_HZ / 1e6

    print("\n" + "=" * 70)
    print("SIGNAL TEMPLATE DIAGNOSTICS")
    print("=" * 70)
    print(f"  Axion mass:       {template.mass_ev:.4e} eV = {template.mass_ev*1e6:.3f} ueV")
    print(f"  Rest frequency:   {nu_0_mhz:.4f} MHz")
    print(f"  Coupling:         {template.coupling_gev:.2e} GeV^-1")
    print(f"  Population:       {template.pop_model} (Pop {template.pop_idx})")
    print(f"  Channel width:    {chan_width_mhz*1e3:.3f} kHz = {CHAN_WIDTH_HZ:.1f} Hz")
    print()
    print(f"  Active channels:  {len(template.channel_indices)}")
    print(f"  NSs contributing: {template.n_ns_contributing}")
    print(f"  Total flux:       {template.total_flux_jy_hz:.4e} Jy-Hz")

    if len(template.flux_density_jy) > 0:
        peak_flux = np.max(template.flux_density_jy)
        peak_chan = template.channel_indices[np.argmax(template.flux_density_jy)]
        peak_freq = channel_to_freq_hz(peak_chan) / 1e6
        median_flux = np.median(template.flux_density_jy)
        total_flux_density = np.sum(template.flux_density_jy) * CHAN_WIDTH_HZ

        print(f"\n  Peak flux density:   {peak_flux:.4e} Jy = {peak_flux*1e6:.4f} uJy")
        print(f"    at channel {peak_chan} ({peak_freq:.4f} MHz)")
        print(f"  Median flux density: {median_flux:.4e} Jy = {median_flux*1e6:.4f} uJy")
        print(f"  Integrated flux:     {total_flux_density:.4e} Jy*Hz")

        # Frequency spread of signal
        freq_min = template.freq_hz.min() / 1e6
        freq_max = template.freq_hz.max() / 1e6
        bw = freq_max - freq_min
        print(f"\n  Signal frequency range: {freq_min:.4f} -- {freq_max:.4f} MHz")
        print(f"  Signal bandwidth:       {bw:.4f} MHz = {bw*1e3:.1f} kHz")
        print(f"  Signal spans:           {len(template.channel_indices)} channels")

        # Fractional bandwidth
        frac_bw = bw / nu_0_mhz
        print(f"  Fractional bandwidth:   {frac_bw:.6f} = {frac_bw*1e3:.3f} x 10^-3")

        # Check if signal is within MeerKAT band
        in_band = (nu_0_mhz >= 856.0) and (nu_0_mhz <= LAST_FREQ_HZ / 1e6)
        print(f"\n  Rest freq in MeerKAT band: {in_band}")
        if not in_band:
            print("  WARNING: Axion rest frequency is OUTSIDE the MeerKAT L-band!")

    print("=" * 70)


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Phase 7: Generate signal templates from NS population models"
    )
    parser.add_argument('--mass', type=float, required=True,
                        help='Axion mass in eV (e.g., 4.13e-6)')
    parser.add_argument('--coupling', type=float, default=1e-12,
                        help='Axion-photon coupling in GeV^-1 (default: 1e-12)')
    parser.add_argument('--pop_model', choices=['young', 'old'], default='young',
                        help='Population model (default: young)')
    parser.add_argument('--pop_idx', type=int, default=0,
                        help='Population realization index (default: 0)')
    parser.add_argument('--all_realizations', action='store_true',
                        help='Average over all available realizations')
    parser.add_argument('--all_masses', action='store_true',
                        help='Generate templates for all available masses')
    parser.add_argument('--no_plot', action='store_true',
                        help='Skip diagnostic plots')
    args = parser.parse_args()

    print("=" * 70)
    print("Phase 7: Signal Template Generation")
    print("=" * 70)

    if args.all_masses:
        # Generate templates for all masses
        overlay_templates = {}

        for mass_ev in sorted(AVAILABLE_MASSES.keys()):
            # Check if this mass exists for the selected pop model
            try:
                _ = find_combined_flux_file(mass_ev, args.pop_model, 0)
            except FileNotFoundError:
                print(f"\nSkipping mass {mass_ev:.3e} eV: not available for {args.pop_model}")
                continue

            print(f"\n{'='*70}")
            print(f"Processing mass = {mass_ev:.3e} eV")
            print(f"{'='*70}")

            if args.all_realizations:
                avg_t, indiv = build_averaged_template(
                    mass_ev, args.pop_model, args.coupling
                )
                template = avg_t
                individuals = indiv
            else:
                filepath = find_combined_flux_file(mass_ev, args.pop_model, args.pop_idx)
                data = load_combined_flux(filepath)
                template = build_template(data, mass_ev, args.coupling,
                                          args.pop_model, args.pop_idx)
                individuals = None

            print_template_diagnostics(template)
            save_template(template)

            label = r'$m_a = %.2f\;\mu$eV' % (mass_ev * 1e6)
            overlay_templates[label] = template

            if not args.no_plot:
                plot_template_spectrum(template, individuals)

        if not args.no_plot and len(overlay_templates) > 1:
            plot_meerkat_band_overlay(overlay_templates)

    else:
        # Single mass
        mass_ev = args.mass

        # Validate mass
        closest_mass = min(AVAILABLE_MASSES.keys(),
                           key=lambda m: abs(m - mass_ev))
        if abs(closest_mass - mass_ev) / closest_mass > 0.01:
            print(f"WARNING: Requested mass {mass_ev:.3e} not in pre-computed set.")
            print(f"  Available: {sorted(AVAILABLE_MASSES.keys())}")
            print(f"  Using closest: {closest_mass:.3e}")
        mass_ev = closest_mass

        if args.all_realizations:
            avg_template, individual_templates = build_averaged_template(
                mass_ev, args.pop_model, args.coupling
            )
            template = avg_template
            individuals = individual_templates
        else:
            print(f"\nLoading {args.pop_model} Pop_{args.pop_idx} for mass={mass_ev:.3e} eV")
            filepath = find_combined_flux_file(mass_ev, args.pop_model, args.pop_idx)
            data = load_combined_flux(filepath)
            template = build_template(data, mass_ev, args.coupling,
                                      args.pop_model, args.pop_idx)
            individuals = None

        print_template_diagnostics(template)
        save_template(template)

        if not args.no_plot:
            plot_template_spectrum(template, individuals)

    print("\nPhase 7 signal template generation complete.")


if __name__ == '__main__':
    main()
