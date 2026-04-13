#!/usr/bin/env python3
"""
Phase 7: Axion Search via Sideband Background Subtraction

Searches for axion-photon conversion signals from neutron star magnetospheres
in cleaned MeerKAT channel images. For each target NS at a predicted axion
frequency, the pipeline:

1. Locates the closest channel FITS image to the predicted frequency
2. Builds a background model from sideband channels (nearby in frequency,
   excluding a guard band around the signal channel)
3. Subtracts the background from the signal channel
4. Measures the residual at the NS position and computes significance
5. Flags candidates above a discovery threshold (default 5 sigma)

The sideband approach naturally removes continuum emission and any
frequency-smooth backgrounds without requiring visibility-domain
continuum subtraction.

Usage:
    python phase7_axion_search.py [--n_sideband 50] [--n_guard 5] [--threshold 5.0]

    Or via SLURM: sbatch phase7_submit.sh
"""

import os
import sys
import glob
import argparse
import warnings
import csv
import numpy as np
from dataclasses import dataclass, field, asdict
from typing import List, Tuple, Optional, Dict

from astropy.io import fits
from astropy.wcs import WCS

# === Paths ===
PROJECT_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project'
)
IMAGES_DIR = os.path.join(PROJECT_DIR, 'images')
RESULTS_DIR = os.path.join(PROJECT_DIR, 'results', 'phase7_axion_search')
RFI_FLAG_FILE = os.path.join(PROJECT_DIR, 'rfi_channel_flags.csv')

# === Dataset parameters ===
N_SUBBANDS = 86
CHANS_PER_SUBBAND = 383
TOTAL_CHANS = 32768
CHAN_WIDTH_KHZ = 26.123
FIRST_FREQ_MHZ = 856.000  # Channel 0 frequency (TOPO)
IMSIZE = 512
PIXEL_SCALE_ARCSEC = 2.0


# ============================================================================
# Data classes
# ============================================================================

@dataclass
class NSTarget:
    """A neutron star target from the template bank."""
    name: str
    ra_deg: float
    dec_deg: float
    predicted_freq_mhz: float
    signal_extent_arcsec: float  # Expected spatial extent of axion signal


@dataclass
class SearchResult:
    """Result of the sideband search for one NS target."""
    ns_name: str
    ra_deg: float
    dec_deg: float
    freq_mhz: float
    subband_idx: int
    channel_idx: int
    peak_residual_flux: float
    local_rms: float
    significance: float
    is_candidate: bool
    n_sideband_used: int
    n_sideband_flagged: int
    residual_fits_path: str


# ============================================================================
# Template bank loader
# ============================================================================

def load_ns_template_bank(filepath: Optional[str] = None) -> List[NSTarget]:
    """
    Load the neutron star template bank.

    Each entry defines a target NS with position, predicted axion frequency,
    and expected signal extent.

    Parameters
    ----------
    filepath : str, optional
        Path to a CSV file with columns: name, ra_deg, dec_deg,
        predicted_freq_mhz, signal_extent_arcsec.
        If None, returns dummy targets for testing.

    Returns
    -------
    list of NSTarget
    """
    if filepath is not None and os.path.exists(filepath):
        targets = []
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                targets.append(NSTarget(
                    name=row['name'],
                    ra_deg=float(row['ra_deg']),
                    dec_deg=float(row['dec_deg']),
                    predicted_freq_mhz=float(row['predicted_freq_mhz']),
                    signal_extent_arcsec=float(row['signal_extent_arcsec']),
                ))
        return targets

    # TODO: Replace with real implementation -- load Sam Witte's template bank
    # The real template bank will come as a file from Sam with NS positions,
    # predicted axion conversion frequencies, and signal morphology parameters.
    print("WARNING: Using dummy NS template bank for testing.", flush=True)
    return [
        NSTarget("PSR_J1745-2900", 266.41684, -29.00781, 920.5, 12.0),
        NSTarget("PSR_J1746-2850", 266.56250, -28.83889, 1050.3, 8.0),
        NSTarget("SGR_J1745-2900", 266.41700, -29.00806, 1180.7, 15.0),
        NSTarget("PSR_J1747-2958", 266.82500, -29.97500, 1320.1, 10.0),
        NSTarget("PSR_J1748-2446", 267.02083, -24.77778, 1450.0, 6.0),
    ]


# ============================================================================
# RFI flag map
# ============================================================================

def load_rfi_flags(filepath: str = RFI_FLAG_FILE) -> Dict[Tuple[int, int], bool]:
    """
    Load RFI channel flag map.

    Parameters
    ----------
    filepath : str
        Path to CSV with columns: subband, channel, frequency_mhz, rfi_flag.

    Returns
    -------
    dict mapping (subband_idx, channel_idx) -> True if flagged as RFI.
    """
    if not os.path.exists(filepath):
        # TODO: Replace with real implementation once Phase 5 produces this file
        print(f"WARNING: RFI flag file not found at {filepath}. "
              f"Proceeding with no RFI flags.", flush=True)
        return {}

    flags = {}
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sb = int(row['subband'])
            ch = int(row['channel'])
            flagged = row['rfi_flag'].strip().lower() in ('true', '1', 'yes')
            if flagged:
                flags[(sb, ch)] = True
    return flags


def is_channel_flagged(subband_idx: int, channel_idx: int,
                       rfi_flags: Dict[Tuple[int, int], bool]) -> bool:
    """Check if a (subband, channel) pair is flagged as RFI."""
    return rfi_flags.get((subband_idx, channel_idx), False)


# ============================================================================
# Frequency <-> channel mapping
# ============================================================================

def freq_to_global_channel(freq_mhz: float) -> int:
    """
    Map a frequency in MHz to the nearest global channel index (0-32767).

    Parameters
    ----------
    freq_mhz : float
        Frequency in MHz.

    Returns
    -------
    int
        Global channel index (clipped to valid range).
    """
    chan = round((freq_mhz - FIRST_FREQ_MHZ) / (CHAN_WIDTH_KHZ / 1e3))
    return int(np.clip(chan, 0, TOTAL_CHANS - 1))


def global_channel_to_subband(global_chan: int) -> Tuple[int, int]:
    """
    Map a global channel index to (subband_index, local_channel_index).

    Parameters
    ----------
    global_chan : int
        Global channel index (0-32767).

    Returns
    -------
    tuple of (int, int)
        (subband_idx, channel_within_subband).
    """
    subband_idx = global_chan // CHANS_PER_SUBBAND
    local_chan = global_chan % CHANS_PER_SUBBAND
    # Clip to valid range
    subband_idx = min(subband_idx, N_SUBBANDS - 1)
    if subband_idx == N_SUBBANDS - 1:
        # Last subband may have fewer channels
        max_local = TOTAL_CHANS - (N_SUBBANDS - 1) * CHANS_PER_SUBBAND - 1
        local_chan = min(local_chan, max_local)
    return subband_idx, local_chan


def channel_to_freq_mhz(subband_idx: int, local_chan: int) -> float:
    """
    Compute frequency in MHz for a given subband and local channel.

    Parameters
    ----------
    subband_idx : int
    local_chan : int

    Returns
    -------
    float
        Frequency in MHz.
    """
    global_chan = subband_idx * CHANS_PER_SUBBAND + local_chan
    return FIRST_FREQ_MHZ + global_chan * (CHAN_WIDTH_KHZ / 1e3)


def find_channel_fits(subband_idx: int, local_chan: int,
                      cleaned: bool = True) -> Optional[str]:
    """
    Find the FITS file for a given subband and local channel index.

    Uses the naming convention: images/subband_XXX/[cleaned/]chan_NNNN_FREQ.fits
    where FREQ is formatted as e.g. 856.000MHz.

    Parameters
    ----------
    subband_idx : int
    local_chan : int
    cleaned : bool
        If True, look in the 'cleaned' subdirectory.

    Returns
    -------
    str or None
        Path to the FITS file, or None if not found.
    """
    freq_mhz = channel_to_freq_mhz(subband_idx, local_chan)
    if cleaned:
        subdir = os.path.join(IMAGES_DIR, f'subband_{subband_idx:03d}', 'cleaned')
    else:
        subdir = os.path.join(IMAGES_DIR, f'subband_{subband_idx:03d}')

    # Try exact filename match first
    fname = f'chan_{local_chan:04d}_{freq_mhz:.3f}MHz.fits'
    exact_path = os.path.join(subdir, fname)
    if os.path.exists(exact_path):
        return exact_path

    # Fall back to glob matching the channel number (freq may differ slightly)
    pattern = os.path.join(subdir, f'chan_{local_chan:04d}_*.fits')
    matches = glob.glob(pattern)
    if matches:
        return matches[0]

    return None


def freq_to_fits_path(freq_mhz: float, cleaned: bool = True) -> Optional[str]:
    """
    Map a frequency in MHz to the corresponding FITS file path.

    Parameters
    ----------
    freq_mhz : float
        Target frequency in MHz.
    cleaned : bool
        Whether to use cleaned images.

    Returns
    -------
    str or None
        Path to the FITS file, or None if not found.
    """
    global_chan = freq_to_global_channel(freq_mhz)
    subband_idx, local_chan = global_channel_to_subband(global_chan)
    return find_channel_fits(subband_idx, local_chan, cleaned=cleaned)


# ============================================================================
# Image I/O helpers
# ============================================================================

def load_fits_image(filepath: str) -> Tuple[np.ndarray, fits.Header]:
    """
    Load a 2D image from a FITS file.

    Handles the case where CASA outputs 4D cubes (Stokes, Freq, RA, Dec)
    by squeezing to 2D.

    Parameters
    ----------
    filepath : str

    Returns
    -------
    tuple of (ndarray, Header)
        2D image array and FITS header.
    """
    with fits.open(filepath) as hdul:
        data = hdul[0].data.copy()
        header = hdul[0].header.copy()
    # Squeeze singleton dimensions (CASA outputs [1, 1, NY, NX])
    data = np.squeeze(data)
    if data.ndim != 2:
        raise ValueError(f"Expected 2D image after squeeze, got shape {data.shape} "
                         f"from {filepath}")
    return data, header


def save_fits_image(data: np.ndarray, header: fits.Header,
                    filepath: str) -> None:
    """
    Save a 2D image to a FITS file.

    Parameters
    ----------
    data : ndarray
        2D image array.
    header : Header
        FITS header (will be modified to reflect 2D output).
    filepath : str
    """
    # Ensure we write a 2D FITS (strip NAXIS3, NAXIS4 etc.)
    out_header = header.copy()
    out_header['NAXIS'] = 2
    for key in ['NAXIS3', 'NAXIS4', 'CRPIX3', 'CRPIX4', 'CDELT3', 'CDELT4',
                'CRVAL3', 'CRVAL4', 'CTYPE3', 'CTYPE4', 'CUNIT3', 'CUNIT4']:
        if key in out_header:
            del out_header[key]
    out_header['COMMENT'] = 'Sideband-subtracted residual (Phase 7 axion search)'

    hdu = fits.PrimaryHDU(data=data.astype(np.float32), header=out_header)
    hdu.writeto(filepath, overwrite=True)


# ============================================================================
# Coordinate helpers
# ============================================================================

def radec_to_pixel(ra_deg: float, dec_deg: float,
                   header: fits.Header) -> Tuple[int, int]:
    """
    Convert RA, Dec (degrees) to pixel coordinates using WCS.

    Parameters
    ----------
    ra_deg, dec_deg : float
        Sky coordinates in degrees.
    header : Header
        FITS header with WCS information.

    Returns
    -------
    tuple of (int, int)
        (x_pixel, y_pixel) -- integer pixel coordinates.
    """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        wcs = WCS(header, naxis=2)
    px, py = wcs.world_to_pixel_values(ra_deg, dec_deg)
    return int(round(px)), int(round(py))


def extract_cutout(image: np.ndarray, x_center: int, y_center: int,
                   half_size: int) -> np.ndarray:
    """
    Extract a square cutout from the image, handling edge clipping.

    Parameters
    ----------
    image : ndarray
        2D image.
    x_center, y_center : int
        Center pixel coordinates.
    half_size : int
        Half-width of the cutout box in pixels.

    Returns
    -------
    ndarray
        2D cutout (may be smaller than requested if near edge).
    """
    ny, nx = image.shape
    y_lo = max(0, y_center - half_size)
    y_hi = min(ny, y_center + half_size + 1)
    x_lo = max(0, x_center - half_size)
    x_hi = min(nx, x_center + half_size + 1)
    return image[y_lo:y_hi, x_lo:x_hi]


# ============================================================================
# Core sideband background subtraction
# ============================================================================

def get_sideband_channels(subband_idx: int, local_chan: int,
                          n_sideband: int, n_guard: int,
                          rfi_flags: Dict[Tuple[int, int], bool]
                          ) -> List[Tuple[int, int]]:
    """
    Determine the sideband channel list for background estimation.

    The sideband consists of channels within +/- n_sideband of the signal
    channel, excluding:
    - The guard band (+/- n_guard channels around signal)
    - Channels flagged as RFI
    - Channels outside the valid range for this subband or neighboring subbands

    Sideband channels can span adjacent subbands if the signal channel is
    near a subband boundary.

    Parameters
    ----------
    subband_idx : int
        Subband index of the signal channel.
    local_chan : int
        Local channel index within the subband.
    n_sideband : int
        Number of channels on each side to include in the sideband.
    n_guard : int
        Number of channels on each side to exclude as guard band.
    rfi_flags : dict
        RFI flag map.

    Returns
    -------
    list of (subband_idx, local_channel_idx) tuples for valid sideband channels.
    """
    global_signal = subband_idx * CHANS_PER_SUBBAND + local_chan

    sideband = []
    for offset in range(-n_sideband, n_sideband + 1):
        # Skip guard band and signal channel itself
        if abs(offset) <= n_guard:
            continue

        global_chan = global_signal + offset
        if global_chan < 0 or global_chan >= TOTAL_CHANS:
            continue

        sb, lc = global_channel_to_subband(global_chan)

        # Skip RFI-flagged channels
        if is_channel_flagged(sb, lc, rfi_flags):
            continue

        sideband.append((sb, lc))

    return sideband


def build_background_model(sideband_channels: List[Tuple[int, int]],
                           cleaned: bool = True,
                           method: str = 'median'
                           ) -> Tuple[Optional[np.ndarray], int]:
    """
    Build a pixel-by-pixel background model from sideband channel images.

    Parameters
    ----------
    sideband_channels : list of (subband_idx, local_chan) tuples
        Channels to include in the background model.
    cleaned : bool
        Whether to use cleaned images.
    method : str
        'median' or 'mean'. Median is more robust to outliers.

    Returns
    -------
    tuple of (ndarray or None, int)
        Background model image (2D) and number of channels actually loaded.
        Returns None if no channels could be loaded.
    """
    images = []
    for sb, lc in sideband_channels:
        fpath = find_channel_fits(sb, lc, cleaned=cleaned)
        if fpath is None:
            continue
        try:
            img, _ = load_fits_image(fpath)
            images.append(img)
        except Exception as e:
            print(f"  WARNING: Could not load {fpath}: {e}", flush=True)
            continue

    if len(images) == 0:
        return None, 0

    stack = np.array(images)  # shape: (N_chan, NY, NX)

    if method == 'median':
        background = np.nanmedian(stack, axis=0)
    elif method == 'mean':
        background = np.nanmean(stack, axis=0)
    else:
        raise ValueError(f"Unknown method: {method}")

    return background, len(images)


def sideband_subtract(signal_image: np.ndarray,
                      background_model: np.ndarray) -> np.ndarray:
    """
    Subtract the sideband background model from the signal channel image.

    The residual should contain only frequency-dependent emission:
    anything that varies between the signal channel and the sideband average,
    including potential axion conversion signals.

    Parameters
    ----------
    signal_image : ndarray
        2D image of the signal channel.
    background_model : ndarray
        2D pixel-by-pixel background model from sideband channels.

    Returns
    -------
    ndarray
        2D residual image.
    """
    return signal_image - background_model


# ============================================================================
# Detection / significance measurement
# ============================================================================

def measure_local_rms(residual: np.ndarray, x_center: int, y_center: int,
                      exclude_radius_pix: int = 20,
                      annulus_width_pix: int = 50) -> float:
    """
    Compute local RMS from a residual image in an annulus around the target,
    excluding the target region.

    Parameters
    ----------
    residual : ndarray
        2D residual image.
    x_center, y_center : int
        Pixel position of the NS target.
    exclude_radius_pix : int
        Radius (pixels) around the target to exclude from RMS calculation.
    annulus_width_pix : int
        Width of the annulus beyond the exclusion zone.

    Returns
    -------
    float
        Local RMS in the annulus region.
    """
    ny, nx = residual.shape
    yy, xx = np.ogrid[:ny, :nx]
    r2 = (xx - x_center)**2 + (yy - y_center)**2

    inner_r2 = exclude_radius_pix**2
    outer_r2 = (exclude_radius_pix + annulus_width_pix)**2

    annulus_mask = (r2 > inner_r2) & (r2 <= outer_r2)
    annulus_pixels = residual[annulus_mask]

    # Remove NaNs
    annulus_pixels = annulus_pixels[np.isfinite(annulus_pixels)]

    if len(annulus_pixels) < 10:
        # Fall back to global RMS if annulus is too small
        valid = residual[np.isfinite(residual)]
        if len(valid) == 0:
            return np.nan
        return np.std(valid)

    return np.std(annulus_pixels)


def measure_peak_flux(residual: np.ndarray, x_center: int, y_center: int,
                      aperture_radius_pix: int = 10) -> float:
    """
    Measure peak flux in the residual within an aperture around the NS position.

    Parameters
    ----------
    residual : ndarray
        2D residual image.
    x_center, y_center : int
        Pixel position of the NS target.
    aperture_radius_pix : int
        Radius (pixels) to search for peak.

    Returns
    -------
    float
        Peak pixel value within the aperture.
    """
    cutout = extract_cutout(residual, x_center, y_center, aperture_radius_pix)
    if cutout.size == 0:
        return np.nan
    valid = cutout[np.isfinite(cutout)]
    if len(valid) == 0:
        return np.nan
    return float(np.max(valid))


def compute_significance(peak_flux: float, local_rms: float) -> float:
    """
    Compute detection significance as peak flux / local RMS.

    Parameters
    ----------
    peak_flux : float
    local_rms : float

    Returns
    -------
    float
        Significance in units of sigma. Returns 0 if rms is invalid.
    """
    if not np.isfinite(local_rms) or local_rms <= 0:
        return 0.0
    return peak_flux / local_rms


# ============================================================================
# Main search pipeline
# ============================================================================

def search_single_target(target: NSTarget,
                         rfi_flags: Dict[Tuple[int, int], bool],
                         n_sideband: int = 50,
                         n_guard: int = 5,
                         threshold_sigma: float = 5.0,
                         bg_method: str = 'median',
                         use_cleaned: bool = True,
                         save_residuals: bool = True,
                         ) -> Optional[SearchResult]:
    """
    Run sideband background subtraction and detection for a single NS target.

    Parameters
    ----------
    target : NSTarget
        The neutron star target.
    rfi_flags : dict
        RFI flag map.
    n_sideband : int
        Number of sideband channels on each side.
    n_guard : int
        Number of guard channels on each side.
    threshold_sigma : float
        Discovery threshold in sigma.
    bg_method : str
        Background estimation method ('median' or 'mean').
    use_cleaned : bool
        Whether to use cleaned images.
    save_residuals : bool
        Whether to save residual FITS files.

    Returns
    -------
    SearchResult or None
        Result of the search, or None if the signal channel is unavailable.
    """
    print(f"\n--- Searching target: {target.name} ---", flush=True)
    print(f"  Position: RA={target.ra_deg:.5f}, Dec={target.dec_deg:.5f}", flush=True)
    print(f"  Predicted freq: {target.predicted_freq_mhz:.3f} MHz", flush=True)

    # Map frequency to channel
    global_chan = freq_to_global_channel(target.predicted_freq_mhz)
    subband_idx, local_chan = global_channel_to_subband(global_chan)
    actual_freq = channel_to_freq_mhz(subband_idx, local_chan)
    print(f"  Mapped to subband {subband_idx}, channel {local_chan} "
          f"(actual freq: {actual_freq:.3f} MHz)", flush=True)

    # Check if signal channel is RFI-flagged
    if is_channel_flagged(subband_idx, local_chan, rfi_flags):
        print(f"  WARNING: Signal channel is RFI-flagged! Skipping.", flush=True)
        return None

    # Load signal channel image
    signal_path = find_channel_fits(subband_idx, local_chan, cleaned=use_cleaned)
    if signal_path is None:
        print(f"  ERROR: Signal channel FITS not found.", flush=True)
        return None
    print(f"  Signal image: {signal_path}", flush=True)

    signal_image, header = load_fits_image(signal_path)

    # Get NS pixel position
    x_ns, y_ns = radec_to_pixel(target.ra_deg, target.dec_deg, header)
    print(f"  NS pixel position: ({x_ns}, {y_ns})", flush=True)

    # Check if target is within image bounds
    if x_ns < 0 or x_ns >= IMSIZE or y_ns < 0 or y_ns >= IMSIZE:
        print(f"  WARNING: Target outside image bounds ({IMSIZE}x{IMSIZE}). "
              f"Skipping.", flush=True)
        return None

    # Build sideband channel list
    sideband_channels = get_sideband_channels(
        subband_idx, local_chan, n_sideband, n_guard, rfi_flags
    )
    n_total_sideband = 2 * (n_sideband - n_guard)  # max possible
    n_flagged = n_total_sideband - len(sideband_channels)
    print(f"  Sideband: {len(sideband_channels)} channels "
          f"({n_flagged} flagged/excluded)", flush=True)

    if len(sideband_channels) < 10:
        print(f"  WARNING: Too few sideband channels (<10). "
              f"Result may be unreliable.", flush=True)

    # Build background model
    background, n_loaded = build_background_model(
        sideband_channels, cleaned=use_cleaned, method=bg_method
    )
    if background is None:
        print(f"  ERROR: Could not build background model (no images loaded).",
              flush=True)
        return None
    print(f"  Background model built from {n_loaded} channels "
          f"(method: {bg_method})", flush=True)

    # Sideband subtraction
    residual = sideband_subtract(signal_image, background)

    # Measure detection significance
    signal_extent_pix = max(1, int(round(
        target.signal_extent_arcsec / PIXEL_SCALE_ARCSEC
    )))
    peak_flux = measure_peak_flux(residual, x_ns, y_ns,
                                  aperture_radius_pix=signal_extent_pix)
    local_rms = measure_local_rms(residual, x_ns, y_ns,
                                  exclude_radius_pix=signal_extent_pix + 5)
    significance = compute_significance(peak_flux, local_rms)

    print(f"  Peak residual flux: {peak_flux:.6f} Jy/beam", flush=True)
    print(f"  Local RMS: {local_rms:.6f} Jy/beam", flush=True)
    print(f"  Significance: {significance:.2f} sigma", flush=True)

    is_candidate = significance >= threshold_sigma
    if is_candidate:
        print(f"  *** CANDIDATE DETECTED (>{threshold_sigma} sigma)! ***", flush=True)

    # Save residual image
    residual_path = ""
    if save_residuals:
        os.makedirs(RESULTS_DIR, exist_ok=True)
        residual_fname = (f"residual_{target.name}_sb{subband_idx:03d}_"
                          f"ch{local_chan:04d}_{actual_freq:.3f}MHz.fits")
        residual_path = os.path.join(RESULTS_DIR, residual_fname)
        save_fits_image(residual, header, residual_path)
        print(f"  Residual saved: {residual_path}", flush=True)

    return SearchResult(
        ns_name=target.name,
        ra_deg=target.ra_deg,
        dec_deg=target.dec_deg,
        freq_mhz=actual_freq,
        subband_idx=subband_idx,
        channel_idx=local_chan,
        peak_residual_flux=peak_flux,
        local_rms=local_rms,
        significance=significance,
        is_candidate=is_candidate,
        n_sideband_used=n_loaded,
        n_sideband_flagged=n_flagged,
        residual_fits_path=residual_path,
    )


def run_search(targets: List[NSTarget],
               rfi_flags: Dict[Tuple[int, int], bool],
               n_sideband: int = 50,
               n_guard: int = 5,
               threshold_sigma: float = 5.0,
               bg_method: str = 'median',
               use_cleaned: bool = True,
               ) -> List:
    """
    Run the full sideband search pipeline across all NS targets.

    Returns list of SearchResult objects.
    """
    results = []
    for target in targets:
        result = search_single_target(
            target, rfi_flags,
            n_sideband=n_sideband,
            n_guard=n_guard,
            threshold_sigma=threshold_sigma,
            bg_method=bg_method,
            use_cleaned=use_cleaned,
        )
        if result is not None:
            results.append(result)

    if not results:
        print("\nNo results to report (all targets skipped or failed).", flush=True)

    return results


# ============================================================================
# Entry point
# ============================================================================

def main():
    """
    Run the Phase 7 axion search pipeline.

    Loads the NS template bank and RFI flags, runs sideband background
    subtraction for each target, and saves results.
    """
    parser = argparse.ArgumentParser(
        description="Phase 7: Axion search via sideband background subtraction"
    )
    parser.add_argument('--n_sideband', type=int, default=50,
                        help='Number of sideband channels on each side (default: 50)')
    parser.add_argument('--n_guard', type=int, default=5,
                        help='Number of guard channels on each side (default: 5)')
    parser.add_argument('--threshold', type=float, default=5.0,
                        help='Discovery threshold in sigma (default: 5.0)')
    parser.add_argument('--bg_method', choices=['median', 'mean'], default='median',
                        help='Background estimation method (default: median)')
    parser.add_argument('--template_bank', type=str, default=None,
                        help='Path to NS template bank CSV (default: dummy targets)')
    parser.add_argument('--use_dirty', action='store_true',
                        help='Use dirty images instead of cleaned')
    args = parser.parse_args()

    print("=" * 70, flush=True)
    print("Phase 7: Axion Search via Sideband Background Subtraction", flush=True)
    print("=" * 70, flush=True)
    print(f"  n_sideband: {args.n_sideband}", flush=True)
    print(f"  n_guard:    {args.n_guard}", flush=True)
    print(f"  threshold:  {args.threshold} sigma", flush=True)
    print(f"  bg_method:  {args.bg_method}", flush=True)
    print(f"  images:     {'dirty' if args.use_dirty else 'cleaned'}", flush=True)
    print(f"  results:    {RESULTS_DIR}", flush=True)
    print("", flush=True)

    # Load inputs
    targets = load_ns_template_bank(args.template_bank)
    print(f"Loaded {len(targets)} NS targets.", flush=True)

    rfi_flags = load_rfi_flags()
    n_rfi = len(rfi_flags)
    print(f"Loaded {n_rfi} RFI-flagged channels.", flush=True)

    # Run search
    use_cleaned = not args.use_dirty
    results = run_search(
        targets, rfi_flags,
        n_sideband=args.n_sideband,
        n_guard=args.n_guard,
        threshold_sigma=args.threshold,
        bg_method=args.bg_method,
        use_cleaned=use_cleaned,
    )

    # Save results
    if results:
        os.makedirs(RESULTS_DIR, exist_ok=True)
        results_csv = os.path.join(RESULTS_DIR, 'search_results.csv')
        fieldnames = list(asdict(results[0]).keys())
        with open(results_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in results:
                writer.writerow(asdict(r))

        print(f"\n{'=' * 70}", flush=True)
        print("RESULTS SUMMARY", flush=True)
        print(f"{'=' * 70}", flush=True)
        print(f"{'NS Name':<25} {'Freq MHz':>10} {'Peak Flux':>12} "
              f"{'Local RMS':>12} {'Signif':>8} {'Cand?':>6}", flush=True)
        print("-" * 75, flush=True)
        for r in results:
            print(f"{r.ns_name:<25} {r.freq_mhz:>10.3f} "
                  f"{r.peak_residual_flux:>12.6f} {r.local_rms:>12.6f} "
                  f"{r.significance:>8.2f} {str(r.is_candidate):>6}", flush=True)
        print(f"\nResults saved to: {results_csv}", flush=True)

        n_candidates = sum(1 for r in results if r.is_candidate)
        print(f"\nCandidates (>{args.threshold} sigma): {n_candidates} / "
              f"{len(results)}", flush=True)
    else:
        print("\nNo results produced.", flush=True)

    print(f"\nPhase 7 complete.", flush=True)


if __name__ == '__main__':
    main()
