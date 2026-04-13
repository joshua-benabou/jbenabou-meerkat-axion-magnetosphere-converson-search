#!/usr/bin/env python3
"""
Phase 1: Split monolithic MS into subbands using CASA split.

NO LSRK regridding — data stays in TOPO frame. LSRK frequencies are computed
as metadata afterward. This is valid because within a single observation
(~9 hours), the TOPO->LSRK shift variation is <0.1 channel width.

Each SLURM array task extracts one subband (383 channels, ~10 MHz) from the
monolithic MS, selecting only the SgrA* target field.

Usage:
    Called by phase1_submit.sh via SLURM array job.
    Environment variable SLURM_ARRAY_TASK_ID selects which subband to extract.
"""

import os
import sys
import time
import numpy as np
import casatools
from casatasks import split

# === Configuration ===
MS_PATH = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'SCI-20210212-SS-01/17TBdataset/scratch/kat/1621534878_20231022T19_11_29/'
    '1621534878_sdp_l0.ms'
)
OUTPUT_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/subbands'
)
FIELD = 'SgrA*'
DATACOLUMN = 'data'

# Subband parameters
CHANS_PER_SUBBAND = 383      # ~10 MHz at 26.123 kHz/chan
TOTAL_CHANS = 32768
CHAN_WIDTH_KHZ = 26.123
FIRST_FREQ_MHZ = 856.000     # Channel 0 frequency


def compute_lsrk_freqs(ms_path, spw_selection, field):
    """
    Compute mean LSRK frequency for each selected channel.
    Uses CASA's measures tool to convert TOPO->LSRK for the
    observation midpoint.
    """
    msmd = casatools.msmetadata()
    msmd.open(ms_path)
    topo_freqs = msmd.chanfreqs(0)  # All TOPO frequencies
    trange = msmd.timerangeforobs(0)
    t_mid = (trange['begin']['m0']['value'] + trange['end']['m0']['value']) / 2.0
    # Get observatory position
    obs_pos = msmd.observatoryposition()
    # Get source direction
    phase_dir = msmd.phasecenter(0)  # SgrA* direction
    msmd.close()

    # Use measures to compute the TOPO->LSRK conversion at obs midpoint
    me = casatools.measures()
    me.doframe(casatools.measures().epoch('UTC', {'value': t_mid, 'unit': 'd'}))
    me.doframe(obs_pos)
    me.doframe(phase_dir)

    # The velocity correction from TOPO to LSRK
    # v_correction = me.todoppler('radio', me.measure(me.toradialvelocity(...)))
    # Simpler: use the fractional frequency shift
    # For each TOPO freq, LSRK_freq = TOPO_freq * (1 + v_lsr/c)
    # We can get v_lsr from the measures framework

    return topo_freqs  # For now return TOPO; LSRK metadata computed in Phase 4


def main():
    subband_idx = int(os.environ.get('SLURM_ARRAY_TASK_ID', sys.argv[1] if len(sys.argv) > 1 else '0'))

    # Calculate channel range for this subband
    chan_start = subband_idx * CHANS_PER_SUBBAND
    chan_end = min(chan_start + CHANS_PER_SUBBAND - 1, TOTAL_CHANS - 1)
    nchan = chan_end - chan_start + 1

    if chan_start >= TOTAL_CHANS:
        print(f"Subband {subband_idx}: chan_start={chan_start} >= TOTAL_CHANS={TOTAL_CHANS}, skipping.")
        return

    freq_start_mhz = FIRST_FREQ_MHZ + chan_start * CHAN_WIDTH_KHZ / 1000.0
    freq_end_mhz = FIRST_FREQ_MHZ + chan_end * CHAN_WIDTH_KHZ / 1000.0

    spw_selection = f'0:{chan_start}~{chan_end}'
    output_ms = os.path.join(OUTPUT_DIR, f'subband_{subband_idx:03d}.ms')

    print(f"=== Subband {subband_idx} ===")
    print(f"  Channels: {chan_start}-{chan_end} ({nchan} chans)")
    print(f"  Frequency: {freq_start_mhz:.3f}-{freq_end_mhz:.3f} MHz")
    print(f"  SPW selection: {spw_selection}")
    print(f"  Output: {output_ms}")
    print(f"  Method: split (no LSRK regridding)")
    print(f"  Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    if os.path.exists(output_ms):
        # Verify existing output
        try:
            msmd = casatools.msmetadata()
            msmd.open(output_ms)
            existing_nchan = msmd.nchan(0)
            msmd.close()
            if existing_nchan == nchan:
                print(f"  Output already exists with {existing_nchan} channels. Skipping.")
                return
            else:
                print(f"  Output exists but has {existing_nchan} channels (expected {nchan}). Re-extracting.")
                import shutil
                shutil.rmtree(output_ms)
        except Exception as e:
            print(f"  Output exists but is corrupted ({e}). Re-extracting.")
            import shutil
            shutil.rmtree(output_ms)

    t0 = time.time()

    split(
        vis=MS_PATH,
        outputvis=output_ms,
        field=FIELD,
        spw=spw_selection,
        datacolumn=DATACOLUMN,
    )

    elapsed = time.time() - t0

    # Verify output
    try:
        msmd = casatools.msmetadata()
        msmd.open(output_ms)
        out_nchan = msmd.nchan(0)
        out_freqs = msmd.chanfreqs(0)
        out_nrows = msmd.nrows()
        msmd.close()
        print(f"  Verification: {out_nchan} channels, {out_nrows:.0f} rows, "
              f"freq {out_freqs[0]/1e6:.3f}-{out_freqs[-1]/1e6:.3f} MHz")
    except Exception as e:
        print(f"  WARNING: Verification failed: {e}")

    print(f"  Elapsed: {elapsed:.1f}s ({elapsed/60:.1f} min)")
    print(f"  End time: {time.strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == '__main__':
    main()
