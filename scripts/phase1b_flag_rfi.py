#!/usr/bin/env python3
"""
Phase 1b: RFI flagging on subband MS files.

Runs after Phase 1 (splitting) and before Phase 2 (continuum subtraction).
Applies known MeerKAT L-band RFI masks and optionally runs automated flagging.

Known MeerKAT L-band RFI:
  - 900-960 MHz: GSM cellular
  - 1160-1310 MHz: GPS L2, GLONASS, Galileo, aircraft DME
  - 1525-1559 MHz: Inmarsat, Thuraya satellite downlinks
  - 1559-1610 MHz: GPS L1, GLONASS L1, Galileo E1
  - 1610-1618 MHz: Iridium satellite uplinks

Usage:
    Called by phase1b_submit.sh via SLURM array job.
    SLURM_ARRAY_TASK_ID selects which subband to process.
"""

import os
import sys
import time
import numpy as np
import casatools
from casatasks import flagdata

# === Configuration ===
SUBBANDS_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/subbands'
)

# Known MeerKAT L-band RFI frequency ranges (MHz)
# Conservative ranges — adjust based on inspection
KNOWN_RFI_RANGES_MHZ = [
    (900.0, 960.0),    # GSM cellular
    (1160.0, 1200.0),  # GPS L5, aircraft DME
    (1217.0, 1237.0),  # GPS L2
    (1525.0, 1559.0),  # Satellite downlinks
    (1559.0, 1610.0),  # GPS L1 / GLONASS L1
    (1610.0, 1618.0),  # Iridium
]


def main():
    subband_idx = int(os.environ.get('SLURM_ARRAY_TASK_ID', sys.argv[1] if len(sys.argv) > 1 else '0'))
    input_ms = os.path.join(SUBBANDS_DIR, f'subband_{subband_idx:03d}.ms')

    if not os.path.exists(input_ms):
        print(f"Subband {subband_idx}: {input_ms} does not exist, skipping.")
        return

    # Get frequency range of this subband
    msmd = casatools.msmetadata()
    msmd.open(input_ms)
    freqs = msmd.chanfreqs(0)
    nchan = msmd.nchan(0)
    msmd.close()

    freq_lo_mhz = freqs[0] / 1e6
    freq_hi_mhz = freqs[-1] / 1e6

    print(f"=== Phase 1b: RFI flagging for subband {subband_idx} ===")
    print(f"  Input: {input_ms}")
    print(f"  Channels: {nchan}, freq: {freq_lo_mhz:.3f}-{freq_hi_mhz:.3f} MHz")
    print(f"  Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    t0 = time.time()

    # 1) Save current flag state as backup
    flagdata(vis=input_ms, mode='save', versionname='before_rfi_flagging')

    # 2) Apply known RFI masks
    flagged_ranges = 0
    for rfi_lo, rfi_hi in KNOWN_RFI_RANGES_MHZ:
        # Check if this RFI range overlaps with our subband
        if rfi_hi < freq_lo_mhz or rfi_lo > freq_hi_mhz:
            continue

        # Find channels within the RFI range
        rfi_mask = (freqs / 1e6 >= rfi_lo) & (freqs / 1e6 <= rfi_hi)
        rfi_chans = np.where(rfi_mask)[0]

        if len(rfi_chans) == 0:
            continue

        chan_sel = f'0:{rfi_chans[0]}~{rfi_chans[-1]}'
        print(f"  Flagging known RFI: {rfi_lo}-{rfi_hi} MHz -> channels {rfi_chans[0]}-{rfi_chans[-1]} ({len(rfi_chans)} chans)")

        flagdata(
            vis=input_ms,
            mode='manual',
            spw=chan_sel,
            action='apply',
        )
        flagged_ranges += 1

    # 3) Run tfcrop (time-frequency automatic flagger) for remaining RFI
    print(f"  Running tfcrop automatic flagging...")
    flagdata(
        vis=input_ms,
        mode='tfcrop',
        datacolumn='data',
        action='apply',
        display='',
        flagbackup=False,
    )

    # 4) Report flag statistics
    flag_summary = flagdata(vis=input_ms, mode='summary')
    flagged_frac = flag_summary['flagged'] / flag_summary['total'] * 100
    print(f"  Flag summary: {flagged_frac:.1f}% flagged ({flag_summary['flagged']:.0f}/{flag_summary['total']:.0f})")

    elapsed = time.time() - t0
    print(f"  Known RFI ranges flagged: {flagged_ranges}")
    print(f"  Elapsed: {elapsed:.1f}s ({elapsed/60:.1f} min)")
    print(f"  End time: {time.strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == '__main__':
    main()
