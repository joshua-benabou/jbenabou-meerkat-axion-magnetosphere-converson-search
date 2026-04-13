# MeerKAT Galactic Center Data Reduction: Processing Summary

**Date**: 2026-04-12 (updated 2026-04-12 evening)
**PI**: Josh Benabou (UC Berkeley, Safdi group)
**Dataset**: MeerKAT L-band, 800-1400 MHz, ~32,768 channels at 26.123 kHz native resolution
**Target**: Sgr A* / Galactic Center
**Science goal**: Axion-photon conversion search in neutron star magnetospheres

---

## Pipeline Overview

The 17 TB monolithic measurement set (MS) is being processed through a multi-phase pipeline on the Lawrencium HPC cluster (UC Berkeley). The core challenge — CASA's file-locking preventing parallel access to the monolithic MS — was solved by serial subband splitting followed by parallel imaging.

```
Phase 0: Reconnaissance & Setup .................. COMPLETE
Phase 1: Subband Splitting (serial) .............. 86/86 COMPLETE (re-splits running for 5 corrupted)
Phase 2: Continuum Subtraction ................... SKIPPED (for now)
Phase 3a: Dirty Channel Imaging (parallel) ....... 73/86 COMPLETE, 13 RUNNING
Phase 3a: FITS Export ............................ 45/86 COMPLETE, 24 recovery jobs submitted
Phase 3b: Cleaned Channel Imaging ................ READY (script written, submit after 3a)
Phase 4: FITS Assembly & Catalog ................. NOT STARTED
```

---

## Phase 0: Reconnaissance & Setup

- Confirmed single SPW with 32,768 channels (856.0-1711.8 MHz)
- Channel width: 26,123.05 Hz (~26.1 kHz)
- Data is pre-calibrated from the MeerKAT archive (DATA column used)
- Software environment: CASA 6.x (modular) via conda (`casa_cookbook` env)
- Imaging: tclean in cube mode (512x512 pixels, 2 arcsec, Briggs robust=0.5)

## Phase 1: Subband Splitting

**Strategy**: Split the monolithic MS into 86 subbands of 383 channels each (~10 MHz per subband), processed serially to avoid file-locking. Once split, all downstream work is fully parallel.

**Method**: `casatasks.split()` with `spw='0:start~end'` channel selection, `datacolumn='data'`.

| Metric | Value |
|--------|-------|
| Total subbands | 86 (indices 0-85) |
| Channels per subband | 383 |
| Completed | 86 |
| Corrupted (re-splitting) | 44, 52, 56, 62, 70 (FilebufIO read errors) |
| Re-split running | 44, 52, 56, 61, 62, 70, 72 |
| Time per subband | 80-113 min (median ~88 min) |
| Total split time | ~125 hours (serial) |
| Subband disk usage | ~12 TB total |

**Note**: No TOPO-to-LSRK regridding was performed. The TOPO-LSRK shift is <0.1 channel width within a single observation, so native TOPO frequencies are preserved. LSRK conversion can be done as metadata post-processing.

## Phase 3a: Dirty Imaging (niter=0)

**Strategy**: Use `tclean` in `specmode='cube'` mode to produce a 383-channel image cube per subband in a single call, then export each channel to individual FITS files.

**Parameters**:
- Image size: 512 x 512 pixels
- Cell size: 2 arcsec
- Weighting: Briggs, robust=0.5
- niter: 0 (dirty images)
- savemodel: 'none' (critical — avoids write locks)
- SLURM: 1 node per subband, 56 cores, 240 GB RAM, lr7 partition

### Imaging Progress

| Status | Count | Subbands |
|--------|-------|----------|
| FITS exported (383/383) | 45 | 002-006, 008, 014-017, 021-025, 032, 035, 041-043, 045-051, 053-055, 057-060, 063-069, 071 |
| Cube exists, no FITS (recovery submitted) | 24 | 000, 001, 007, 009-013, 018-020, 026-031, 033, 034, 036-040 |
| Partial FITS (still imaging) | 3 | 073 (302), 074 (274), 085 (17) |
| Corrupted MS (re-splitting) | 5 | 044, 052, 056, 062, 070 |
| Currently imaging (SLURM) | 13 | 073-085 on lr7 |
| Re-splitting on lr6 | 7 | 44, 52, 56, 61, 62, 70, 72 |

**Total FITS exported**: ~16,567 / 32,768 channels (51%)

### Timing

| Step | Time |
|------|------|
| tclean cube (383 channels) | 37-60 min (median ~50 min) |
| FITS export (383 channels) | 2.1-2.4 min |
| Total per subband | ~52 min |

### Known Issue: exportfits 'channel' Keyword Bug

The initial implementation attempted to use `exportfits(channel=N)` to extract individual channels from the cube. This parameter does not exist in the installed CASA version, causing silent export failure — tclean succeeded but no FITS were produced.

**Fix**: Two-step export using `imsubimage(chans=str(N))` to extract a single-channel image, then `exportfits()` on that. A recovery script (`phase3_export_recovery.py`) was created to re-export FITS from cubes that had already been produced.

### FITS Validation

Spot-checked subbands 023 and 050:

| Check | Result |
|-------|--------|
| Pixel data distinct across channels | PASS |
| No NaN values | PASS |
| No all-zero images | PASS |
| Shape | (1, 1, 512, 512) — correct |
| Units | Jy/beam — correct |
| Noise level (mid-channel) | ~0.09 Jy/beam |
| Channel 0 anomaly | Peak 3x higher than other channels (edge effect) |

**Frequency header note**: FITS `CRVAL3` is set to the tclean cube reference frequency (first channel), not the per-channel frequency. The actual channel frequency is encoded in the filename (`chan_NNNN_FFFF.FFFMHz.fits`) which is taken directly from the MS metadata. The WCS can recover frequency via `CRVAL3 + (1 - CRPIX3) * CDELT3`, but this gives the tclean grid frequency which is offset by ~88 kHz from the true MS channel frequency. For science use, **read frequency from the filename**.

---

## Disk Usage

| Component | Size |
|-----------|------|
| Original MS | ~17 TB |
| Compressed archive | ~12 TB |
| Split subbands | ~12 TB |
| Images (FITS + CASA products) | ~73 GB |
| **Total** | **~41 TB** |

## Outstanding Work

1. **Finish dirty imaging**: 13 jobs running (subbands 73-85), ~30 min remaining
2. **Export recovery**: Submitted (job 21918030) for 24 subbands with cubes but no FITS
3. **Re-split corrupted subbands**: Jobs running/submitted for 44, 52, 56, 61, 62, 70, 72 — then need re-imaging
4. **Cleaned imaging (Phase 3b)**: Script ready (`phase3_clean_cube.py`, niter=10000, auto-multithresh). Submit after all dirty imaging + export is done: `sbatch scripts/phase3_clean_submit.sh`
5. **Continuum subtraction (Phase 2)**: Skipped for dirty images; may be needed before cleaned imaging
6. **FITS assembly & frequency catalog (Phase 4)**: Create master mapping of filename to frequency
7. **Investigate channel 0 anomaly**: First channel of each subband has elevated flux — may need flagging

## Cluster Resources Used

- **Partition**: lr7 (56 cores, ~240 GB RAM per node)
- **Concurrent jobs**: Up to 86 array tasks (one per subband)
- **Account**: pc_heptheory
- **Software**: CASA 6.x via `/global/home/users/osning/.conda/envs/casa_cookbook`

## Key Scripts

All in `scripts/`:

| Script | Purpose |
|--------|---------|
| `phase1_split_subbands.py` | Serial subband splitting from monolithic MS |
| `phase1_submit.sh` | SLURM submission for Phase 1 |
| `phase3_image_cube.py` | Cube-mode tclean + FITS export per subband |
| `phase3_cube_submit.sh` | SLURM array submission for Phase 3 |
| `phase3_export_recovery.py` | Re-export FITS from existing cubes (bug fix) |
| `phase3_export_recovery_submit.sh` | SLURM submission for export recovery |
| `phase3_image_channels.py` | Per-channel tclean with multiprocessing (alternative approach) |
| `phase3_clean_cube.py` | Cleaned imaging (niter=10000, auto-multithresh) |
| `phase3_clean_submit.sh` | SLURM submission for cleaned imaging |
| `validate_fits.py` | FITS validation script |
