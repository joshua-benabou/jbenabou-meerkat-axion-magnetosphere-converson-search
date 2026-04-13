# Session Handoff

**Last updated**: 2026-04-13 (end-of-session update)

## Current Pipeline State

### Phase 1: Subband Splitting -- COMPLETE
- All 86 subbands split from the monolithic 17 TB MS (subband_000 through subband_085)
- Located in `subbands/subband_XXX.ms`
- 383 channels per subband, ~26 kHz resolution
- Total: ~12 TB of split subbands

### Phase 3a: Dirty Imaging -- COMPLETE
- All 86 subbands imaged with niter=0 (dirty images)
- 32,768 total channel FITS files across all subbands
- Located in `images/subband_XXX/chan_NNNN_FREQ.fits`
- 512x512 pixels, 2 arcsec/pixel, Briggs robust=0.5
- Channel 0 of each subband has ~3x higher peak flux (edge channel artifact)

### Phase 3b: Cleaned Imaging -- FUNDAMENTAL PROBLEM IDENTIFIED
- Previous attempt: auto-multithresh with noisethreshold=5.0 → zero clean components
- **Test cleaning job 21925201** (still running, ~30 min in) revealed:
  - MFS continuum (383 ch averaged): peak=138 mJy, RMS=52.9 mJy, **SNR=2.6**
  - The RMS is dominated by **sidelobes from the bright GC field**, not thermal noise
  - Theoretical sensitivity is 75.6 µJy but actual RMS is 700x higher
  - Even CASA's auto-multithresh mask is **empty** on the continuum
  - Manual sigma-based masks also empty (3σ=159 mJy > peak=138 mJy)
  - Absolute-threshold masks work (10% peak → 31% of pixels) but are not useful for cleaning
  - **Root cause**: The GC is so bright and extended that sidelobes dominate the entire field at 512x512/2" with standard gridding. Self-calibration is needed first.
- **Options**:
  1. Self-cal the continuum first (proper approach, adds complexity)
  2. Accept dirty images for axion search (sideband subtraction removes frequency-smooth sidelobes)
  3. Try WSClean (different auto-masking may work better)
  4. Use natural weighting or larger cell to improve surface brightness sensitivity
- Output in `test_cleaning/subband_030/`

### Phase 5: RFI Flagging -- JOB SUBMITTED
- SLURM job 21925202 on lr6
- Script: `scripts/phase5_rfi_flagging.py`
- Reads all 32,768 FITS (read-only), computes per-channel stats
- Flags by: RMS outliers (3x local median), known RFI bands (GSM-900, GPS L1), band edges
- Output: `rfi_channel_flags.csv` + `plots/rfi_overview.png`

### Phase 6: Sanity Checks -- NOTEBOOK IN PROGRESS
- `phase6_sanity_checks.ipynb` being created (background agent)
- Point source injection, noise assessment, 5-sigma discovery threshold

### Phase 7: Axion Search Pipeline -- SKELETON IN PROGRESS  
- `scripts/phase7_axion_search.py` being created (background agent)
- Real sideband subtraction code, dummy NS template bank
- Waiting on Sam Witte's actual template bank

## Meeting Decisions (2026-04-13, Josh + collaborators)

1. **No continuum subtraction** -- sideband analysis replaces it
2. **Axion analysis uses CLEANED images**
3. **Sideband analysis**: nearby frequency channels' spatial morphology = background model
4. **Sam Witte has template bank of converting neutron stars**
5. **Axion signals resolved by multiple pixels**
6. **RFI flagging**: boolean flags only, don't modify images
7. **Sanity checks**: point source injection, 5-sigma discovery threshold, cleaning fidelity

## GitHub
- Repo: https://github.com/joshua-benabou/jbenabou-meerkat-axion-magnetosphere-converson-search.git
- Auth: HTTPS with personal access token
- `reference/` directory removed from history (private emails/Slack)

## Collaborator Code References
- **Josh Foster SBI**: https://github.com/joshwfoster/AxionLinePop (appears private). Public repo: https://github.com/joshwfoster/RadioAxionSearch
- **Sam Witte population modeling**: delivered via Slack uploads (Analysis.zip, BF_Pop1.zip). No persistent repo.
- **Overleaf**: https://www.overleaf.com/4891248482zfcvpgppyjpp#b1ece9

## Key Technical Details

- **Cluster**: Lawrencium, lr7 (56 cores, 240 GB) for imaging, lr6 for lighter tasks
- **Conda env**: `/global/home/users/osning/.conda/envs/casa_cookbook` (CASA 6.6.0.20)
- **CASA activation**: `eval "$(conda shell.bash hook)" && conda activate /global/home/users/osning/.conda/envs/casa_cookbook`
- **Timing**: ~50 min per dirty cube, ~56 min per clean cube on lr7
- **Export bug**: use `imsubimage(chans=str(N))` + `exportfits()`, not `exportfits(channel=N)`

## Scripts

All in `scripts/`:
- `phase1_split_subbands.py` / `phase1_submit.sh` -- subband splitting
- `phase3_image_cube.py` / `phase3_cube_submit.sh` -- dirty imaging
- `phase3_clean_cube.py` / `phase3_clean_submit.sh` -- cleaned imaging (broken auto-multithresh)
- `phase3_export_recovery.py` / `phase3_export_recovery_submit.sh` -- FITS export fix
- `phase5_rfi_flagging.py` / `phase5_submit.sh` -- RFI channel quality map
- `phase7_axion_search.py` (in progress) -- sideband analysis pipeline skeleton
- `test_cleaning_strategy.py` / `test_cleaning_submit.sh` -- cleaning test (3 strategies)
- `validate_fits.py` -- FITS validation

## Session 2026-04-13: What Was Done

1. Created `handoff.md` and added Session Handoff section to CLAUDE.md
2. Updated `master_plan.md` with meeting decisions (dropped Phase 2, added Phases 5-7)
3. Pushed to GitHub: https://github.com/joshua-benabou/jbenabou-meerkat-axion-magnetosphere-converson-search.git
4. Scrubbed `reference/` from git history (private emails/Slack)
5. **Phase 5 RFI flagging**: Complete. 1,656/32,768 channels flagged (5.1%). Output: `rfi_channel_flags.csv`, `plots/rfi_overview.png`
6. **Phase 6 sanity checks notebook**: Created `phase6_sanity_checks.ipynb` (not yet run)
7. **Phase 7 axion search skeleton**: Created `scripts/phase7_axion_search.py` with real sideband subtraction code
8. **Cleaning test** (job 21925201): Tested 3 strategies on subband_030 -- ALL produced zero clean components
9. Cloned AxionLinePop to `reference/AxionLinePop/` and analyzed the SBI/CASCAL pipeline
10. Analyzed Josh Foster's public RadioAxionSearch repo

## Critical Finding: Cleaning Doesn't Work

**The GC field cannot be deconvolved at native 26 kHz resolution without self-calibration.**

Test results from subband_030 (1156-1166 MHz):
- MFS continuum: peak=138 mJy, RMS=52.9 mJy (sidelobe-dominated), SNR=2.6
- Per-channel: peaks=268-276 mJy, RMS=94.8 mJy, 3σ threshold=284 mJy
- All masks (5σ, 3σ, absolute) are empty -- sidelobes dominate everywhere
- All 3 cleaning strategies (mask+threshold, mask+threshold, no-mask+threshold) → modelFlux=0

**Decision needed**: Proceed with dirty images for sideband analysis (sidelobes are frequency-smooth, will subtract), or invest in self-calibration first.

## Gotchas / Lessons Learned

1. **Never run parallel jobs reading the monolithic MS** -- 2-year blocker
2. **exportfits channel keyword doesn't exist** -- use imsubimage+exportfits
3. **CRVAL3 in FITS headers** is cube reference freq, not per-channel; read from filenames
4. **Band-edge subbands** (0, 85) and **RFI subbands** (~8) take 2-4x longer
5. **auto-multithresh doesn't work** for per-channel imaging -- SNR too low AND sidelobes dominate
6. **Cleaning fundamentally limited** by dynamic range of the GC field without self-cal
7. **Per-channel tclean is slow** (~7 min/call) because it scans the full 383-channel MS for each channel; cube mode is much faster
8. **imstat RMS** over the full image is dominated by sidelobes, not thermal noise
9. **AxionLinePop** cloned to `reference/AxionLinePop/` -- private repo, CARL/CASCAL SBI with TensorFlow
10. **pandas not available** in casa_cookbook conda env -- scripts use csv/numpy instead
11. **Sam's data** only exists as Slack uploads (may have expired WeTransfer links)
