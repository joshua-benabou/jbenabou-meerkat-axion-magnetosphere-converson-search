# Session Handoff

**Last updated**: 2026-04-13

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

### Phase 3b: Cleaned Imaging -- COMPLETED BUT INEFFECTIVE
- Ran tclean with niter=10000 and auto-multithresh on all 86 subbands
- 32,768 cleaned FITS in `images/subband_XXX/cleaned/chan_NNNN_FREQ.fits`
- **PROBLEM**: auto-multithresh with noisethreshold=5.0 found NO emission in any channel
  - Per-channel SNR is only ~2.7, well below the threshold
  - Cleaned images are identical to dirty images (zero clean components produced)
  - This approach fundamentally doesn't work for 26 kHz channel-width imaging

## Meeting Decisions (2026-04-13, Josh + collaborators)

Key decisions from the collaborator meeting:

1. **No continuum subtraction** -- the sideband analysis approach makes it unnecessary
2. **Axion analysis uses CLEANED images** -- cleaning must be redone properly
3. **Sideband analysis**: for each axion frequency, use nearby frequency channels' spatial morphology as background model
4. **Sam Witte has a template bank of converting neutron stars** -- these define the search targets
5. **Axion signals are resolved by multiple pixels** (not point sources)
6. **RFI flagging**: mark bad channels with a boolean flag, do NOT modify existing dirty images
7. **Sanity checks needed**:
   - Inject point source into dirty image, clean it, verify recovery (Jupyter notebook)
   - Assess 5-sigma threshold in cleaned images
   - Inject signal into dirty data, re-clean, verify cleaning doesn't destroy the signal

## What Needs to Happen Next

### Priority 1: Fix Cleaning (Phase 3b redo)
The auto-multithresh approach failed (noisethreshold=5.0 too high for per-channel SNR ~2.7). Need a working cleaning strategy:
- Continuum-derived clean mask
- WSClean with lower auto-mask threshold
- Must be validated via signal injection (Phase 6) before production

### Priority 2: Sanity Checks (Phase 6, Jupyter notebooks)
- Point source injection into dirty image → clean → verify recovery
- 5-sigma threshold assessment in cleaned images
- Signal injection into dirty data → re-clean → verify signal survives

### Priority 3: RFI Channel Quality Map (Phase 5)
- Mark known RFI bands (GSM-900, GPS L1, etc.)
- Statistical anomaly detection across channels
- Output boolean flag map -- don't touch existing images

### Priority 4: Sideband Analysis Pipeline (Phase 7)
- Get Sam's template bank of NS targets
- Implement sideband background subtraction
- Search for excess emission at predicted axion frequencies/positions

## Key Technical Details

- **Cluster**: Lawrencium, partition lr7 (56 cores, 240 GB) for imaging, lr6 for splitting
- **Conda env**: `/global/home/users/osning/.conda/envs/casa_cookbook`
- **Timing**: ~50 min per dirty cube, ~56 min per clean cube (on lr7)
- **Export bug fix**: CASA `exportfits(channel=N)` doesn't work; must use `imsubimage(chans=str(N))` then `exportfits()`

## Scripts

All in `scripts/`:
- `phase1_split_subbands.py` / `phase1_submit.sh` -- subband splitting
- `phase3_image_cube.py` / `phase3_cube_submit.sh` -- dirty imaging
- `phase3_clean_cube.py` / `phase3_clean_submit.sh` -- cleaned imaging
- `phase3_export_recovery.py` / `phase3_export_recovery_submit.sh` -- FITS export fix
- `validate_fits.py` -- FITS validation

## Gotchas / Lessons Learned

1. **Never run parallel jobs reading the monolithic MS** -- this was the original 2-year blocker
2. **exportfits channel keyword doesn't exist** in this CASA version -- use imsubimage+exportfits
3. **CRVAL3 in FITS headers** is the cube reference frequency, not the per-channel frequency; read freq from filenames
4. **Band-edge subbands** (0, 85) and **RFI subbands** (~8, GSM-900) take 2-4x longer to image
5. **auto-multithresh doesn't work** for per-channel imaging at native resolution (SNR too low)
