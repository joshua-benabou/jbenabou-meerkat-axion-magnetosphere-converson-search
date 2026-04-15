# Session Handoff

**Last updated**: 2026-04-15

## Current Pipeline State

### Phase 1: Subband Splitting -- COMPLETE
- All 86 subbands split from the monolithic 17 TB MS
- Located in `subbands/subband_XXX.ms`, 383 channels each
- Total: ~12 TB

### Phase 3a: Dirty Imaging -- COMPLETE
- All 86 subbands, 32,768 total channel FITS files
- Located in `images/subband_XXX/chan_NNNN_FREQ.fits`
- 512x512 pixels, 2 arcsec/pixel, Briggs robust=0.5

### Phase 3b: Cleaned Imaging -- PRODUCTION RUNNING (job 21933772)
- **Fixed**: Previous auto-multithresh approach failed (used wrong RMS). Now uses 10% of peak threshold.
- **niter=100000**, threshold=10% of peak (adaptive per subband), no mask
- Validated: per-channel model flux matches single-channel debug test (15 Jy/chan)
- Wall time: ~115 min per subband, ~165 node-hours total for all 86
- Output: `images/subband_XXX/cleaned/chan_NNNN_FREQ.fits`
- Old broken cleaned FITS are automatically removed by the script

### Phase 5: RFI Flagging -- COMPLETE
- 1,656/32,768 channels flagged (5.1%)
- Output: `rfi_channel_flags.csv`, `plots/rfi_overview.png`

### Phase 7: Axion Search -- SKELETON WRITTEN
- `scripts/phase7_axion_search.py` with sideband subtraction code
- Waiting on integration with Sam's NS template bank

## Session 2026-04-15: What Was Done

### Cleaning Fix
1. **Debug cleaning test** (`debug_cleaning.py`) had already run -- showed cleaning WORKS with correct threshold
   - Root cause of failure: full-image RMS (95 mJy, sidelobe-dominated) used instead of corner RMS (40 mJy)
   - Fix: threshold = 10% of peak flux, no mask
   - 1%, 5%, 10% peak all converge to same solution
2. **Generated comparison plots** with LaTeX/plotting defaults: `plots/cleaning_comparison_sb030_ch100.png`, `plots/cleaning_strategies_comparison.png`
3. **niter investigation**: In cube mode, niter is TOTAL across all channels, not per-channel
   - niter=10,000 on 383 channels = only 26 iter/channel (not enough)
   - niter=100,000 on 383 channels = ~261 iter/channel (model flux matches single-channel test)
   - niter test results (20 channels): 10k/100k/1M all tested
4. **Full subband validation** (niter=100000, subband 30): 115 min, model flux 5,782 Jy (15.1 Jy/chan -- matches debug)
5. **Production launch**: job 21933772, array 0-85, ~165 node-hours total

### Sam Witte's NS Population Code
1. **Unzipped** `Real_Analysis.zip` (2 GB) to `reference/Real_Analysis/`
2. **Explored code structure**:
   - `Create_Population.py`: generates NS populations (B, P, chi, position, age, etc.)
   - `Analyze_Population.py`: compiles per-NS ray-tracing into `Combined_Flux.dat`
   - `Mix_Match.py`: resamples populations (new positions/Doppler) for fast Monte Carlo
3. **Available data**:
   - Young pop (tau_ohm=10 Myr): 10 realizations, 3 masses (3.54, 4.13, 7.08 µeV)
   - Old pop (tau_ohm=1 Tyr): 4 realizations, 2 masses (3.54, 4.13 µeV)
4. **Key findings**:
   - Old pop dominates: 0.77 Jy vs 0.044 Jy total flux at g=1e-12
   - 45 NSs contribute 50% of flux, 264 contribute 90%
   - Signal spectrum shows sharp lines near nu=ma*c^2/h with Doppler broadening
   - Frequency range 856-1712 MHz matches exactly ma=3.54-7.08 µeV
5. **Generated 5 plots** with LaTeX defaults in `plots/sam_*.png`
6. **Created notebook** `notebooks/explore_sam_populations.ipynb` (not executed via nbconvert due to kernel issue -- run manually in Jupyter)

### PDF & Plots
- Regenerated all plots with LaTeX + Josh's `plotting_defaults.py`
- Updated `processing_summary.md` with cleaning results, Sam's population section
- Rebuilt PDF: `processing_summary.pdf` (6.4 MB, all plots embedded)
- Removed 128 stale CASA log files + 10 TempLattice dirs from data directory

## Key Files Created/Modified

- `scripts/phase3_clean_cube.py` -- updated: niter=100000, threshold=10% peak, no mask
- `scripts/phase3_clean_submit.sh` -- updated comments
- `scripts/debug_cleaning.py` / `debug_cleaning_submit.sh` -- cleaning threshold debug
- `scripts/test_niter.py` / `test_niter_submit.sh` -- niter scaling test
- `scripts/test_full_subband_clean.py` / `test_full_subband_submit.sh` -- full subband validation
- `scripts/plot_cleaning_comparison.py` -- dirty vs cleaned comparison plots
- `scripts/plot_sam_populations.py` -- Sam's NS population plots
- `notebooks/explore_sam_populations.ipynb` -- interactive exploration of Sam's data
- `reference/Real_Analysis/` -- Sam Witte's code and pre-computed populations
- `processing_summary.md` / `processing_summary.pdf` -- updated

## Next Steps

1. **Monitor cleaning jobs** (21933772): Should complete in ~2 hours as nodes become available
2. **Validate cleaned images**: Compare dirty vs cleaned across multiple subbands after jobs finish
3. **Run the Sam exploration notebook** in Jupyter to interact with plots
4. **Integrate Sam's Mix_Match output** into Phase 7 axion search pipeline
5. **Ask Sam** to generate more axion masses (finer grid between 3.54 and 7.08 µeV)
6. **Phase 6 sanity checks**: Run `phase6_sanity_checks.ipynb`

## Gotchas / Lessons Learned (cumulative)

1. **Never run parallel jobs reading the monolithic MS** -- 2-year blocker
2. **exportfits channel keyword doesn't exist** -- use imsubimage+exportfits
3. **auto-multithresh doesn't work** at per-channel SNR -- use percentage-of-peak threshold
4. **Full-image RMS != thermal noise** for GC field -- measure in off-source corners
5. **tclean cube niter is TOTAL across all channels** -- need niter >> nchan for proper cleaning
6. **Cleaning convergence** confirmed by model flux matching, not just RMS metrics
7. **pandas not available** in casa_cookbook -- use csv/numpy
8. **LaTeX on Lawrencium**: set PATH to texlive/2016 before matplotlib import
9. **Plotting defaults**: import from `/global/scratch/projects/pc_heptheory/jbenabou/plotting_defaults.py`
10. **PDF rebuild**: use xelatex (not pdflatex), avoid unicode chars in markdown
11. **Sam's old pop** dominates signal by ~17x over young pop
