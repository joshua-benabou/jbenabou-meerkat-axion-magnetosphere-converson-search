# Master Plan: MeerKAT Galactic Center Data Reduction for Axion Search

## Project Overview

**Goal**: Reduce a 17 TB MeerKAT L-band measurement set (MS) of the Galactic Center into ~30,000 individual FITS images (one per spectral channel at native ~26 kHz resolution) covering 800-1400 MHz. These will be used to search for spectral line signatures of axion-photon conversion in neutron star magnetospheres.

**Dataset**: MeerKAT L-band, ~30,000 spectral channels, stored as a single monolithic measurement set on the Lawrencium cluster at `/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/SCI-20210212-SS-01/17TBdataset/`.

**Prior work**: ~81 channels at 1.3 GHz have been processed (both niter=0 and niter=10000 cleaned versions) as proof of concept. Stored as `SgrA*_niter_X_chan_Y.image.fits`.

---

## The Core Problem

CASA (tclean, split, etc.) places read locks on the monolithic MS, making parallel processing impossible. Every attempt to run multiple SLURM array jobs in parallel (whether tclean or split) results in serial execution because jobs block each other waiting for locks. This has been confirmed by both NRAO (Urvashi Rao Venkata) and Aritra Basu (Tautenburg).

Additional findings from Srikrishna Sekhar (NRAO, co-author of the IDIA frocc pipeline):
- tclean's `savemodel` parameter sets a **write lock** — must be set to `"none"`
- Using `specmode='mfs'` with `selectdata=True, spw='0:i_chan'` avoids the cube gridder overhead
- The frocc pipeline is SLURM-based and could potentially be adapted to run on Lawrencium
- Both suggestions were tried but the parallelism issue persisted, suggesting the core problem is MS-level read locking rather than tclean configuration

---

## Pipeline Overview

The pipeline has 7 phases. Phases 0-1 and 3a are complete. Phase 2 (continuum subtraction) has been **dropped** -- the sideband analysis approach makes it unnecessary.

```
Phase 0: Reconnaissance & Setup                              [COMPLETE]
    |
Phase 1: Subband Splitting (SERIAL)                          [COMPLETE]
    |
Phase 2: Continuum Subtraction                               [DROPPED - not needed]
    |
Phase 3a: Dirty Imaging (PARALLEL across subbands)           [COMPLETE]
    |
Phase 3b: Cleaned Imaging (PARALLEL across subbands)         [NEEDS REDO - auto-multithresh failed]
    |
Phase 4: FITS Export & Assembly                              [COMPLETE for dirty; redo after 3b]
    |
Phase 5: RFI Flagging / Channel Quality Map                  [NEW - mark bad channels]
    |
Phase 6: Sanity Checks & Signal Injection Tests              [NEW]
    |
Phase 7: Axion Search via Sideband Analysis                  [NEW]
```

---

## Phase 0: Reconnaissance & Setup

**Goal**: Understand the MS metadata, verify the data, and set up the software environment.

### Tasks

1. **Inspect the MS with `listobs` and `msmd`**:
   - Confirm total number of channels, channel width, frequency range
   - Determine current reference frame (should be TOPO)
   - Check number of spectral windows (SPWs)
   - Get integration time, number of antennas, polarizations
   - Identify any existing flagging or calibration applied

2. **Install software**:
   - CASA 6.x (modular casatools + casatasks via pip, NOT monolithic CASA)
   - WSClean (latest version, compiled with support for MeerKAT)
   - Tricolour (for additional RFI flagging if needed): https://github.com/ratt-ru/tricolour
   - Python environment with astropy, numpy, etc.

3. **Determine if data is already calibrated**:
   - MeerKAT archive data may come with calibration applied. Check the DATA vs CORRECTED_DATA columns.
   - If not calibrated, calibration must be done first (out of scope of this plan; would need a separate effort with a radio astronomer).

4. **Small test extraction**:
   - Extract a single subband (~10-20 MHz, ~500 channels) from the MS
   - Time how long this takes
   - Use this test subband for all subsequent phase testing

### Key commands
```python
import casatools
msmd = casatools.msmetadata()
msmd.open('path/to/ms')
print(msmd.nchans(0))        # channels per SPW
print(msmd.chanfreqs(0))     # frequency array
print(msmd.reffreq(0))       # reference frequency
print(msmd.exposuretime())   # integration time
msmd.close()
```

```python
from casatasks import listobs
listobs(vis='path/to/ms', listfile='listobs_output.txt')
```

---

## Phase 1: TOPO to LSRK Conversion + Subband Splitting

**Goal**: Convert the frequency reference frame from TOPO to LSRK and split the monolithic MS into ~30-60 independent subband MS files (each 10-20 MHz wide at native channel resolution).

**Why this must be serial**: Any parallel access to the monolithic MS causes read locks. This phase is the bottleneck but only needs to be done once.

### Approach A: `mstransform` (Recommended first attempt)

`mstransform` can do the TOPO->LSRK regridding AND output to a new MS simultaneously. Process one subband at a time.

```python
from casatasks import mstransform

# Example: extract subband starting at 1.0 GHz, width 10 MHz
mstransform(
    vis='path/to/original.ms',
    outputvis='subbands/subband_1000MHz.ms',
    regridms=True,
    outframe='LSRK',
    mode='frequency',
    start='1.0GHz',
    width='26.123kHz',    # native channel width
    nchan=383,            # ~10 MHz / 26.123 kHz
    datacolumn='data',    # or 'corrected' if calibrated
    spw='0'               # adjust based on listobs output
)
```

Write a loop script that processes subbands sequentially:
```bash
# submit_subband_split.sh
# Run ONE subband at a time to avoid MS locking
for i in $(seq 0 $N_SUBBANDS); do
    python split_subband.py --subband_index $i
done
```

**Estimated time**: If each subband split takes ~1-2 hours (based on prior experience with single-channel split taking ~1h40m, but subbands are larger reads with fewer seeks), then ~30-60 subbands = ~30-120 hours of serial processing. This is a one-time cost.

**Disk space**: The subbands together will be comparable in size to the original MS (~17 TB). Ensure sufficient scratch space (~35 TB total during this phase).

### Approach B: casatools direct conversion (Aritra's suggestion)

If `mstransform` is too slow or problematic, use the lower-level casatools approach:
https://gist.github.com/tnakazato/af0496d4c3a25c828f240579f4e3c050

This gives more control but requires more manual implementation.

### Approach C: Skip LSRK conversion, use WSClean directly

WSClean may be able to handle TOPO data directly for narrow subbands where the TOPO-LSRK shift is small (a few channels at most). This should be tested. If viable, Phase 1 simplifies to just splitting the MS into subbands without regridding.

### Validation
- Check that the output subband MS has the expected number of channels
- Verify the reference frame is LSRK (using msmd)
- Compare a quick dirty image of one channel from the subband vs the original MS to confirm data integrity

---

## Phase 2: Continuum Subtraction -- DROPPED

**Status**: Not needed. The sideband analysis approach (Phase 7) removes the need for visibility-domain continuum subtraction. For each axion frequency, the background model is derived from nearby frequency channels (the "sideband"), which naturally captures and removes the continuum spatial morphology. This is more robust than uvsub for our science case.

---

## Phase 3: Channel Imaging

### Phase 3a: Dirty Imaging -- COMPLETE

**Status**: All 86 subbands imaged with niter=0. 32,768 channel FITS files produced.
- Located in `images/subband_XXX/chan_NNNN_FREQ.fits`
- 512x512 pixels, 2 arcsec/pixel, Briggs robust=0.5
- Used CASA tclean (not WSClean) with cube mode, then imsubimage+exportfits per channel
- Timing: ~50 min per subband cube on lr7 (56 cores, 240 GB)

### Phase 3b: Cleaned Imaging -- NEEDS REDO

**Status**: First attempt with auto-multithresh FAILED. The noisethreshold=5.0 was too high for per-channel SNR (~2.7), so tclean produced zero clean components. Cleaned images are identical to dirty.

**Next approach**: Need a cleaning strategy that works at low per-channel SNR. Options:
1. **Continuum-derived clean mask**: Make a wideband continuum image (high SNR), create a mask from it, apply to per-channel tclean
2. **WSClean with auto-mask at lower threshold**: WSClean may handle this better than CASA's auto-multithresh
3. **Interactive mask from dirty continuum**: Use the existing dirty images averaged in frequency to define a mask

**The axion search requires CLEANED images** (per collaborator meeting 2026-04-13). The cleaning step must be validated via signal injection tests (Phase 6) before production.

### Original WSClean approach (for reference)

WSClean does NOT place file locks on the MS, is faster than tclean, and outputs FITS natively. It remains a viable alternative if CASA tclean cleaning continues to be problematic.

```bash
wsclean -name subband_XX_channels \
    -size 512 512 \
    -scale 2asec \
    -channels-out 383 \
    -niter 10000 \
    -auto-threshold 3 \
    -auto-mask 5 \
    -j 16 \
    -mem 80 \
    -no-mf-weighting \
    subband_XX.ms
```

### Important considerations

- **Image size**: Currently using 512x512. Full MeerKAT FOV at 2" would be ~4096x4096 but not needed for axion search if targets are near field center.
- **Memory**: WSClean with large cubes can require 1 TB+ RAM. 512x512 is manageable.
- **Single-node constraint**: WSClean is not optimized for multi-node MPI. Each SLURM job should request a single node.

---

## Phase 4: FITS Export & Assembly

**Goal**: Collect all channel images into a usable format for the axion analysis pipeline.

### WSClean output format
WSClean outputs FITS files by default. Each channel image will be named like:
```
subband_XX-0000-image.fits
subband_XX-0001-image.fits
...
```

### Tasks
1. Rename/organize FITS files with a consistent naming scheme that maps to absolute frequency
2. Verify FITS headers contain correct frequency information
3. Create a master catalog mapping filename -> frequency (MHz)
4. Optionally combine into a single FITS cube (may not be necessary for analysis)

---

## Resource Estimates

| Phase | Time (estimated) | Disk Space | Parallelizable |
|-------|-----------------|------------|----------------|
| 0: Recon | 1-2 hours | Minimal | N/A |
| 1: Split | 30-120 hours (serial) | ~17 TB (subbands) | No |
| 2: Continuum sub | 2-4 hours per subband | ~1 TB (models) | Yes |
| 3: Imaging | 4-12 hours per subband | ~2-5 TB (FITS) | Yes |
| 4: Assembly | 1-2 hours | Minimal | N/A |

**Total wall time (optimistic)**: ~1-2 weeks with good parallelization of Phases 2-3.
**Total disk space needed**: ~35-40 TB peak during Phase 1, settling to ~20-25 TB.

---

## Risk Factors & Fallbacks

1. **mstransform is also slow/locks**: If serial subband splitting is too slow, try:
   - Splitting without LSRK conversion (just frequency selection), then handle frame conversion in WSClean or post-processing
   - Using `taql` (Table Query Language) for raw MS splitting, which may bypass CASA's locking

2. **WSClean can't handle TOPO data**: If LSRK conversion is mandatory before WSClean, Phase 1 must include regridding.

3. **Memory limits**: If WSClean needs >256 GB per subband, reduce image size or number of channels per WSClean call.

4. **Calibration issues**: If the archival data is not well-calibrated, this will show up as imaging artifacts. Self-calibration in Phase 2 can help, but severe issues may require a full re-calibration.

5. **RFI**: MeerKAT L-band has known RFI bands. Channels heavily affected by RFI should be flagged before imaging. Tricolour can help with this.

---

## Phase 5: RFI Flagging / Channel Quality Map

**Goal**: Identify and mark channels affected by RFI. **Do NOT modify the existing dirty or cleaned FITS files.** Instead, produce a boolean channel quality map that downstream analysis can use to exclude bad channels.

### Approach

1. **Known MeerKAT L-band RFI**: GSM-900 (~925-960 MHz, subbands ~7-10), GPS L1 (1575.42 MHz, subband ~72), and others documented in MeerKAT RFI environment reports.
2. **Statistical detection**: For each channel FITS, compute image statistics (RMS, peak, kurtosis). Channels with anomalous statistics (e.g., RMS >> median RMS of neighboring channels) are flagged.
3. **Output**: A CSV or FITS table mapping `(subband, channel, frequency_MHz) -> rfi_flag (bool)` plus summary statistics.
4. **Preserve data**: The dirty images in `images/subband_XXX/` must remain untouched. The RFI map is metadata only.

### Known problematic subbands
- Subband 0 (856 MHz): band edge, channel 0 artifact
- Subbands ~7-10 (~920-960 MHz): GSM-900 mobile RFI
- Subband ~72 (~1575 MHz): GPS L1

---

## Phase 6: Sanity Checks & Signal Injection Tests

**Goal**: Validate the imaging and cleaning pipeline before trusting it for the axion search. This phase is done in Jupyter notebooks for interactive exploration.

### 6a: Point Source Injection Baseline

1. Take a dirty image (any clean channel with no RFI)
2. Add a synthetic point source (known flux, known position) to the dirty image
3. Clean the modified dirty image
4. Verify the cleaned version recovers the injected source with correct flux and position
5. **Note**: Axion signal point sources will be **resolved by multiple pixels** (not unresolved), so also test with extended Gaussian sources matching expected axion signal morphology

### 6b: Cleaning Fidelity Test (5-sigma threshold)

1. In cleaned images, measure the noise (RMS) and establish the 5-sigma detection threshold
2. Inject a signal at ~5-sigma into the **dirty** data (visibilities or image-plane)
3. Re-clean the data containing the injected signal
4. Verify that the cleaning process does **not** absorb, suppress, or distort the injected signal
5. This confirms that the cleaning step is safe for the axion search -- i.e., a real signal above threshold would survive cleaning

### Implementation
- Jupyter notebook(s) in the project root
- Use a few representative subbands (low-RFI, mid-band channels)
- Document results with plots showing before/after injection and cleaning

---

## Phase 7: Axion Search via Sideband Analysis

**Goal**: Search for axion-photon conversion signals from neutron star magnetospheres using cleaned channel images.

### Source Catalog
- **Sam Witte has a template bank of converting neutron stars** -- these are the target positions and expected signal properties
- Each target NS has a predicted axion signal frequency (or frequency range) and spatial extent
- Axion signals are **resolved by multiple pixels** in MeerKAT images

### Sideband Background Model

For each candidate axion frequency channel:

1. Define a **sideband** -- a set of nearby frequency channels (excluding the signal channel and a small guard band)
2. Use the spatial morphology of the sideband channels as the **background model**
   - The sideband naturally captures continuum emission, diffuse structure, and any frequency-smooth backgrounds
   - This is why continuum subtraction (Phase 2) is unnecessary
3. Subtract (or divide) the background model from the signal channel
4. Search the residual for excess emission at the expected NS position and morphology

### Detection Pipeline

1. For each NS in Sam's template bank:
   - Identify the predicted axion frequency → map to the corresponding channel FITS
   - Construct sideband from neighboring channels
   - Build background model from sideband spatial morphology
   - Compute residual at the NS position
   - Assess significance against the local noise
2. Catalog all candidates above threshold
3. Cross-check with RFI flag map (Phase 5) to reject false positives from RFI channels

### Validation (ties back to Phase 6)
- The signal injection test (Phase 6b) must confirm that the cleaning + sideband subtraction pipeline preserves injected signals above the 5-sigma threshold

---

## References

- WSClean documentation: https://wsclean.readthedocs.io/en/latest/
- WSClean image cubes: https://wsclean.readthedocs.io/en/latest/making_image_cubes.html
- CASA mstransform: https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.manipulation.mstransform.html
- casatools TOPO->LSRK gist: https://gist.github.com/tnakazato/af0496d4c3a25c828f240579f4e3c050
- IDIA MeerKAT pipeline (frocc): https://github.com/idia-astro/frocc
- MALS survey pipeline: https://mals.iucaa.in/ and https://iopscience.iop.org/article/10.3847/1538-4357/abcb85/pdf
- Tricolour RFI flagger: https://github.com/ratt-ru/tricolour
- MeerKAT IDIA pipelines: https://idia-pipelines.github.io/
