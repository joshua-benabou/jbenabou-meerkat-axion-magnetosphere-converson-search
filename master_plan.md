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

The pipeline has 5 phases. Each phase must be validated on a small test before scaling up.

```
Phase 0: Reconnaissance & Setup
    |
Phase 1: TOPO -> LSRK Conversion + Subband Splitting (SERIAL)
    |
Phase 2: Continuum Subtraction (PARALLEL across subbands)
    |
Phase 3: Channel Imaging with WSClean (PARALLEL across subbands)
    |
Phase 4: FITS Export & Assembly
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

## Phase 2: Continuum Subtraction

**Goal**: Remove the bright continuum emission (especially from the Galactic Center) from the visibilities. This is critical because:
1. The GC has extremely high sky temperature at L-band
2. Continuum subtraction reduces spectral noise
3. It prevents continuum sidelobes from contaminating channel maps

### Approach

For each subband:

1. **Make a continuum image** of the subband (averaging all channels):
   ```bash
   wsclean -name subband_XX_cont \
       -size 4096 4096 \
       -scale 2asec \
       -niter 50000 \
       -auto-threshold 3 \
       -auto-mask 5 \
       -channels-out 1 \
       -j 8 \
       -mem 80 \
       subband_XX.ms
   ```

2. **Predict the continuum model back into the MS** and subtract:
   ```bash
   wsclean -name subband_XX_cont \
       -size 4096 4096 \
       -scale 2asec \
       -predict \
       subband_XX.ms
   ```
   Then use CASA or taql to subtract:
   ```python
   # In CASA
   uvsub(vis='subband_XX.ms')
   ```

3. **Optional: Self-calibration** on the continuum to improve dynamic range before subtraction. This is important for the GC field which has very bright sources.

### Parallelization
- Each subband is an independent MS file, so all subbands can be processed in parallel via SLURM array jobs.
- Allocate one compute node per subband.

### Validation
- Image a few channels before and after continuum subtraction
- Verify that bright point sources are removed
- Check that noise levels are reduced as expected

---

## Phase 3: Channel Imaging with WSClean

**Goal**: Produce individual channel images for every spectral channel in every subband.

### Why WSClean over tclean
- WSClean does NOT place file locks on the MS (confirmed by Aritra)
- WSClean is significantly faster than tclean for this type of work
- WSClean can produce image cubes or individual channel maps natively

### Approach

For each subband (parallelized across subbands via SLURM):

```bash
wsclean -name subband_XX_channels \
    -size 4096 4096 \
    -scale 2asec \
    -channels-out 383 \          # one output per channel
    -niter 10000 \
    -auto-threshold 3 \
    -auto-mask 5 \
    -j 8 \
    -mem 80 \
    -no-mf-weighting \
    subband_XX.ms
```

### Important considerations

- **Image size**: Aritra suggests starting with smaller images (512x512 or even 128x128) for testing, then scaling to the full FOV. The full MeerKAT L-band primary beam is large (~1 degree at 1.4 GHz), requiring ~4096x4096 at 2" pixels.
- **Memory**: WSClean with large cubes can require 1 TB+ RAM (per Aritra). For 4096x4096 images with ~500 channels per subband, memory may be the limiting factor.
  - Mitigation: Process in smaller frequency chunks within each subband, or use smaller image sizes and mosaic.
- **Single-node constraint**: WSClean is not optimized for multi-node parallelism. Each SLURM job should request a single node with many cores and large RAM.
- **Phase-center shifting**: If the full FOV is too large for memory, image smaller patches (e.g., 1024x1024) at different phase centers and stitch. For the axion search, it may be sufficient to image only the central region around Sgr A* and known magnetars.

### SLURM job template
```bash
#!/bin/bash
#SBATCH --job-name=wsclean_subband
#SBATCH --array=0-59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G          # adjust based on testing
#SBATCH --time=24:00:00

SUBBAND_IDX=$SLURM_ARRAY_TASK_ID
INPUT_MS="subbands/subband_${SUBBAND_IDX}.ms"
OUTPUT_PREFIX="images/subband_${SUBBAND_IDX}"

wsclean -name ${OUTPUT_PREFIX} \
    -size 4096 4096 \
    -scale 2asec \
    -channels-out 383 \
    -niter 10000 \
    -auto-threshold 3 \
    -auto-mask 5 \
    -j 16 \
    -mem 80 \
    ${INPUT_MS}
```

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

## References

- WSClean documentation: https://wsclean.readthedocs.io/en/latest/
- WSClean image cubes: https://wsclean.readthedocs.io/en/latest/making_image_cubes.html
- CASA mstransform: https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.manipulation.mstransform.html
- casatools TOPO->LSRK gist: https://gist.github.com/tnakazato/af0496d4c3a25c828f240579f4e3c050
- IDIA MeerKAT pipeline (frocc): https://github.com/idia-astro/frocc
- MALS survey pipeline: https://mals.iucaa.in/ and https://iopscience.iop.org/article/10.3847/1538-4357/abcb85/pdf
- Tricolour RFI flagger: https://github.com/ratt-ru/tricolour
- MeerKAT IDIA pipelines: https://idia-pipelines.github.io/
