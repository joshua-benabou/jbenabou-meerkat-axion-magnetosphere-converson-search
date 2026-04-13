# CLAUDE.md - MeerKAT Data Reduction Project

## CRITICAL DATA SAFETY WARNING

**THE 17 TB MEASUREMENT SET TOOK 1 MONTH TO TRANSFER FROM SOUTH AFRICA. IT IS IRREPLACEABLE.**

**AT ALL COSTS, DO NOT ACCIDENTALLY DELETE OR CORRUPT THE DATA.**

Before running ANY command that touches the MS or its parent directories:
1. **Triple-check** that the command is READ-ONLY or writes to a SEPARATE output path
2. **NEVER** use `rm`, `mv`, or any destructive command on or near the MS directory
3. **NEVER** run a command that could modify the MS in-place (e.g., `flagdata` without a backup, `split` with `outputvis` pointing inside the MS)
4. **NEVER** use wildcard deletions (`rm -rf *`) anywhere in the data tree
5. **Always verify paths** before executing — a typo could destroy months of work
6. When in doubt, **ASK THE USER** before proceeding

The MS path is:
`/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/SCI-20210212-SS-01/17TBdataset/scratch/kat/1621534878_20231022T19_11_29/1621534878_sdp_l0.ms`

The compressed archive is also preserved at:
`/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/SCI-20210212-SS-01/17TBdataset/1621534878_sdp_l0.ms.tar.gz`

---

## Project Context

You are helping Josh Benabou (UC Berkeley, Safdi group) reduce a 17 TB MeerKAT L-band measurement set of the Galactic Center. The ultimate science goal is searching for axion-photon conversion signals from neutron star magnetospheres across 800-1400 MHz at native ~26 kHz channel resolution (~30,000 channels).

This project has been stalled for ~2 years because of technical difficulties with CASA parallelization. See `master_plan.md` for the full pipeline design, and `reference/` for email exchanges and Slack history that document the problem and expert advice received.

## What Has Been Tried (and Failed)

1. **tclean with SLURM array jobs**: Each job processes one channel. Jobs interfere via MS read locks, running effectively serially. Tried both `specmode='cube'` and `specmode='mfs'` with `selectdata`.
2. **CASA split in parallel**: Same locking problem as tclean.
3. **savemodel='none'** was already set; did not fix the issue.
4. **specmode='mfs' with spw selection** (Srikrishna's suggestion): Same cross-talk / slowdown.

## What Has NOT Been Tried Yet

These are the remaining suggestions from radio astronomers (Aritra Basu, Urvashi Rao Venkata, Srikrishna Sekhar):

1. **Serial subband splitting with `mstransform`**: Split the monolithic MS into ~30-60 subbands (10-20 MHz each) ONE AT A TIME (serial). This avoids locking because only one process touches the MS. Once split, downstream work can be fully parallel.
2. **TOPO to LSRK conversion via casatools** (not tclean): https://gist.github.com/tnakazato/af0496d4c3a25c828f240579f4e3c050
3. **WSClean instead of tclean**: No file.lock issue. Must run on single node (not multi-node MPI). Allocate via `-j` and `-mem` flags.
4. **Continuum subtraction**: Make ~16 subband continuum images (~50 MHz width), self-cal, subtract from visibilities before channel imaging.
5. **Smaller image cubes**: Start with 128x128 or 512x512 pixels instead of 4096x4096 for testing.
6. **Adapting the frocc pipeline** (https://github.com/idia-astro/frocc) to run on Lawrencium. It's SLURM-based, written for Ilifu but potentially portable.
7. **MALS pipeline philosophy** (https://mals.iucaa.in/) for blind spectral line searches.

## The Dataset

- **Location on Lawrencium (extracted MS)**: `/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/SCI-20210212-SS-01/17TBdataset/scratch/kat/1621534878_20231022T19_11_29/1621534878_sdp_l0.ms`
- **Compressed archive**: `/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/SCI-20210212-SS-01/17TBdataset/1621534878_sdp_l0.ms.tar.gz` (12 TB compressed)
- **Size**: ~17 TB monolithic measurement set
- **Band**: L-band, 800-1400 MHz
- **Channels**: ~30,000 at ~26.123 kHz native resolution
- **Target**: Sgr A* / Galactic Center
- **Proposal ID**: SCI-20210212-SS-01 (Galactic Centre Legacy Survey)
- **Prior processed data**: ~81 channels at 1.3 GHz in `scratch/` directory, files named `SgrA*_niter_X_chan_Y.image.fits`

## Cluster Environment

- **Cluster**: Lawrencium (UC Berkeley HPC)
- **Scheduler**: SLURM
- **Project allocation**: `pc_heptheory`
- **Scratch path**: `/global/scratch/projects/pc_heptheory/jbenabou/`

## How to Proceed

Follow `master_plan.md` phases 0-4 in order. Key principles:

1. **Always test small first**: Before any large operation, test on a tiny subset (1 subband, 64x64 pixels, ~10 channels).
2. **Serial access to the monolithic MS**: NEVER run parallel jobs that all read from the original 17 TB MS simultaneously. This is the #1 lesson from 2 years of debugging.
3. **Parallel access to split subbands is fine**: Once the MS is split into independent subband MS files, parallelism works normally.
4. **Prefer WSClean over tclean for imaging**: WSClean is faster, doesn't lock, and outputs FITS natively.
5. **Monitor disk space**: The cluster scratch space may be limited. Phase 1 temporarily doubles storage requirements (~35 TB).

## Software Requirements

- CASA 6.x (modular: `pip install casatools casatasks`)
- WSClean (compile from source or load via module if available)
- Tricolour for RFI flagging: `pip install tricolour` or https://github.com/ratt-ru/tricolour
- Python 3.8+ with astropy, numpy

### Conda Environment

The existing conda environment with CASA installed is from collaborator Orion (osning):
```
/global/home/users/osning/.conda/envs/casa_cookbook
```
This was used in all prior SLURM scripts. It may or may not still work after cluster updates — needs to be tested.

Check what modules are available on Lawrencium:
```bash
module avail casa
module avail wsclean
```

## Existing Code

Josh has existing Python and SLURM scripts from prior attempts. Look in the dataset directory and nearby directories for:
- `MS_to_FITS_1chan_mfs.py` - tclean script for single channel
- `run_casa_MS_to_FITS_mfs_1chan.sh` - SLURM submission script
- `split_MS.py` - CASA split script
- `run_split_MS.sh` - SLURM submission for split

These scripts document what was tried. They can be adapted for the new pipeline but should NOT be run as-is (they hit the locking problem).

## Key Contacts (for context only)

- **Aritra Basu** (Tautenburg): Radio astronomer, wrote MeerKAT Galactic Plane Survey pipeline. Gave detailed advice on WSClean, subband splitting, continuum subtraction.
- **Urvashi Rao Venkata** (NRAO): CASA developer. Confirmed MS locking issue, suggested mstransform + split approach.
- **Srikrishna Sekhar** (NRAO): Co-author of frocc/IDIA pipeline. Suggested savemodel='none', specmode='mfs' approach.
- **Sam Witte** (ICC Barcelona): Collaborator, connected Josh with Aritra.
- **Joshua Foster** (MIT): Collaborator, working on SBI analysis pipeline.
- **Itay Bloch**: Collaborator, working on ML analysis.
- **Ben Safdi** (Berkeley): PI.

## Science Requirement: Native Frequency Resolution

**The axion search requires one FITS image per spectral channel at the native 26.123 kHz resolution (~32,768 channels across 800-1400 MHz). Do NOT downbin in frequency.** The axion line width may be comparable to the channel width, so preserving native resolution is essential. Every channel gets its own dirty and cleaned image.

## Performance & Parallelization

**Always maximize use of cluster resources.** Lawrencium has 410 nodes with 32 cores and ~93 GB RAM each.

- **Use SLURM for all non-trivial compute** — never run long jobs interactively on login/shared nodes. Submit via `sbatch`.
- **Parallelize across nodes**: Use SLURM array jobs (`--array=0-85`) for independent tasks. Each subband is independent after Phase 1.
- **Parallelize within nodes**: Use Python `multiprocessing` to utilize all 32 cores per node (e.g., imaging multiple channels simultaneously).
- **I/O bottleneck awareness**: The monolithic MS is striped across 8 Lustre HDD OSTs. More than ~8 simultaneous readers of the original MS will saturate I/O. But split subbands are independent files — no I/O contention.
- **LSRK regridding is unnecessary**: TOPO->LSRK shift is <0.1 channel width within a single observation. Use `split` (not `mstransform` with `regridms=True`) for Phase 1. Compute LSRK frequencies as metadata.
- **Test small first, then scale**: Run 1-channel and 1-subband tests before full 32,768-channel production runs. Use SLURM for timing tests too.

## Important Notes

- The data may already be calibrated (from the archive). Check DATA vs CORRECTED_DATA columns.
- Pixel size is 2 arcsec (between 1/3 and 1/5 of synthesized beam).
- FOV changes with frequency across the L-band.
- Units are Jy/beam.
- There are known RFI bands in MeerKAT L-band; channels around 1.3 GHz are relatively RFI-quiet.
- The 81 previously processed channels start at 1.3 GHz.

## Session Handoff

**At the end of every session**, create or update `handoff.md` in the project root with a summary of what happened during the session. This file serves as context for the next session so continuity is not lost. Include:

1. **What was done** — key actions, commands run, files created/modified
2. **Current state** — what's working, what's broken, where things stand
3. **Next steps** — what should be done next, any blockers or open questions
4. **Gotchas / lessons learned** — anything surprising that came up

At the **start of every session**, read `handoff.md` before doing anything else to pick up where the last session left off.
