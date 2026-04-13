#!/usr/bin/env python3
"""Generate a PDF with paper-ready data processing description and embedded plots."""

import os
from fpdf import FPDF

PROJECT_DIR = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project'
)
PLOTS_DIR = os.path.join(PROJECT_DIR, 'plots')
PDF_PATH = os.path.join(PROJECT_DIR, 'processing_summary.pdf')


class Report(FPDF):
    def header(self):
        if self.page_no() > 1:
            self.set_font('Helvetica', 'I', 8)
            self.cell(0, 5, 'MeerKAT GC Data Reduction - Benabou et al.', align='L')
            self.cell(0, 5, f'Page {self.page_no()}', align='R', new_x='LMARGIN', new_y='NEXT')
            self.ln(2)

    def section(self, title):
        self.set_font('Helvetica', 'B', 12)
        self.cell(0, 8, title, new_x='LMARGIN', new_y='NEXT')
        self.ln(1)

    def body_text(self, text):
        self.set_font('Times', '', 10)
        self.multi_cell(0, 5, text)
        self.ln(2)

    def figure(self, img_path, caption, width=180):
        if os.path.exists(img_path):
            # Center the image
            x = (self.w - width) / 2
            # Check if we need a new page (image + caption ~ 100mm)
            if self.get_y() + 110 > self.h - 20:
                self.add_page()
            self.image(img_path, x=x, w=width)
            self.ln(2)
            self.set_font('Helvetica', 'I', 9)
            self.multi_cell(0, 4, caption)
            self.ln(4)
        else:
            self.set_font('Helvetica', 'I', 9)
            self.cell(0, 5, f'[Plot not found: {os.path.basename(img_path)}]',
                      new_x='LMARGIN', new_y='NEXT')
            self.ln(2)

    def table_row(self, label, value):
        self.set_font('Helvetica', 'B', 9)
        self.cell(55, 5, label)
        self.set_font('Times', '', 9)
        self.cell(0, 5, value, new_x='LMARGIN', new_y='NEXT')


pdf = Report()
pdf.set_auto_page_break(auto=True, margin=20)
pdf.add_page()

# === Title ===
pdf.set_font('Helvetica', 'B', 18)
pdf.cell(0, 12, 'MeerKAT Galactic Center', align='C', new_x='LMARGIN', new_y='NEXT')
pdf.cell(0, 12, 'Data Reduction Summary', align='C', new_x='LMARGIN', new_y='NEXT')
pdf.ln(3)
pdf.set_font('Helvetica', '', 11)
pdf.cell(0, 7, 'J. Benabou, UC Berkeley / Safdi Group', align='C', new_x='LMARGIN', new_y='NEXT')
pdf.cell(0, 7, 'April 2026', align='C', new_x='LMARGIN', new_y='NEXT')
pdf.ln(8)

# === 1. Observation ===
pdf.section('1. Observation')
pdf.body_text(
    'We use data from the MeerKAT Galactic Centre Legacy Survey '
    '(proposal SCI-20210212-SS-01). The observation targets Sgr A* in '
    'L-band (856-1712 MHz) with 32,768 spectral channels at a native '
    'channel width of 26.123 kHz. The total measurement set (MS) is '
    '~17 TB. The data were pre-calibrated by the MeerKAT Science Data '
    'Processor (SDP) pipeline prior to archive delivery.'
)

# === 2. Subband Splitting ===
pdf.section('2. Subband Splitting')
pdf.body_text(
    'The monolithic MS was split into 86 independent subbands of 383 '
    'channels each (~10 MHz per subband) using the CASA split task '
    '(casatasks.split, CASA 6.x). Splitting was performed serially to '
    'avoid file-locking contention on the Lustre parallel filesystem. '
    'Each split operation selected a contiguous range of channels via '
    'the spw parameter (spw="0:start~end") and wrote the DATA column '
    'to a new MS. The typical split time was ~90 minutes per subband. '
    'No frequency regridding was applied; the TOPO-to-LSRK Doppler '
    'correction is <0.1 channel widths for a single observation and is '
    'applied as a metadata correction in post-processing.'
)

# === 3. Imaging ===
pdf.section('3. Imaging')
pdf.body_text(
    'Channel images were produced using the CASA tclean task in cube '
    'mode (specmode="cube"), generating a 383-channel image cube per '
    'subband in a single tclean call. Each cube was then decomposed '
    'into individual single-channel FITS files using imsubimage and '
    'exportfits, yielding one 512x512 pixel image per spectral channel '
    'at 2 arcsec pixel scale (between 1/3 and 1/5 of the synthesized '
    'beam across the band). Briggs weighting with robust=0.5 was used. '
    'The full dataset produces 32,768 channel images spanning '
    '856-1712 MHz at the native 26.123 kHz resolution.'
)
pdf.body_text(
    'Dirty images (niter=0) are the direct Fourier inversion '
    'of the calibrated visibilities with no deconvolution, providing an '
    'unbiased representation of the sky convolved with the synthesized '
    'beam (PSF). The dirty image noise is ~0.09 Jy/beam per 26 kHz '
    'channel in the mid-band (~1.1 GHz).'
)

# Figure: example dirty image
pdf.figure(
    os.path.join(PLOTS_DIR, 'dirty_image_example.png'),
    'Figure 1: Representative dirty channel image at ~1.1 GHz. '
    'Left: linear scale showing the thermal noise and source structure. '
    'Right: log scale revealing faint emission and PSF sidelobes across the field.'
)

# Figure: zoomed Sgr A*
pdf.figure(
    os.path.join(PLOTS_DIR, 'sgra_zoom.png'),
    'Figure 2: Zoomed view of the central 200x200 arcsec (Sgr A* region) at three '
    'frequencies within one subband. The dominant compact source is clearly detected '
    'in each 26 kHz channel.'
)

# === 4. Computational Resources ===
pdf.section('4. Computational Resources')
pdf.body_text(
    'All processing was performed on the Lawrencium HPC cluster at '
    'UC Berkeley, managed by SLURM. Subband splitting was submitted as '
    'a serial array job (max 8 concurrent tasks to match the 8 Lustre '
    'OSTs over which the MS is striped). Imaging was fully parallelized '
    'as a SLURM array job with one subband per compute node (56 cores, '
    '240 GB RAM per node, lr7 partition). Dirty imaging required ~50 '
    'minutes per subband wall-clock time.'
)

# === 5. Data Products ===
pdf.section('5. Data Products')
pdf.body_text(
    'The primary data products are individual FITS images for each of '
    'the 32,768 spectral channels. '
    'Each FITS file contains a single Stokes I image of 512x512 pixels '
    'at 2 arcsec resolution, in units of Jy/beam. Channel frequencies '
    'are encoded in the FITS filename (chan_NNNN_FFFF.FFFMHz.fits) '
    'using the exact channel frequency from the MS metadata.'
)

# Figure: cross-subband comparison
pdf.figure(
    os.path.join(PLOTS_DIR, 'cross_subband.png'),
    'Figure 3: Mid-channel dirty images from six subbands spanning the L-band. '
    'The source morphology evolves with frequency as expected from the '
    'frequency-dependent synthesized beam and sky brightness distribution.',
    width=170
)

# === 6. Frequency Coverage and RFI ===
pdf.section('6. Frequency Coverage and RFI')
pdf.body_text(
    'The MeerKAT L-band receiver covers 856-1712 MHz. Several known '
    'radio frequency interference (RFI) sources affect portions of the '
    'band, including GSM-900 downlink (925-960 MHz), aircraft DME/SSR '
    '(1030-1090 MHz), GPS L5/Galileo E5 (1164-1214 MHz), GPS L2 '
    '(1217-1237 MHz), GPS L1/Galileo E1 (1559-1591 MHz), Glonass L1 '
    '(1598-1610 MHz), and Iridium (1610-1618 MHz). No RFI flagging '
    'has been applied to the current images; the data were delivered '
    'with the MeerKAT SDP calibration but without post-delivery RFI '
    'excision. The relatively RFI-quiet region around 1.3 GHz is of '
    'primary interest for the axion search.'
)

# Figure: RMS and peak vs frequency
pdf.figure(
    os.path.join(PLOTS_DIR, 'rms_peak_vs_freq.png'),
    'Figure 4: Peak flux density (top) and RMS noise (bottom) of the mid-channel '
    'dirty image for each completed subband across the L-band. Shaded regions indicate '
    'known MeerKAT RFI bands. Elevated noise is visible near the GSM-900 and '
    'aircraft DME bands.'
)

# Figure: waterfall
pdf.figure(
    os.path.join(PLOTS_DIR, 'waterfall.png'),
    'Figure 5: Waterfall plot showing flux density along a horizontal slice through '
    'the image center as a function of frequency for one subband. Left: full dynamic range. '
    'Right: stretched scale (5th-95th percentile). Horizontal stripes would indicate '
    'RFI-affected channels; the smooth vertical structure traces the central source.'
)

# === 7. Systematics ===
pdf.section('7. Known Systematics')
pdf.body_text(
    'The first channel of each subband exhibits elevated peak flux '
    '(approximately 3x the median), consistent with a bandpass edge '
    'effect introduced by the channel splitting. These edge channels '
    'are excluded from the science analysis. The FITS WCS frequency '
    'axis (CRVAL3) records the tclean cube reference frequency rather '
    'than the per-channel frequency; the true channel frequency is '
    'obtained from the filename or computed via '
    'CRVAL3 + (1 - CRPIX3) * CDELT3.'
)

# === 8. Validation ===
pdf.section('8. Image Validation')
pdf.body_text(
    'Several diagnostic checks were performed on the dirty images to '
    'verify data quality and pipeline correctness.'
)

# Figure: noise histogram
pdf.figure(
    os.path.join(PLOTS_DIR, 'noise_histogram.png'),
    'Figure 6: Pixel value distribution for a representative mid-band channel '
    'compared to a Gaussian fit. Left: linear scale showing the noise core is '
    'well-described by thermal noise. Right: log scale revealing the positive '
    'tail from real astronomical emission. The Gaussian noise core confirms '
    'proper calibration.'
)

# Figure: beam vs frequency
pdf.figure(
    os.path.join(PLOTS_DIR, 'beam_vs_freq.png'),
    'Figure 7: Synthesized beam major and minor axis as a function of frequency. '
    'The beam size decreases at higher frequencies as expected (beam ~ wavelength/baseline). '
    'Shaded bands indicate known RFI regions.',
    width=170
)

# === 9. Science Application ===
pdf.section('9. Science Application')
pdf.body_text(
    'The native 26.123 kHz channel resolution is preserved throughout '
    'the imaging pipeline, as required for the axion-photon conversion '
    'search. The expected axion signal line width may be comparable to '
    'the channel width, making frequency downbinning unacceptable. '
    'Each channel image will be analyzed independently for narrow '
    'spectral features consistent with axion dark matter conversion in '
    'neutron star magnetospheres across the 856-1712 MHz band '
    '(corresponding to axion masses of ~3.5-7.1 ueV).'
)

# === Summary Table ===
pdf.section('Summary Table')
pdf.ln(2)
rows = [
    ('Telescope',           'MeerKAT (64 dishes)'),
    ('Band',                'L-band, 856-1712 MHz'),
    ('Channels',            '32,768 @ 26.123 kHz'),
    ('Target',              'Sgr A* / Galactic Center'),
    ('Calibration',         'MeerKAT SDP pipeline (pre-calibrated)'),
    ('Image size',          '512 x 512 pixels, 2 arcsec/pixel'),
    ('Weighting',           'Briggs, robust = 0.5'),
    ('Dirty image noise',   '~0.09 Jy/beam per channel (mid-band)'),
    ('Data products',       '32,768 FITS channel images'),
    ('Total MS size',       '~17 TB'),
    ('Processing cluster',  'Lawrencium (UC Berkeley HPC)'),
    ('Time per subband',    '~50 min (dirty imaging)'),
]
for label, value in rows:
    pdf.table_row(label, value)

pdf.ln(5)

# === Software ===
pdf.section('10. Software')
pdf.body_text(
    'CASA 6.x (modular installation via pip: casatools, casatasks). '
    'Image analysis and visualization: astropy, matplotlib, numpy. '
    'Cluster job management: SLURM. All scripts are available in the '
    'project repository under scripts/.'
)

pdf.output(PDF_PATH)
print(f'PDF written to: {PDF_PATH}')
print(f'Pages: {pdf.page_no()}')
