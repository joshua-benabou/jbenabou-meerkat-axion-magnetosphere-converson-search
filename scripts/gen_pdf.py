#!/usr/bin/env python3
"""Generate a PDF with paper-ready data processing description."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import textwrap

pdf_path = (
    '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/'
    'meerkat_reduction_project/processing_summary.pdf'
)


class PageWriter:
    def __init__(self, ax):
        self.ax = ax
        self.y = 0.95

    def title(self, s, size=14, dy=0.045):
        self.ax.text(0.0, self.y, s, fontsize=size, fontweight='bold',
                     transform=self.ax.transAxes, verticalalignment='top',
                     fontfamily='serif')
        self.y -= dy

    def heading(self, s, size=11, dy=0.035):
        self.ax.text(0.0, self.y, s, fontsize=size, fontweight='bold',
                     transform=self.ax.transAxes, verticalalignment='top',
                     fontfamily='serif')
        self.y -= dy

    def para(self, s, size=9.5, width=95, dy=0.022):
        """Write wrapped paragraph text."""
        lines = textwrap.wrap(s, width)
        for line in lines:
            self.ax.text(0.0, self.y, line, fontsize=size,
                         transform=self.ax.transAxes, verticalalignment='top',
                         fontfamily='serif')
            self.y -= dy

    def gap(self, dy=0.015):
        self.y -= dy


def make_page():
    fig = plt.figure(figsize=(8.5, 11))
    ax = fig.add_axes([0.08, 0.08, 0.84, 0.86])
    ax.axis('off')
    return fig, PageWriter(ax)


with PdfPages(pdf_path) as pdf:
    # === Page 1 ===
    fig, pw = make_page()

    pw.title('MeerKAT Galactic Center Data Reduction', size=15, dy=0.05)
    pw.para('J. Benabou, UC Berkeley / Safdi Group', size=10)
    pw.gap(dy=0.03)

    pw.heading('1. Observation')
    pw.para(
        'We use data from the MeerKAT Galactic Centre Legacy Survey '
        '(proposal SCI-20210212-SS-01). The observation targets Sgr A* in '
        'L-band (856\u20131712 MHz) with 32,768 spectral channels at a native '
        'channel width of 26.123 kHz. The total measurement set (MS) is '
        '~17 TB. The data were pre-calibrated by the MeerKAT Science Data '
        'Processor (SDP) pipeline prior to archive delivery.'
    )
    pw.gap()

    pw.heading('2. Subband Splitting')
    pw.para(
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
    pw.gap()

    pw.heading('3. Imaging')
    pw.para(
        'Channel images were produced using the CASA tclean task in cube '
        'mode (specmode="cube"), generating a 383-channel image cube per '
        'subband in a single tclean call. Each cube was then decomposed '
        'into individual single-channel FITS files using imsubimage and '
        'exportfits, yielding one 512x512 pixel image per spectral channel '
        'at 2 arcsec pixel scale (between 1/3 and 1/5 of the synthesized '
        'beam across the band). Briggs weighting with robust=0.5 was used. '
        'The full dataset produces 32,768 channel images spanning '
        '856\u20131712 MHz at the native 26.123 kHz resolution.'
    )
    pw.gap()
    pw.para(
        'Two sets of images are produced:'
    )
    pw.gap(dy=0.005)
    pw.para(
        '(a) Dirty images (niter=0): These are the direct Fourier inversion '
        'of the calibrated visibilities with no deconvolution, providing an '
        'unbiased representation of the sky convolved with the synthesized '
        'beam (PSF). The dirty image noise is ~0.09 Jy/beam per 26 kHz '
        'channel in the mid-band (~1.1 GHz).'
    )
    pw.gap(dy=0.005)
    pw.para(
        '(b) Cleaned images (niter=10,000): Deconvolved images using the '
        'Hogbom CLEAN algorithm with auto-multithresh masking '
        '(sidelobethreshold=2.0, noisethreshold=5.0, lownoisethreshold=1.5, '
        'minbeamfrac=0.3). The auto-masking identifies emission regions '
        'without manual intervention, appropriate for the complex Galactic '
        'Center field.'
    )
    pw.gap()

    pw.heading('4. Computational Resources')
    pw.para(
        'All processing was performed on the Lawrencium HPC cluster at '
        'UC Berkeley, managed by SLURM. Subband splitting was submitted as '
        'a serial array job (max 8 concurrent tasks to match the 8 Lustre '
        'OSTs over which the MS is striped). Imaging was fully parallelized '
        'as a SLURM array job with one subband per compute node (56 cores, '
        '240 GB RAM per node, lr7 partition). Dirty imaging required ~50 '
        'minutes per subband; cleaned imaging requires approximately 2\u20134 '
        'hours per subband.'
    )
    pw.gap()

    pw.heading('5. Data Products')
    pw.para(
        'The primary data products are individual FITS images for each of '
        'the 32,768 spectral channels, in both dirty and cleaned versions. '
        'Each FITS file contains a single Stokes I image of 512x512 pixels '
        'at 2 arcsec resolution, in units of Jy/beam. Channel frequencies '
        'are encoded in the FITS filename (chan_NNNN_FFFF.FFFMHz.fits) '
        'using the exact channel frequency from the MS metadata. Dirty '
        'images are stored in images/subband_XXX/ and cleaned images in '
        'images/subband_XXX/cleaned/.'
    )

    pdf.savefig(fig)
    plt.close(fig)

    # === Page 2 ===
    fig, pw = make_page()

    pw.heading('6. Frequency Coverage and RFI')
    pw.para(
        'The MeerKAT L-band receiver covers 856\u20131712 MHz. Several known '
        'radio frequency interference (RFI) sources affect portions of the '
        'band, including GSM-900 downlink (925\u2013960 MHz), aircraft DME/SSR '
        '(1030\u20131090 MHz), GPS L5/Galileo E5 (1164\u20131214 MHz), GPS L2 '
        '(1217\u20131237 MHz), GPS L1/Galileo E1 (1559\u20131591 MHz), Glonass L1 '
        '(1598\u20131610 MHz), and Iridium (1610\u20131618 MHz). No RFI flagging '
        'has been applied to the current images; the data were delivered '
        'with the MeerKAT SDP calibration but without post-delivery RFI '
        'excision. The relatively RFI-quiet region around 1.3 GHz is of '
        'primary interest for the axion search.'
    )
    pw.gap()

    pw.heading('7. Known Systematics')
    pw.para(
        'The first channel of each subband exhibits elevated peak flux '
        '(approximately 3x the median), consistent with a bandpass edge '
        'effect introduced by the channel splitting. These edge channels '
        'are excluded from the science analysis. The FITS WCS frequency '
        'axis (CRVAL3) records the tclean cube reference frequency rather '
        'than the per-channel frequency; the true channel frequency is '
        'obtained from the filename or computed via '
        'CRVAL3 + (1 - CRPIX3) * CDELT3.'
    )
    pw.gap()

    pw.heading('8. Science Application')
    pw.para(
        'The native 26.123 kHz channel resolution is preserved throughout '
        'the imaging pipeline, as required for the axion-photon conversion '
        'search. The expected axion signal line width may be comparable to '
        'the channel width, making frequency downbinning unacceptable. '
        'Each channel image will be analyzed independently for narrow '
        'spectral features consistent with axion dark matter conversion in '
        'neutron star magnetospheres across the 856\u20131712 MHz band '
        '(corresponding to axion masses of ~3.5\u20137.1 \u03bceV).'
    )
    pw.gap()

    pw.heading('9. Software')
    pw.para(
        'CASA 6.x (modular installation via pip: casatools, casatasks). '
        'Image analysis and visualization: astropy, matplotlib, numpy. '
        'Cluster job management: SLURM. All scripts are available in the '
        'project repository under scripts/.'
    )
    pw.gap()

    pw.heading('Summary Table')
    pw.gap(dy=0.005)
    rows = [
        ('Telescope',           'MeerKAT (64 dishes)'),
        ('Band',                'L-band, 856\u20131712 MHz'),
        ('Channels',            '32,768 @ 26.123 kHz'),
        ('Target',              'Sgr A* / Galactic Center'),
        ('Calibration',         'MeerKAT SDP pipeline (pre-calibrated)'),
        ('Image size',          '512 x 512 pixels, 2 arcsec/pixel'),
        ('Weighting',           'Briggs, robust = 0.5'),
        ('Dirty image noise',   '~0.09 Jy/beam per channel (mid-band)'),
        ('Clean niter',         '10,000 (auto-multithresh mask)'),
        ('Data products',       '32,768 FITS per image type (dirty + cleaned)'),
        ('Total MS size',       '~17 TB'),
        ('Processing cluster',  'Lawrencium (UC Berkeley HPC)'),
    ]
    for label, value in rows:
        pw.ax.text(0.02, pw.y, label, fontsize=9, fontweight='bold',
                   transform=pw.ax.transAxes, verticalalignment='top',
                   fontfamily='serif')
        pw.ax.text(0.30, pw.y, value, fontsize=9,
                   transform=pw.ax.transAxes, verticalalignment='top',
                   fontfamily='serif')
        pw.y -= 0.024

    pdf.savefig(fig)
    plt.close(fig)

print(f'PDF written to: {pdf_path}')
