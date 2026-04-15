#!/usr/bin/env python3
"""Generate all Sam Witte NS population plots."""

import os
os.environ["PATH"] = '/global/software/sl-7.x86_64/modules/tools/texlive/2016/bin/x86_64-linux/:' + os.environ["PATH"]

import sys
sys.path.insert(0, '/global/scratch/projects/pc_heptheory/jbenabou')
import plotting_defaults

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import warnings

BASE = '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/reference/Real_Analysis/Real_Analysis'
OUT_DIR = BASE + '/Output_Files'
PLOT_DIR = '/global/scratch/projects/pc_heptheory/jbenabou/NS_megaproject/MeerKAT_data/meerkat_reduction_project/plots'

POP_YOUNG = 'PopYoung_WB_TauO_1.00e+07_B_12.79_sB_0.51_P_1.07_sP_14.56_'
POP_OLD = 'PopOld_WB_TauO_1.00e+12_B_12.82_sB_0.56_P_0.19_sP_22.19_'

h_eV = 4.1357e-15  # eV*s

os.makedirs(PLOT_DIR, exist_ok=True)

# === 1. Population Properties ===
pop_young_0 = np.loadtxt(OUT_DIR + '/' + POP_YOUNG + '/Population_Details_Pop_0_.txt')
xE = np.array([0.0, 8.3, 0.0])
dists = np.sqrt(np.sum((pop_young_0[:,10:13] - xE)**2, axis=1))

print(f"Young Pop 0: {len(pop_young_0)} NSs")
print(f"B_final: {pop_young_0[:,4].min():.2e} -- {pop_young_0[:,4].max():.2e} G")
print(f"P_final: {pop_young_0[:,5].min():.4f} -- {pop_young_0[:,5].max():.4f} s")
print(f"Ages: {pop_young_0[:,7].min():.2f} -- {pop_young_0[:,7].max():.2f} Myr")

fig, axes = plt.subplots(2, 3, figsize=(20, 12))

axes[0,0].hist(np.log10(pop_young_0[:,4]), bins=50, color='cornflowerblue', alpha=0.8)
axes[0,0].set_xlabel(r'$\log_{10}(B_{\rm final}$ / G)')
axes[0,0].set_ylabel(r'Count')
axes[0,0].set_title(r'Magnetic Field (evolved)')

axes[0,1].hist(np.log10(pop_young_0[:,5]), bins=50, color='forestgreen', alpha=0.8)
axes[0,1].set_xlabel(r'$\log_{10}(P_{\rm final}$ / s)')
axes[0,1].set_ylabel(r'Count')
axes[0,1].set_title(r'Spin Period (evolved)')

axes[0,2].hist(pop_young_0[:,7], bins=50, color='maroon', alpha=0.8)
axes[0,2].set_xlabel(r'Age [Myr]')
axes[0,2].set_ylabel(r'Count')
axes[0,2].set_title(r'Age Distribution')

axes[1,0].scatter(pop_young_0[:,10], pop_young_0[:,11], s=1, alpha=0.3, c='cornflowerblue')
axes[1,0].plot(0, 8.3, 'r*', markersize=15, label=r'Earth')
axes[1,0].plot(0, 0, 'k+', markersize=15, label=r'GC')
axes[1,0].set_xlabel(r'$x$ [kpc]')
axes[1,0].set_ylabel(r'$y$ [kpc]')
axes[1,0].set_title(r'Spatial Distribution')
axes[1,0].legend(fontsize=14)
axes[1,0].set_aspect('equal')

axes[1,1].hist(dists, bins=50, color='goldenrod', alpha=0.8)
axes[1,1].set_xlabel(r'Distance from Earth [kpc]')
axes[1,1].set_ylabel(r'Count')
axes[1,1].set_title(r'Distance Distribution')

r_gc = np.sqrt(pop_young_0[:,10]**2 + pop_young_0[:,11]**2 + pop_young_0[:,12]**2)
axes[1,2].scatter(r_gc, pop_young_0[:,8], s=1, alpha=0.3, c='firebrick')
axes[1,2].set_xlabel(r'$r_{\rm GC}$ [kpc]')
axes[1,2].set_ylabel(r'$\rho_{\rm DM}$ [GeV/cm$^3$]')
axes[1,2].set_title(r'DM Density vs GC Distance')
axes[1,2].set_yscale('log')
axes[1,2].set_xscale('log')

plt.suptitle(r'Young NS Population (Pop 0, $\tau_{\rm ohm} = 10$ Myr, %d NSs)' % len(pop_young_0), fontsize=18)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(PLOT_DIR + '/sam_young_pop_properties.png', dpi=150)
plt.close()
print("Saved: sam_young_pop_properties.png")

# === 2. Flux Spectrum at Ma=4.13e-6 ===
flux_file = OUT_DIR + '/' + POP_YOUNG + '/Ma_4.130e-06/Pop_0/Combined_Flux.dat'
flux_data = np.loadtxt(flux_file)
freqs_mhz = flux_data[:,2] / h_eV / 1e6
MassA = 4.13e-6

print(f"\nCombined flux: {len(flux_data)} photon samples, {len(np.unique(flux_data[:,0]))} unique NSs")
print(f"Freq range: {freqs_mhz.min():.3f} -- {freqs_mhz.max():.3f} MHz")

bins = np.linspace(freqs_mhz.min(), freqs_mhz.max(), 200)
bw = (bins[1] - bins[0]) * 1e6
spectrum, edges = np.histogram(freqs_mhz, bins=bins, weights=flux_data[:,1])
spectrum /= bw
centers = 0.5 * (edges[:-1] + edges[1:])

fig, axes = plt.subplots(1, 2, figsize=(20, 7))

axes[0].plot(centers, spectrum * 1e6, color='cornflowerblue', lw=0.8)
axes[0].set_xlabel(r'Frequency [MHz]')
axes[0].set_ylabel(r'Flux Density [$\mu$Jy] at $g_{a\gamma\gamma} = 10^{-12}$ GeV$^{-1}$')
axes[0].set_title(r'Young Pop 0 Spectrum ($m_a = 4.13\;\mu$eV)')
axes[0].axvline(MassA / h_eV / 1e6, color='maroon', ls='--', alpha=0.5, label=r'$\nu = m_a / h$')
axes[0].legend()

core_mask = np.abs(centers - np.median(freqs_mhz)) < 0.5
if np.any(core_mask):
    axes[1].plot(centers[core_mask], spectrum[core_mask] * 1e6, color='cornflowerblue', lw=1)
    axes[1].set_xlabel(r'Frequency [MHz]')
    axes[1].set_ylabel(r'Flux Density [$\mu$Jy]')
    axes[1].set_title(r'Zoom: Central $\pm 0.5$ MHz')
    axes[1].axvline(MassA / h_eV / 1e6, color='maroon', ls='--', alpha=0.5, label=r'$\nu = m_a / h$')
    axes[1].legend()

plt.tight_layout()
plt.savefig(PLOT_DIR + '/sam_flux_spectrum_ma4p13.png', dpi=150)
plt.close()
print("Saved: sam_flux_spectrum_ma4p13.png")

# === 3. Per-NS Flux Contributions ===
ns_indices = np.unique(flux_data[:,0]).astype(int)
ns_total_flux = np.zeros(len(ns_indices))
for i, ns_idx in enumerate(ns_indices):
    mask = flux_data[:,0] == ns_idx
    ns_total_flux[i] = np.sum(flux_data[mask, 1])

sort_idx = np.argsort(ns_total_flux)[::-1]
ns_indices_sorted = ns_indices[sort_idx]
ns_flux_sorted = ns_total_flux[sort_idx]
cumflux = np.cumsum(ns_flux_sorted) / np.sum(ns_flux_sorted)

fig, axes = plt.subplots(1, 3, figsize=(20, 6))

axes[0].hist(np.log10(ns_total_flux[ns_total_flux > 0]), bins=50, color='cornflowerblue', alpha=0.8)
axes[0].set_xlabel(r'$\log_{10}$(Total flux weight / Jy-Hz)')
axes[0].set_ylabel(r'Count')
axes[0].set_title(r'Per-NS Flux Distribution')

axes[1].plot(np.arange(1, len(cumflux)+1), cumflux * 100, color='forestgreen', lw=2)
axes[1].set_xlabel(r'Number of brightest NSs')
axes[1].set_ylabel(r'Cumulative flux [\%]')
axes[1].set_title(r'Cumulative Flux (brightest first)')
axes[1].axhline(90, color='maroon', ls='--', alpha=0.5)
axes[1].axhline(50, color='goldenrod', ls='--', alpha=0.5)
n_50 = np.searchsorted(cumflux, 0.5) + 1
n_90 = np.searchsorted(cumflux, 0.9) + 1
axes[1].text(n_50 * 1.5, 52, r'%d NSs = 50\%%' % n_50, fontsize=14)
axes[1].text(n_90 * 1.1, 92, r'%d NSs = 90\%%' % n_90, fontsize=14)

top_ns_locs = pop_young_0[ns_indices_sorted[:50], 10:13]
sc = axes[2].scatter(top_ns_locs[:,0], top_ns_locs[:,1],
                c=np.log10(ns_flux_sorted[:50]), cmap='hot', s=50, edgecolors='k', lw=0.5)
axes[2].plot(0, 8.3, 'c*', markersize=15, label=r'Earth')
axes[2].plot(0, 0, 'g+', markersize=15, label=r'GC')
axes[2].set_xlabel(r'$x$ [kpc]')
axes[2].set_ylabel(r'$y$ [kpc]')
axes[2].set_title(r'Top 50 brightest NSs')
axes[2].legend(fontsize=14)
axes[2].set_aspect('equal')
plt.colorbar(sc, ax=axes[2], label=r'$\log_{10}$(flux / Jy-Hz)')

plt.tight_layout()
plt.savefig(PLOT_DIR + '/sam_per_ns_flux.png', dpi=150)
plt.close()
print("Saved: sam_per_ns_flux.png")
print(f"  {n_50} NSs contribute 50% of flux, {n_90} contribute 90%")

# === 4. Three masses comparison ===
masses = {'3.540e-06': 3.54e-6, '4.130e-06': 4.13e-6, '7.080e-06': 7.08e-6}
colors_m = {'3.540e-06': 'cornflowerblue', '4.130e-06': 'forestgreen', '7.080e-06': 'maroon'}

fig, axes = plt.subplots(1, 2, figsize=(20, 7))

for mass_label, mass_val in masses.items():
    ff = OUT_DIR + '/' + POP_YOUNG + '/Ma_' + mass_label + '/Pop_0/Combined_Flux.dat'
    if not os.path.exists(ff):
        print(f"  {mass_label}: no Combined_Flux.dat")
        continue
    fd = np.loadtxt(ff)
    freqs = fd[:,2] / h_eV / 1e6
    freq_center = mass_val / h_eV / 1e6
    b = np.linspace(freqs.min(), freqs.max(), 200)
    bw2 = (b[1] - b[0]) * 1e6
    sp, ed = np.histogram(freqs, bins=b, weights=fd[:,1])
    sp /= bw2
    ct = 0.5 * (ed[:-1] + ed[1:])
    n_ns = len(np.unique(fd[:,0]))
    total_f = np.sum(fd[:,1]) / (mass_val * 1e-5 / 6.58e-16)
    lab = r'$m_a = %.2f\;\mu$eV ($\nu_0 = %.0f$ MHz)' % (mass_val*1e6, freq_center)
    axes[0].plot(ct, sp * 1e6, color=colors_m[mass_label], lw=1, label=lab)
    delta_f = (ct - freq_center) / freq_center * 1e3
    axes[1].plot(delta_f, sp * 1e6, color=colors_m[mass_label], lw=1, label=lab)
    print(f"  {mass_label}: {n_ns} NSs, total flux ~ {total_f:.3e} Jy")

axes[0].set_xlabel(r'Frequency [MHz]')
axes[0].set_ylabel(r'Flux Density [$\mu$Jy] at $g = 10^{-12}$ GeV$^{-1}$')
axes[0].set_title(r'Young Pop 0: Three Axion Masses')
axes[0].legend(fontsize=13)
axes[1].set_xlabel(r'$(\nu - \nu_0) / \nu_0 \times 10^3$')
axes[1].set_ylabel(r'Flux Density [$\mu$Jy]')
axes[1].set_title(r'Normalized Frequency Offset')
axes[1].legend(fontsize=13)

plt.tight_layout()
plt.savefig(PLOT_DIR + '/sam_three_masses_comparison.png', dpi=150)
plt.close()
print("Saved: sam_three_masses_comparison.png")

# === 5. Young vs Old ===
fig, axes = plt.subplots(1, 2, figsize=(20, 7))
for pop_tag, pop_dir, color, label in [
    ('Young', POP_YOUNG, 'cornflowerblue', r'Young ($\tau_{\rm ohm} = 10$ Myr)'),
    ('Old', POP_OLD, 'maroon', r'Old ($\tau_{\rm ohm} = 1$ Tyr)'),
]:
    ff = OUT_DIR + '/' + pop_dir + '/Ma_4.130e-06/Pop_0/Combined_Flux.dat'
    if not os.path.exists(ff):
        print(f"  {pop_tag}: missing")
        continue
    fd = np.loadtxt(ff)
    freqs = fd[:,2] / h_eV / 1e6
    b = np.linspace(freqs.min(), freqs.max(), 200)
    bw2 = (b[1] - b[0]) * 1e6
    sp, ed = np.histogram(freqs, bins=b, weights=fd[:,1])
    sp /= bw2
    ct = 0.5 * (ed[:-1] + ed[1:])
    n_ns = len(np.unique(fd[:,0]))
    total_f = np.sum(fd[:,1]) / (4.13e-6 * 1e-5 / 6.58e-16)
    axes[0].plot(ct, sp * 1e6, color=color, lw=1,
                 label=label + r' (%d NSs)' % n_ns)
    pop_file = OUT_DIR + '/' + pop_dir + '/Population_Details_Pop_0_.txt'
    pop = np.loadtxt(pop_file)
    axes[1].hist(np.log10(pop[:,4]), bins=50, color=color, alpha=0.5,
                 label=label + ' (%d NSs)' % len(pop))
    print(f"{pop_tag}: {n_ns} converting / {len(pop)} total, flux ~ {total_f:.3e} Jy")

axes[0].set_xlabel(r'Frequency [MHz]')
axes[0].set_ylabel(r'Flux Density [$\mu$Jy] at $g = 10^{-12}$ GeV$^{-1}$')
axes[0].set_title(r'Young vs Old at $m_a = 4.13\;\mu$eV')
axes[0].legend(fontsize=13)
axes[1].set_xlabel(r'$\log_{10}(B_{\rm final}$ / G)')
axes[1].set_ylabel(r'Count')
axes[1].set_title(r'Magnetic Field Distributions')
axes[1].legend(fontsize=13)

plt.tight_layout()
plt.savefig(PLOT_DIR + '/sam_young_vs_old.png', dpi=150)
plt.close()
print("Saved: sam_young_vs_old.png")

print("\nAll plots saved to:", PLOT_DIR)
