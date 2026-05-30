import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# read csv fluxes
fluxes = pd.read_csv('../HATS-26_transit/results/fluxes.csv')
results_folder = '../HATS-26_transit/results/'
target_name = 'HATS-26'

# compute relative flux and error
cols = [c for c in fluxes.columns if c.startswith("flux_star_") and c != "flux_star_0"]
fluxes_for_median = fluxes[cols]
median_flux = np.median(fluxes_for_median, axis=1)
norm = np.median(fluxes['flux_star_0']/(median_flux))
fluxes['relative_flux'] = fluxes['flux_star_0']/(norm * median_flux)

# Error propagation: r = f0 / (N * m), so sigma_r/r = sqrt((sigma_f0/f0)^2 + (sigma_m/m)^2)
# Uncertainty in median of n comparison stars: sigma_m ~ sqrt(pi/2) * sqrt(sum(sigma_i^2)) / n
err_cols = ['err_' + c for c in cols]
err_matrix = fluxes[err_cols].values
n_comp = len(cols)
sigma_median = np.sqrt(np.pi / 2) * np.sqrt(np.sum(err_matrix**2, axis=1)) / n_comp

fluxes['relative_err'] = fluxes['relative_flux'] * np.sqrt(
    (fluxes['err_flux_star_0'] / fluxes['flux_star_0'])**2 +
    (sigma_median / median_flux)**2
)


# plot relative_flux
rms = np.sqrt(np.mean((fluxes['relative_flux'] - 1)**2))
t_margin = (fluxes['time'].max() - fluxes['time'].min()) * 0.05
t_fill = [fluxes['time'].min() - t_margin, fluxes['time'].max() + t_margin]
plt.errorbar(fluxes['time'], fluxes['relative_flux'], yerr=fluxes['relative_err'], fmt='o--', capsize=3)
plt.fill_between(t_fill, 1 - rms, 1 + rms, alpha=0.1, label=f'RMS = {rms:.4f}')
plt.xlim(t_fill)
plt.title('HATS-26')
plt.xlabel('Time (JD)')
plt.ylabel('Relative Flux')
plt.legend()
plt.savefig(results_folder + 'relative_flux.png')
plt.show()

# plot every star
stars = [key for key in fluxes.keys() if key != 'time' and 'err' not in key.lower() and 'relative' not in key.lower()]
for star in stars:
    plt.errorbar(fluxes['time'], fluxes[star], yerr=fluxes['err_'+star], fmt='o--', capsize=3, label='Star ' + star ) 
    plt.xlabel('Time (JD)')
    plt.ylabel('Flux (counts)')
plt.title('Fluxes')
plt.legend(fontsize=8)
plt.savefig(results_folder + 'stars_flux.png')
plt.show()

stars = [key for key in fluxes.keys() if key != 'time' and 'err' not in key.lower() and 'relative' not in key.lower()]
for star in stars:
    plt.errorbar(fluxes['time'], fluxes[star] - median_flux, yerr=fluxes['err_'+star], fmt='o--', capsize=3, label='Star ' + star ) 
    plt.xlabel('Time (JD)')
    plt.ylabel('Flux (counts)')
plt.title('Fluxes - Median Flux')
plt.legend(fontsize=8)
plt.savefig(results_folder + 'stars_flux_miuns_median.png')
plt.show()

# binning
num_bins = int(len(fluxes)/4) 
binned = fluxes.copy()
binned['time_bin'] = pd.cut(binned['time'], bins=num_bins)

bin_stats = binned.groupby('time_bin').agg({
    'time': 'mean',
    'relative_flux': ['mean', 'std', 'count']
}).reset_index()

# plot binning
rms_binned = np.sqrt(np.mean((bin_stats[('relative_flux', 'mean')] - 1)**2))
t_bin = bin_stats[('time', 'mean')]
t_margin_bin = (t_bin.max() - t_bin.min()) * 0.05
t_fill_bin = [t_bin.min() - t_margin_bin, t_bin.max() + t_margin_bin]
plt.figure(figsize=(10, 6))
plt.errorbar( bin_stats[('time', 'mean')], bin_stats[('relative_flux', 'mean')], yerr=bin_stats[('relative_flux', 'std')], fmt='o--', capsize=3)
plt.fill_between(t_fill_bin, 1 - rms_binned, 1 + rms_binned, alpha=0.1, label=f'RMS = {rms_binned:.4f}')
plt.xlim(t_fill_bin)
plt.title('HATS-26 - Binned')
plt.xlabel('Time (JD)')
plt.ylabel('Relative Flux')
plt.legend()
plt.savefig(results_folder +'relative_flux_binned.png')
plt.show()

