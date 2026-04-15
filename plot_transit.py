import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# read csv fluxes
fluxes = pd.read_csv('RRPic_2/fluxes.csv')

# compute relative flux and error
cols = [c for c in fluxes.columns if c.startswith("flux_star_") and c != "flux_star_0"]
fluxes_for_median = fluxes[cols]
median_flux = np.median(fluxes_for_median, axis=1)
norm = np.median(fluxes['flux_star_0']/(median_flux))
fluxes['relative_flux'] = fluxes['flux_star_0']/(norm * median_flux)
fluxes['relative_err'] = fluxes['err_flux_star_0'] / (norm * median_flux)


# plot relative_flux
plt.errorbar(fluxes['time'], fluxes['relative_flux'], yerr=fluxes['relative_err'], fmt='o--', capsize=3)
plt.title('RRPic')
plt.xlabel('Time (JD)')
plt.ylabel('Relative Flux')
plt.savefig('relative_flux.png')
plt.show()

# plot every star
stars = [key for key in fluxes.keys() if key != 'time' and 'err' not in key.lower() and 'relative' not in key.lower()]
for star in stars:
    plt.errorbar(fluxes['time'], fluxes[star], yerr=fluxes['err_'+star], fmt='o--', capsize=3)
    plt.title(star)
    plt.xlabel('Time (JD)')
    plt.ylabel('Flux (counts)')
    plt.savefig(star+'.png')
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
plt.figure(figsize=(10, 6))
plt.errorbar( bin_stats[('time', 'mean')], bin_stats[('relative_flux', 'mean')], yerr=bin_stats[('relative_flux', 'std')], fmt='o--', capsize=3)
plt.title('RRPic- Binned')
plt.xlabel('Time (JD)')
plt.ylabel('Relative Flux')
plt.savefig('relative_flux_binned.png')
plt.show()

