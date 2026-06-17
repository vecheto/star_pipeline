
from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def read_table(fits_path, skip_hdus, save=True):
    rows = []
    with fits.open(fits_path) as hdul:
        for i, hdu in enumerate(hdul[1:], start=1):
            if i in skip_hdus:
                continue
            if hdu.data is None:
                continue

            #epoch = hdu.header.get("MJD-OBS")
            data = hdu.data
            for target_id, (obstime, X, Y, ra, dec, mag, errmag, flux, flux_err) in enumerate(zip(
                data["obstime"],
                data["X"],
                data["Y"],
                data["RA"],
                data["DEC"],
                data["mag"],
                data["errmag"],
                data['flux'],
                data['flux_err']
            )):
                rows.append({
                    "epoch": obstime,
                    "star_id": target_id,
                    "X": X,
                    "Y": Y,
                    "ra": (360/(2*np.pi)) * float(ra),
                    "dec": (360/(2*np.pi)) * float(dec),
                    "mag": mag,
                    "errmag": errmag,
                    "flux": flux,
                    "flux_err": flux_err
                })

    # data
    data_full = pd.DataFrame(rows)
    if save:
        data_full.to_csv('data_table.csv', index=False)
    return data_full

def find_brightest(data, n, drop_stars_id, target_star_id):
    top_star_medians = data.groupby("star_id")["mag_median_star"].first()
    brightest_ids = top_star_medians.drop(index=target_star_id, errors="ignore").nsmallest(n).index
    comparison_stars = data[ (data["star_id"].isin(brightest_ids )) & (~data["star_id"].isin(drop_stars_id)) &
                            (~(data["star_id"]==target_star_id))].copy()
    print('comparison stars:', comparison_stars.star_id.unique())
    return comparison_stars

def epoch_comp_stats(group):
    n = len(group)
    # σ_mediana = √(π/2) · √(Σ σᵢ²) / n
    return pd.Series({
        'mag_median':        group['mag'].median(),
        'flux_median':       group['flux'].median(),
        'sigma_mag_median':  np.sqrt(np.pi / 2) * np.sqrt((group['errmag']**2).sum()) / n,
        'sigma_flux_median': np.sqrt(np.pi / 2) * np.sqrt((group['flux_err']**2).sum()) / n,
    })

def relative_phot(data, medians, star_id, istarget=False, plot= True):
    star = data[data["star_id"] == star_id].set_index("epoch")
    rel = pd.DataFrame(index=star.index)    
    
    # RELATIVE MAGNITUDE
    rel['mag'] = star['mag'] - medians['mag_median']

    # MAG ERROR
    # σ_Δm = √(σ_estrella² + σ_mediana²)
    rel['errmag'] = np.sqrt(star['errmag']**2 + medians['sigma_mag_median']**2)
    rel = rel.dropna()

    # RELATIVE FLUX
    rel['flux'] = star['flux'] / medians['flux_median']
    rel = rel.dropna()
    norm = np.median(rel['flux'])
    rel['flux'] /= norm

    # FLUX ERROR
    # σ_r / r = √( (σ_f/f)² + (σ_m/m)² )
    rel['flux_err'] = rel['flux'] * np.sqrt(
        (star['flux_err'] / star['flux'])**2 +
        (medians['sigma_flux_median'] / medians['flux_median'])**2
    )

    # PLOTS
    if plot:
        plt.figure(figsize=(10,4))
        plt.errorbar(rel.index, rel['mag'], yerr=rel['errmag'], fmt='o')
        plt.xlabel("Epoch")
        plt.ylabel("Magnitude star - median")
        plt.title("Light curve Star " + str(star_id))
        plt.gca().invert_yaxis()
        path = '../figures/mag_target.png' if istarget else '../figures/reference_stars/mag_'+str(star_id)+'_star.png'
        plt.savefig(path)
        plt.show()

        plt.figure(figsize=(10,4))
        plt.errorbar(rel.index, rel['flux'], yerr=rel['flux_err'], fmt='o', capsize=3)
        plt.xlabel("Epoch")
        plt.ylabel("Relative Flux")
        plt.title("Light curve Star " + str(star_id))
        path = '../figures/flux_target.png' if istarget else '../figures/reference_stars/flux_'+str(star_id)+'_star.png'
        plt.savefig(path)
        plt.show()

    return rel