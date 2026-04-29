from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

rows = []

with fits.open("J0559-1404/J0559-1404_NTT_sequence.fits") as hdul:
    for hdu in hdul[1:]:
        if hdu.data is None:
            continue

        #epoch = hdu.header.get("MJD-OBS")
        data = hdu.data

        for obstime, star_id, X, Y, ra, dec, mag, errmag, flux, flux_err in zip(
            data["obstime"],
            data["ID_file"],
            data["X"],
            data["Y"],
            data["RA"],
            data["DEC"],
            data["mag"],
            data["errmag"],
            data['flux'],
            data['flux_err']
        ):
            if np.isnan(star_id):
                continue
            rows.append({
                "epoch": obstime,
                "star_id": int(star_id),
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

# find median mag and median error to characterize star
data_full["mag_median_star"] = data_full.groupby("star_id")["mag"].transform("median")
data_full["errmag_median_star"] = data_full.groupby("star_id")["errmag"].transform("median")

# extract similar stars
target_star_id = 102
target_mag = data_full.loc[data_full["star_id"] == target_star_id, "mag_median_star"].iloc[0]
target_sigma = data_full.loc[data_full["star_id"] == target_star_id, "errmag_median_star"].iloc[0]
df_similar = pd.DataFrame()
star_medians = data_full.groupby("star_id")["mag_median_star"].first()
similar_ids = star_medians[(star_medians - target_mag).abs() <= 5 * target_sigma].index
df_similar = data_full[data_full["star_id"].isin(similar_ids)].copy()

print('similar stars:', df_similar.star_id.unique())

# extract median per epoch
epoch_stats = df_similar.groupby("epoch").agg(
    mag_median=("mag", "median"),
    flux_median=("flux", "median"),
    errmag_median=("errmag", "median"),
    flux_err_median=("flux_err", "median"),
)

plt.errorbar(epoch_stats.index, epoch_stats.mag_median, yerr=epoch_stats['errmag_median'], fmt='--o')
plt.title('mag median')
plt.gca().invert_yaxis()
plt.savefig('figures/magmedian.png')
plt.show()

plt.errorbar(epoch_stats.index, epoch_stats.flux_median, yerr=epoch_stats['flux_err_median'], fmt='--o')
plt.title('flux median')
plt.savefig('figures/fluxmedian.png')
plt.show()

# flujo de las estrellas
for star_id in df_similar['star_id'].unique():
    df_star = df_similar[df_similar["star_id"] == star_id].set_index("epoch")
    dif = pd.DataFrame()
    dif['dif'] = df_star['mag'] - epoch_stats['mag_median']
    dif['errmag'] = df_star['errmag']
    dif = dif.dropna()
    plt.errorbar(dif.index, dif['dif'], yerr=dif['errmag'], fmt='--o')
    plt.xlabel("Epoch")
    plt.ylabel("Magnitude star - median")
    plt.title(f"Light curves")
    plt.gca().invert_yaxis()
    plt.savefig('figures/mag_ref_stars.png')
plt.show()

star_id = 102
df_star = df_similar[df_similar["star_id"] == star_id].set_index("epoch")
dif = pd.DataFrame()
dif['dif'] = df_star['mag'] - epoch_stats['mag_median']
dif['errmag'] = df_star['errmag']
dif = dif.dropna()
plt.errorbar(dif.index, dif['dif'], yerr=dif['errmag'], fmt='--o')
plt.xlabel("Epoch")
plt.ylabel("Magnitude star - median")
plt.title(f"Light curve Target Star")
plt.gca().invert_yaxis()
plt.savefig('figures/mag_target_star.png')
plt.show()

star_id = 102
df_star = df_similar[df_similar["star_id"] == star_id].set_index("epoch")
dif = pd.DataFrame()
dif['dif'] = df_star['flux'] / epoch_stats['flux_median']
dif['flux_err'] = df_star['flux_err']
dif = dif.dropna()
norm = np.median(dif['dif'])
plt.errorbar(dif.index, dif['dif']/norm, yerr=dif['flux_err']/norm, fmt='--o')
plt.xlabel("Epoch")
plt.ylabel("Relative Flux")
plt.title(f"Light curves")
plt.savefig('figures/flux_target_star.png')
plt.show()
