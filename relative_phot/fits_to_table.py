from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import relative

###############################################################
# DATA EXTRACTION
###############################################################

fits_path = "../J0559-1404/J0559-1404_NTT_sequence.fits"
skip_hdus = {17}     
target_star= 102
drop_stars= {55, 75, 167, 195}  

data = relative.read_table(fits_path, skip_hdus, save=True)

# Find median mag and median error to characterize star
data["mag_median_star"] = data.groupby("star_id")["mag"].transform("median")
data["errmag_median_star"] = data.groupby("star_id")["errmag"].transform("median")

# Extract brighter stars
comp = relative.find_brightest(data, 10, drop_stars, target_star)
star =  data[data["star_id"] == target_star].set_index("epoch")

###############################################################
# MEDIAN MAGNITUDE AND FLUX
###############################################################

medians = comp.groupby('epoch').apply(relative.epoch_comp_stats)

plt.errorbar(medians.index, medians.mag_median, yerr=medians['sigma_mag_median'], fmt='--o', label='median mag')
plt.errorbar(star.index, star.mag, yerr=star.errmag, fmt='--o', label='target star')
plt.title('mag median')
plt.gca().invert_yaxis()
plt.legend()
plt.savefig('../figures/mag_median.png')
plt.show()

plt.errorbar(medians.index, medians.flux_median, yerr=medians['sigma_flux_median'], fmt='--o', label='median flux')
plt.errorbar(star.index, star.flux, yerr=star.flux_err, fmt='--o', label='target star')
plt.title('flux median')
plt.legend()
plt.savefig('../figures/flux_median.png')
plt.show()

###############################################################
# RELATIVE MAG OF STARS
###############################################################

for star_id in  np.append(target_star, comp['star_id'].unique()):
    ref_star = data[data["star_id"] == star_id].set_index("epoch")
    rel = pd.DataFrame(index=ref_star.index)

    # relative magnitude
    rel['mag'] = ref_star['mag'] - medians['mag_median']

    # relative error σ_Δm = √(σ_estrella² + σ_mediana²)
    rel['errmag'] = np.sqrt(ref_star['errmag']**2 + medians['sigma_mag_median']**2)

    rel = rel.dropna()
    label = 'target' if star_id == target_star else str(star_id)
    plt.errorbar(rel.index, rel['mag'], yerr=rel['errmag'], fmt='--o', label=label)
plt.xlabel("Epoch")
plt.ylabel("Magnitude star - median")
plt.title("Light curves")
plt.gca().invert_yaxis()
plt.legend()
plt.savefig('../figures/magnitudes_refstars.png')
plt.show()

###############################################################
# RELATIVE TARGET STAR
###############################################################

# RELATIVE MAG
# select target
comparison_stars = [102, 76,  90,  96,  97, 114, 167, 181, 195]
for id in comparison_stars:
    istarget = True if id == 102 else False
    relative.relative_phot(data, medians, id, istarget=istarget)
