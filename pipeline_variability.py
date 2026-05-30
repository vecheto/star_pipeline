import os
import alignement
import calibrations
import photometry
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from glob import glob
from astropy.time import Time
from astropy.io import fits
import logging
import matplotlib.patches as patches
from matplotlib.patches import Circle

######################################################################################
# INITIALIZATION
######################################################################################

target_name = 'HATS-26'
extension = '.fits'
data_dir = '../HATS-26_transit'
results_folder = data_dir + "/results/"
os.makedirs(results_folder, exist_ok=True)
dark100_path = '../RH200Calibs/master_dark_100s.fits'
dark120_path = '../RH200Calibs/master_dark_120s.fits'
flat_path = '../RH200Calibs/master_flat_L.fits'
dark100 = fits.open(dark100_path)[0]
dark120 = fits.open(dark120_path)[0]
flat = fits.open(flat_path)[0]

# stars
stars = [
    [2242, 1792],
    #[2337, 1606],
    [2207, 2027],
    [2350, 1760],
    [2056, 1846],
    [2223, 1870],
    [2183, 1840],
   # [2251, 1744],
    [2209, 1757]
         ]

stamp_size = 10
subframe_size = 400

######################################################################################
# CALIBRATION & ALIGNEMENT
######################################################################################

# load science images
data_files = glob(os.path.join(data_dir, '*'+extension))
data_files = calibrations.sort_by_date(data_files)
sciences = []
times = []
for data_file in data_files:
    print(f"processing image: {data_file}")
    data_fits = fits.open(data_file)[0]
    
    # time indexation
    ref_time = Time(data_fits.header['DATE-OBS'], format='isot').jd
    times.append(ref_time)

    # calibration
    if int(data_fits.header['exptime']) == 100:
        sci = calibrations.calibrate(data_fits, flat, dark100)
        
    if int(data_fits.header['exptime']) == 120:
        sci = calibrations.calibrate(data_fits, flat, dark120)
    
    sci = calibrations.extract_subframe(sci, stars[0], subframe_size)

    sciences.append(sci)

# alignement
print("start alignement")
original_stack = np.array(sciences) # <- this works
aligned_stack = np.array(alignement.phase_correlation_alignment(original_stack))
print("alignement done")

######################################################################################
# PHOTOMETRY
######################################################################################

print("start photometry")

stamps = []

transformed_stars = calibrations.coords_to_subframe(stars, stars[0], subframe_size)

for star in transformed_stars:
    star_stamps = photometry.cut_stamp(aligned_stack, star, stamp_size)
    stamps.append(star_stamps)
    img = np.max(star_stamps, axis=0)
    mask = img > 6e4
    fig, ax = plt.subplots()
    ax.imshow(img, origin="lower", interpolation="nearest")
    for y, x in np.argwhere(mask):
        rect = patches.Rectangle(
            (x - 0.5, y - 0.5),
            1, 1,
            linewidth=0,
            facecolor="red",
            alpha=0.5
        )
        ax.add_patch(rect)
    n_sat = np.sum(img > 6e4)
    plt.title(f"Max stamp, star: {star}, number of saturated pix: {n_sat}")
    plt.show()



stamps = np.array(stamps)
fluxes = np.zeros((len(stars), len(aligned_stack)))
fluxes_err =  np.zeros((len(stars), len(aligned_stack)))

######################################################################################
# OPTIMIZE APERTURE AND SKIES
######################################################################################

params = photometry.optimize_parameters(np.median(stamps[0], axis=0), step=0.5, fwhm_min=0.5, fwhm_max=5)
photometry.plot_optimization(params, results_folder)
radio_apertura = params['best_params']['aperture']
radiosky_in = params['best_params']['sky_in']
radiosky_out = params['best_params']['sky_out']


# photometry for every star and stamp
for star in enumerate(stars):
    star_number = star[0]
    for stamp in enumerate(stamps[star_number]):
        stamp_number = stamp[0]
        print(f"processing star : {star_number}, stamp: {stamp_number}")

        centroide = photometry.centroid(stamp[1])
        flux, flux_err = photometry.aperture_phot(stamp[1], centroide, radio_apertura, radiosky_in, radiosky_out)
        fluxes[star_number, stamp_number] = flux
        fluxes_err[star_number, stamp_number] = flux_err
        
        if star_number==1:
            hdu = fits.PrimaryHDU(stamp[1])
            os.makedirs(results_folder + "/stamps/", exist_ok=True)
            hdu.writeto(results_folder + "/stamps/stamp"+str(stamp_number)+".fits", overwrite=True)
        
    # --- PLOT DE CONTROL: stamp mediana + aperturas ---
    median_stamp = np.median(stamps[star_number], axis=0)
    centroide_med = photometry.centroid(median_stamp)

    fig, ax = plt.subplots()
    ax.imshow(median_stamp, origin="lower", interpolation="nearest")
    ax.set_title(f"Star {star_number} – Aperture check")
    plt.colorbar(ax.images[0], ax=ax, label="Flux")

    # Centroide
    ax.plot(
        centroide_med[0], centroide_med[1],
        marker="+", color="cyan", markersize=10, label="Centroid"
    )
    # Apertura
    ap = Circle(
        centroide_med, radio_apertura,
        edgecolor="lime", facecolor="none", linewidth=1.5, label="Aperture"
    )
    ax.add_patch(ap)
    # Sky inner
    sky_in = Circle(
        centroide_med, radiosky_in,
        edgecolor="orange", facecolor="none", linestyle="--", label="Sky in"
    )
    ax.add_patch(sky_in)
    # Sky outer
    sky_out = Circle(
        centroide_med, radiosky_out,
        edgecolor="red", facecolor="none", linestyle="--", label="Sky out"
    )
    ax.add_patch(sky_out)
    ax.legend(loc="upper right")
    plt.savefig(results_folder + f'{star_number}_aperture_check.png')
    plt.show()
            
print("photometry done.")

######################################################################################
# SAVING
######################################################################################

fluxes = fluxes.transpose()   
fluxes_err = fluxes_err.transpose()
fluxes_df = pd.DataFrame()
fluxes_df['time'] = np.array(times)
for star in enumerate(stars):
    fluxes_df['flux_star_'+str(star[0])] = fluxes[:,star[0]]
    fluxes_df['err_flux_star_'+str(star[0])] = fluxes_err[:,star[0]]
fluxes_df.to_csv(results_folder + '/fluxes.csv', index=False)