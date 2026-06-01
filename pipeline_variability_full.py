import os
import alignement
import calibrations
import photometry
import sources
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
dark100 = fits.open(dark100_path)[0].data
dark120 = fits.open(dark120_path)[0].data
flat = fits.open(flat_path)[0].data

target_star = [2242, 1792]
n_refs = 15

stamp_size = 15
subframe_size = 400

######################################################################################
# CALIBRATION
######################################################################################

# load science images
data_files = glob(os.path.join(data_dir, '*'+extension))
data_files = calibrations.sort_by_date(data_files)
dark100 = calibrations.extract_subframe(dark100, target_star, subframe_size)
dark120 = calibrations.extract_subframe(dark120, target_star, subframe_size)
flat = calibrations.extract_subframe(flat, target_star, subframe_size)
sciences = []
times = []
for data_file in data_files:
    print(f"processing image: {data_file}")
    data_fits = fits.open(data_file)[0]
    sci = calibrations.extract_subframe(data_fits.data, target_star, subframe_size)
    
    # time indexation
    ref_time = Time(data_fits.header['DATE-OBS'], format='isot').jd
    times.append(ref_time)

    # calibration
    if int(data_fits.header['exptime']) == 100:
        sci = calibrations.calibrate(sci, flat, dark100, isfit=False)

    if int(data_fits.header['exptime']) == 120:
        sci = calibrations.calibrate(sci, flat, dark120, isfit=False)

    sciences.append(sci)

######################################################################################
# ALIGNMENT
######################################################################################

print("start alignement")
original_stack = np.array(sciences)
aligned_stack = np.array(alignement.phase_correlation_alignment(original_stack))
print("alignement done")

######################################################################################
# PHOTOMETRY VISUALIZATION
######################################################################################

print("start photometry")

# auto-detect reference stars from first aligned frame
_catalog = sources.extract_sources(
    aligned_stack[0],
    saturation_threshold=60000,
    edge_margin=80,
)
ref_xys = sources.brightest_xy(_catalog, n=n_refs)
target_xy = [float(subframe_size), float(subframe_size)]
transformed_stars = [target_xy] + ref_xys
transformed_stars = np.array(transformed_stars).astype(int)
print(f"target: {target_xy}  |  {len(ref_xys)} reference stars detected")

stamps = []

# SHOW STAMPS
for star in transformed_stars:
    star_stamps = photometry.cut_stamp(aligned_stack, star, stamp_size)
    stamps.append(star_stamps)
    img = np.max(star_stamps, axis=0)
    mask = img > 6e4
    fig, ax = plt.subplots()
    ax.imshow(img, origin="lower", interpolation="nearest")

    # SHOW SATURATED PIXELS
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
fluxes = np.zeros((len(transformed_stars), len(aligned_stack)))
fluxes_err = np.zeros((len(transformed_stars), len(aligned_stack)))

######################################################################################
# OPTIMIZE APERTURE AND SKIES — PER STAR
######################################################################################

star_params = []

for star_number, star in enumerate(transformed_stars):
    label = "target" if star_number == 0 else f"ref_{star_number}"
    print(f"\nOptimizing aperture for star {star_number} ({label})...")

    median_stamp = np.median(stamps[star_number], axis=0)
    params = photometry.optimize_parameters(median_stamp, step=0.5, fwhm_min=0.5, fwhm_max=7)
    star_params.append(params)

    photometry.plot_optimization(params, results_folder + f'star_{star_number}_')
    print(f"  aperture={params['best_params']['aperture']:.2f}  "
          f"sky_in={params['best_params']['sky_in']:.2f}  "
          f"sky_out={params['best_params']['sky_out']:.2f}  "
          f"SNR={params['best_snr']:.1f}")

######################################################################################
# PHOTOMETRY
######################################################################################

for star_number, star in enumerate(transformed_stars):
    ap = star_params[star_number]['best_params']
    radio_apertura = ap['aperture']
    radiosky_in    = ap['sky_in']
    radiosky_out   = ap['sky_out']

    for stamp_number, stamp in enumerate(stamps[star_number]):
        print(f"processing star: {star_number}, stamp: {stamp_number}")
        centroide = photometry.centroid(stamp)
        flux, flux_err = photometry.aperture_phot(stamp, centroide, radio_apertura, radiosky_in, radiosky_out)
        fluxes[star_number, stamp_number] = flux
        fluxes_err[star_number, stamp_number] = flux_err

        if star_number == 1:
            hdu = fits.PrimaryHDU(stamp)
            os.makedirs(results_folder + "/stamps/", exist_ok=True)
            hdu.writeto(results_folder + "/stamps/stamp"+str(stamp_number)+".fits", overwrite=True)

    # --- PLOT DE CONTROL: stamp mediana + aperturas ---
    median_stamp = np.median(stamps[star_number], axis=0)
    centroide_med = photometry.centroid(median_stamp)

    fig, ax = plt.subplots()
    ax.imshow(median_stamp, origin="lower", interpolation="nearest")
    ax.set_title(f"Star {star_number} – Aperture check")
    plt.colorbar(ax.images[0], ax=ax, label="Flux")

    ax.plot(centroide_med[0], centroide_med[1],
            marker="+", color="cyan", markersize=10, label="Centroid")
    ax.add_patch(Circle(centroide_med, radio_apertura,
                        edgecolor="lime", facecolor="none", linewidth=1.5, label="Aperture"))
    ax.add_patch(Circle(centroide_med, radiosky_in,
                        edgecolor="orange", facecolor="none", linestyle="--", label="Sky in"))
    ax.add_patch(Circle(centroide_med, radiosky_out,
                        edgecolor="red", facecolor="none", linestyle="--", label="Sky out"))
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
for i in range(len(transformed_stars)): 
    fluxes_df['flux_star_'+str(i)] = fluxes[:,i] 
    fluxes_df['err_flux_star_'+str(i)] = fluxes_err[:,i]
fluxes_df.to_csv(results_folder + '/fluxes.csv', index=False)
