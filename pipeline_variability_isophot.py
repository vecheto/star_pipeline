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
import matplotlib.patches as patches

######################################################################################
# INITIALIZATION
######################################################################################

target_name  = 'HATS-26'
extension    = '.fits'
data_dir     = '../HATS-26_transit'
results_folder = data_dir + '/results_isophot/'
os.makedirs(results_folder, exist_ok=True)

dark100_path = '../RH200Calibs/master_dark_100s.fits'
dark120_path = '../RH200Calibs/master_dark_120s.fits'
flat_path    = '../RH200Calibs/master_flat_L.fits'
dark100 = fits.open(dark100_path)[0]
dark120 = fits.open(dark120_path)[0]
flat    = fits.open(flat_path)[0]

target_star    = [2242, 1792]
n_refs         = 15
stamp_size     = 10
subframe_size  = 400
saturation     = 60000

######################################################################################
# CALIBRATION
######################################################################################

data_files = glob(os.path.join(data_dir, '*' + extension))
data_files = calibrations.sort_by_date(data_files)

sciences = []
times    = []
for data_file in data_files:
    print(f"processing image: {data_file}")
    data_fits = fits.open(data_file)[0]

    ref_time = Time(data_fits.header['DATE-OBS'], format='isot').jd
    times.append(ref_time)

    if int(data_fits.header['exptime']) == 100:
        sci = calibrations.calibrate(data_fits, flat, dark100)
    if int(data_fits.header['exptime']) == 120:
        sci = calibrations.calibrate(data_fits, flat, dark120)

    sci = calibrations.extract_subframe(sci, target_star, subframe_size)
    sciences.append(sci)

######################################################################################
# ALIGNMENT
######################################################################################

print("start alignement")
original_stack = np.array(sciences)
aligned_stack  = np.array(alignement.phase_correlation_alignment(original_stack))
print("alignement done")

######################################################################################
# REFERENCE STARS
######################################################################################

print('searching reference stars')
_catalog = sources.extract_sources(
    aligned_stack[0],
    saturation_threshold=saturation,
    edge_margin=70,
)
ref_xys    = sources.brightest_xy(_catalog, n=n_refs)
target_xy  = [float(subframe_size), float(subframe_size)]
all_stars  = np.array([target_xy] + ref_xys).astype(int)
print(f"target: {target_xy}  |  {len(ref_xys)} reference stars detected")

######################################################################################
# STAMPS
######################################################################################

print("cutting stamps")
stamps = []
for star in all_stars:
    star_stamps = photometry.cut_stamp(aligned_stack, star, stamp_size)
    stamps.append(star_stamps)

    img  = np.max(star_stamps, axis=0)
    mask = img > saturation
    fig, ax = plt.subplots()
    ax.imshow(img, origin='lower', interpolation='nearest')
    for y, x in np.argwhere(mask):
        rect = patches.Rectangle(
            (x - 0.5, y - 0.5), 1, 1,
            linewidth=0, facecolor='red', alpha=0.5,
        )
        ax.add_patch(rect)
    n_sat = np.sum(mask)
    ax.set_title(f"Max stamp – star {star}, saturated px: {n_sat}")
    plt.show()

stamps = np.array(stamps)

######################################################################################
# ISOPHOTAL PHOTOMETRY
######################################################################################

print("start isophotal photometry")
fluxes     = np.zeros((len(all_stars), len(aligned_stack)))
fluxes_std = np.zeros((len(all_stars), len(aligned_stack)))

for star_idx, star in enumerate(all_stars):
    for frame_idx, stamp in enumerate(stamps[star_idx]):
        print(f"star {star_idx}  frame {frame_idx}")
        flux ,pixel_std = photometry.isophotal_phot(stamp)
        fluxes[star_idx, frame_idx]     = flux
        fluxes_std[star_idx, frame_idx] = pixel_std

print("photometry done.")

######################################################################################
# SAVING
######################################################################################

fluxes     = fluxes.transpose()
fluxes_std = fluxes_std.transpose()

df = pd.DataFrame()
df['time'] = np.array(times)
for i in range(len(all_stars)):
    df[f'flux_star_{i}']     = fluxes[:, i]
    df[f'pixel_std_star_{i}'] = fluxes_std[:, i]

df.to_csv(results_folder + 'isofluxes.csv', index=False)
print(f"saved → {results_folder}isofluxes.csv")
