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

######################################################################################
# INITIALIZATION
######################################################################################

target_name = 'corot-7'
extension = '.fits'
data_dir = 'COROT7'
#dark_path = 'calibration/darks/dark_100.fit'
#flat_path = 'calibration/flat_l.fit'
#dark = fits.open(dark_path)[0]
#flat = fits.open(flat_path)[0]

# stars
stars = [
         [520,511],
         [476, 547],
         [312, 539],
         [939, 726],
         [278, 432]
         ]

stamp_size = 18

######################################################################################
# CALIBRATION & ALIGNEMENT
######################################################################################

# load science images
data_files = glob(os.path.join(data_dir, '*'+extension))
data_files.sort() # TODO: SORT BY DATEOBS
sciences = []
times = []
for data_file in data_files:
    print(f"processing image: {data_file}")
    data_fits = fits.open(data_file)[0]
    sci = calibrations.calibrate(data_fits, None, None, use_subframe=False)
    sciences.append(sci)

    # time indexation
    ref_time = Time(data_fits.header['DATE-OBS'], format='isot').jd
    times.append(ref_time)

# alignement
print("start alignement")
original_stack = np.array(sciences) #<- this works
aligned_stack = np.array(alignement.phase_correlation_alignment(original_stack))
print("alignement done")

######################################################################################
# PHOTOMETRY
######################################################################################

print("start photometry")

stamps = []
for star in stars:
    star_stamps = photometry.cut_stamp(aligned_stack, star, stamp_size)
    stamps.append(star_stamps)
    plt.imshow(np.sum(star_stamps,axis=0))
    plt.title(f"Sum stamp, star: {star}")
    plt.show()

stamps = np.array(stamps)
fluxes = np.zeros((len(stars), len(aligned_stack)))
fluxes_err =  np.zeros((len(stars), len(aligned_stack)))

######################################################################################
# OPTIMIZE APERTURE AND SKIES
######################################################################################

params = photometry.optimize_parameters(np.median(stamps[0], axis=0), step=0.5, fwhm_min=0.5, fwhm_max=4)
photometry.plot_optimization(params, target_name+'/')
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
            hdu.writeto(target_name + "/stamps/stamp"+str(stamp_number)+".fits", overwrite=True)
            
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
fluxes_df.to_csv(target_name+'/fluxes.csv', index=False)