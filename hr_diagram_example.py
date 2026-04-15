import os
import alignement
import calibrations
import photometry
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from glob import glob
from astropy.io import fits
import logging

extension = '.fit'
data_dir = 'omega/B/'
dark_path = 'calibration/darks/dark_60.fit'
flat_path = 'calibration/flat_b.fit'
dark = fits.open(dark_path)[0]
flat = fits.open(flat_path)[0]

data_files = glob(os.path.join(data_dir, '*'+extension))
data_files.sort()
sciences = []
for data_file in data_files:
    print(f"processing image: {data_file}")
    data_fits = fits.open(data_file)[0]
    sci = calibrations.calibrate(data_fits, flat, dark, subframe=[[1000,-1000],[300,-300]])
    sciences.append(sci)

# alignement
print("start alignement")
original_stack = sciences # np.array(sciences) <- this works
aligned_stack = np.array(alignement.phase_correlation_alignment(original_stack))
print("alignement done")

# stack
stacked_median = alignement.stack(aligned_stack, method='median')
photometry.plot_image(stacked_median,'Calibrated & Aligned & Stacked Images')

# background
b,rms = photometry.estimate_background(stacked_median)
photometry.plot_image(b,'Background Map')
photometry.plot_image(rms,'Background Error')
photometry.plot_image(stacked_median-b, 'Background Substracted Image')

# sources
s,c = photometry.extract_sources(stacked_median, b, rms, threshold_sigma=10, npixels=40)

