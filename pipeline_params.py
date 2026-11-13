# Template: copy into the science-data folder and edit.
# Relative paths are resolved from this file's folder.

target_name = 'HATS-26'          # resolved with Sesame/SIMBAD (needs internet)
# offline alternative:
# from astropy.coordinates import SkyCoord
# target_name = SkyCoord(ra_deg, dec_deg, unit='deg')
extension = '.fits'

# flat = None -> flat = 1;  darks = None -> dark = 0 (no calibration)
flat = '../RH200Calibs/master_flat_L.fits'
darks = {                         # EXPTIME [s] -> master dark
    100: '../RH200Calibs/master_dark_100s.fits',
    120: '../RH200Calibs/master_dark_120s.fits',
}

results_dir = 'results'

# astrometry: plate solve around the target, then around the header pointing,
# then fall back to the source closest to the image center. Check results/<mode>/finder.png
search_radius = 2.0              # deg
solve_timeout = 120              # s per attempt
target_xy = None                 # (x, y) in pixels to skip astrometry, e.g. (250, 248)

# reference stars
n_refs = 25                      # > 10, or the 3-sigma clip of the ensemble cannot reject anything
clip_sigma = 3                   # ensemble sigma clip (plot_results.py)
saturation = 60000
edge_margin = 50

# photometry
# 'aperture_common'   : one aperture optimized on the target, used for all stars
# 'aperture_per_star' : aperture optimized for each star
# 'isophotal'         : isophotal flux (fwhm_* ignored)
photometry_mode = 'aperture_common'
stamp_size = 10
fwhm_min = 0.5
fwhm_max = 7
fwhm_step = 0.5

show_plots = True               # True opens every figure interactively
