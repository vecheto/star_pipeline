# Template: copy into the science-data folder and edit.
# Relative paths are resolved from this file's folder.

target_name = 'VW Hyi'          # resolved with Sesame/SIMBAD (needs internet)
# offline alternative:
# from astropy.coordinates import SkyCoord
# target_name = SkyCoord(ra_deg, dec_deg, unit='deg')
extension = '.fits'

# flat = None -> flat = 1;  darks = None -> dark = 0 (no calibration)
flat = None
darks = None

results_dir = 'results'

# astrometry
search_radius = 2.0              # deg around header pointing

# reference stars
n_refs = 25
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

show_plots = False               # True opens every figure interactively
