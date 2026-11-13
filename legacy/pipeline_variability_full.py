"""
Aperture photometry pipeline (per-star aperture optimization).

Usage, from the folder with the science FITS:
    python3 /path/to/star_pipeline/pipeline_variability_full.py pipeline_params.py

Relative paths in the params file are resolved from the params file's folder.
Outputs go to <params folder>/results/.
"""
import os
import sys
import runpy
import alignement
import plate_solve
import calibrations
import photometry
import sources
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from glob import glob
from astropy.time import Time
from astropy.io import fits
import matplotlib.patches as patches
from matplotlib.patches import Circle

######################################################################################
# INITIALIZATION
######################################################################################

DEFAULTS = {
    'extension': '.fits',
    'results_dir': 'results',
    'search_radius': 2.0,
    'n_refs': 10,
    'stamp_size': 10,
    'saturation': 60000,
    'edge_margin': 80,
    'fwhm_min': 0.5,
    'fwhm_max': 7,
    'fwhm_step': 0.5,
    'show_plots': False,
}

if len(sys.argv) != 2:
    sys.exit(__doc__)

params_path = os.path.abspath(sys.argv[1])
work_dir = os.path.dirname(params_path)
p = {**DEFAULTS, **{k: v for k, v in runpy.run_path(params_path).items() if not k.startswith('__')}}


def resolve(path):
    return os.path.join(work_dir, os.path.expanduser(path))


if not p['show_plots']:
    matplotlib.use('Agg')

target_name = p['target_name']
stamp_size = p['stamp_size']
results_folder = os.path.join(resolve(p['results_dir']), '')
os.makedirs(results_folder, exist_ok=True)

flat = fits.getdata(resolve(p['flat']))
darks = {int(t): fits.getdata(resolve(path)) for t, path in p['darks'].items()}


def finish_figure(path):
    plt.savefig(path)
    if p['show_plots']:
        plt.show()
    plt.close()


######################################################################################
# CALIBRATION
######################################################################################

data_files = glob(os.path.join(work_dir, '*' + p['extension']))
if not data_files:
    sys.exit(f"no '*{p['extension']}' files in {work_dir}")
data_files = calibrations.sort_by_date(data_files)
sciences = []
times = []
for data_file in data_files:
    print(f"processing image: {data_file}")
    data_fits = fits.open(data_file)[0]
    sci = data_fits.data.astype(float)
    if not sciences:
        ref_header = data_fits.header

    times.append(Time(data_fits.header['DATE-OBS'], format='isot').jd)

    exptime = int(data_fits.header['EXPTIME'])
    if exptime not in darks:
        sys.exit(f"no master dark for EXPTIME={exptime} s ({data_file})")
    sci = calibrations.calibrate(sci, flat, darks[exptime], isfit=False)

    sciences.append(sci)

######################################################################################
# ALIGNMENT
######################################################################################

print("start alignement")
original_stack = np.array(sciences)
aligned_stack = np.array(alignement.phase_correlation_alignment(original_stack))
print("alignement done")

######################################################################################
# ASTROMETRY — locate target on the reference frame
######################################################################################

print("plate solving reference frame")
pointing = plate_solve.header_pointing(ref_header)
wcs = plate_solve.solve_wcs(aligned_stack[0], pointing, radius=p['search_radius'])
target_xy = list(plate_solve.target_pixel(wcs, target_name))
print(f"target at pixel {target_xy}")

######################################################################################
# PHOTOMETRY VISUALIZATION
######################################################################################

print("start photometry")

# auto-detect reference stars from first aligned frame
_catalog = sources.extract_sources(
    aligned_stack[0],
    saturation_threshold=p['saturation'],
    edge_margin=p['edge_margin'],
)
_dist = np.hypot(_catalog['x'] - target_xy[0], _catalog['y'] - target_xy[1])
ref_xys = sources.brightest_xy(_catalog[_dist > 2 * stamp_size], n=p['n_refs'])
transformed_stars = [target_xy] + ref_xys
transformed_stars = np.array(transformed_stars).astype(int)
print(f"target: {target_xy}  |  {len(ref_xys)} reference stars detected")

stamps = []

# SHOW STAMPS
for star_number, star in enumerate(transformed_stars):
    star_stamps = photometry.cut_stamp(aligned_stack, star, stamp_size)
    stamps.append(star_stamps)
    img = np.max(star_stamps, axis=0)
    mask = img > p['saturation']
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
    plt.title(f"Max stamp, star: {star}, number of saturated pix: {mask.sum()}")
    finish_figure(results_folder + f'{star_number}_max_stamp.png')

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
    params = photometry.optimize_parameters(median_stamp, step=p['fwhm_step'],
                                            fwhm_min=p['fwhm_min'], fwhm_max=p['fwhm_max'])
    star_params.append(params)

    photometry.plot_optimization(params, results_folder + f'star_{star_number}_')
    plt.close('all')
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
    finish_figure(results_folder + f'{star_number}_aperture_check.png')

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
fluxes_df.to_csv(results_folder + 'fluxes.csv', index=False)
print(f"results saved in {results_folder}")
