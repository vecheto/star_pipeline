"""
Differential-photometry pipeline: calibrate -> align -> plate solve -> photometry.

Usage, from the folder with the science FITS:
    python3 /path/to/star_pipeline/pipeline_variability.py pipeline_params.py

photometry_mode (in the params file):
    'aperture_common'    one aperture optimized on the target, applied to every star
    'aperture_per_star'  aperture optimized independently for each star
    'isophotal'          isophotal flux, no fixed aperture

Relative paths in the params file are resolved from the params file's folder.
Outputs go to <params folder>/<results_dir>/<photometry_mode>/.
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
from astropy.visualization import ZScaleInterval
from scipy.spatial import cKDTree
import matplotlib.patches as patches
from matplotlib.patches import Circle

######################################################################################
# INITIALIZATION
######################################################################################

MODES = ('aperture_common', 'aperture_per_star', 'isophotal')
DEFAULTS = {
    'photometry_mode': 'aperture_common',
    'extension': '.fits',
    'results_dir': 'results',
    'target_xy': None,
    'search_radius': 2.0,
    'solve_timeout': 120,
    'plate_solver': 'auto',
    'solve_downsample': None,
    'n_refs': 25,
    'stamp_size': 10,
    'saturation': 60000,
    'edge_margin': 80,
    'fwhm_min': 0.5,
    'fwhm_max': 7,
    'fwhm_step': 0.5,
    'align_rotation': False,
    'show_plots': False,
    'exofop': False,
}

if len(sys.argv) != 2:
    sys.exit(__doc__)

params_path = os.path.abspath(sys.argv[1])
work_dir = os.path.dirname(params_path)
p = {**DEFAULTS, **{k: v for k, v in runpy.run_path(params_path).items() if not k.startswith('__')}}

mode = p['photometry_mode']
if mode not in MODES:
    sys.exit(f"photometry_mode must be one of {MODES}, got '{mode}'")


def resolve(path):
    return os.path.join(work_dir, os.path.expanduser(path))


if not p['show_plots']:
    matplotlib.use('Agg')

stamp_size = p['stamp_size']
results_folder = os.path.join(resolve(p['results_dir']), mode, '')
os.makedirs(results_folder + 'stamps', exist_ok=True)

if p.get('flat') is None:
    print("no flat given: using flat = 1")
    flat = 1.0
else:
    flat = fits.getdata(resolve(p['flat']))

if p.get('darks') is None:
    print("no darks given: using dark = 0")
    darks = None
else:
    darks = {int(t): fits.getdata(resolve(path)) for t, path in p['darks'].items()}


def finish_figure(path):
    plt.savefig(path)
    if p['show_plots']:
        plt.show()
    plt.close('all')


######################################################################################
# CALIBRATION
######################################################################################

data_files = glob(os.path.join(work_dir, '*' + p['extension']))
if not data_files:
    sys.exit(f"no '*{p['extension']}' files in {work_dir}")
data_files = calibrations.sort_by_date(data_files)
sciences = []
times = []
exptimes = []
for data_file in data_files:
    print(f"processing image: {data_file}")
    data_fits = fits.open(data_file)[0]
    sci = data_fits.data.astype(float)
    if not sciences:
        ref_header = data_fits.header

    times.append(Time(data_fits.header['DATE-OBS'], format='isot').jd)
    exptime = float(data_fits.header['EXPTIME'])
    exptimes.append(exptime)

    if darks is None:
        dark = 0.0
    else:
        exptime = int(exptime)
        if exptime not in darks:
            sys.exit(f"no master dark for EXPTIME={exptime} s ({data_file})")
        dark = darks[exptime]
    sciences.append(calibrations.calibrate(sci, flat, dark, isfit=False))

######################################################################################
# ALIGNMENT
######################################################################################

print("start alignement")
if p['align_rotation']:
    aligned, angles, shifts, corrs = alignement.rotation_alignment(sciences)
    aligned_stack = np.array(aligned)
    # pixels not covered by every frame (rotated-out corners): fill so downstream code works,
    # and reject sources there (below)
    covered = np.isfinite(aligned_stack).all(axis=0)
    for frame in aligned_stack:
        frame[~np.isfinite(frame)] = np.nanmedian(frame)
    fig, ax = plt.subplots(2, 1, sharex=True)
    ax[0].plot(times, angles, '.')
    ax[0].set_ylabel('rotation [deg]')
    ax[1].plot(times, corrs, '.')
    ax[1].set_ylabel('correlation')
    ax[1].set_xlabel('JD')
    plt.savefig(results_folder + 'alignment.png')
    plt.close(fig)
else:
    aligned_stack = np.array(alignement.phase_correlation_alignment(np.array(sciences)))
    covered = np.ones(aligned_stack.shape[1:], dtype=bool)
del sciences
print("alignement done")

######################################################################################
# ASTROMETRY — locate target on the reference frame
######################################################################################

# detect on the median stack: hot pixels and cosmics drift relative to the stars after
# alignment and are rejected by the median. Saturation is checked on the max stack.
median_stack = np.median(aligned_stack, axis=0)
max_stack = aligned_stack.max(axis=0)
catalog = sources.extract_sources(median_stack, saturation_threshold=None, edge_margin=p['edge_margin'])
catalog['peak_raw'] = [max_stack[int(round(y)) - 2:int(round(y)) + 3, int(round(x)) - 2:int(round(x)) + 3].max()
                       for x, y in catalog[['x', 'y']].values]
if not covered.all():
    from scipy.ndimage import binary_erosion
    inner = binary_erosion(covered, iterations=2 * stamp_size)
    keep = inner[catalog['y'].round().astype(int), catalog['x'].round().astype(int)]
    print(f"{(~keep).sum()} sources dropped: not covered by every frame after rotation")
    catalog = catalog[keep].reset_index(drop=True)


target_xy, target_method = plate_solve.locate_target(
    aligned_stack[0], ref_header, catalog, p['target_name'], target_xy=p['target_xy'],
    search_radius=p['search_radius'], timeout=p['solve_timeout'], solver=p['plate_solver'],
    cache_dir=work_dir, use_exofop=p['exofop'], downsample=p['solve_downsample'])
print(f"target at pixel ({target_xy[0]:.1f}, {target_xy[1]:.1f}) — {target_method}")

######################################################################################
# REFERENCE STARS AND STAMPS
######################################################################################

# references: unsaturated, isolated (no other detection inside the stamp), away from the target
xy = catalog[['x', 'y']].values
nn_dist = cKDTree(xy).query(xy, k=2)[0][:, 1]
dist = np.hypot(catalog['x'] - target_xy[0], catalog['y'] - target_xy[1])
usable = catalog[(catalog['peak_raw'] < p['saturation']) & (nn_dist > 1.5 * stamp_size) & (dist > 2 * stamp_size)]
ref_xys = sources.brightest_xy(usable, n=p['n_refs'])
stars = np.array([target_xy] + ref_xys).astype(int)
print(f"{len(ref_xys)} reference stars")

vmin, vmax = ZScaleInterval().get_limits(median_stack)
fig, ax = plt.subplots(figsize=(8, 8))
ax.imshow(median_stack, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
for star_number, (x, y) in enumerate(stars):
    color = "red" if star_number == 0 else "cyan"
    ax.add_patch(Circle((x, y), 2 * stamp_size, edgecolor=color, facecolor="none", linewidth=1.2))
    ax.annotate("T" if star_number == 0 else str(star_number), (x + 2 * stamp_size, y + 2 * stamp_size),
                color=color, fontsize=10)
ax.set_title(f"Finder — target (T): {target_method}")
finish_figure(results_folder + 'finder.png')

stamps = []
for star_number, star in enumerate(stars):
    star_stamps = photometry.cut_stamp(aligned_stack, star, stamp_size)
    stamps.append(star_stamps)
    fits.writeto(results_folder + f'stamps/star_{star_number}.fits', star_stamps, overwrite=True)

    img = np.max(star_stamps, axis=0)
    mask = img > p['saturation']
    if mask.any():
        print(f"WARNING: star {star_number} has {mask.sum()} saturated px in at least one frame")
    fig, ax = plt.subplots()
    ax.imshow(img, origin="lower", interpolation="nearest")
    for y, x in np.argwhere(mask):
        ax.add_patch(patches.Rectangle((x - 0.5, y - 0.5), 1, 1,
                                       linewidth=0, facecolor="red", alpha=0.5))
    ax.set_title(f"Max stamp, star {star_number} {star}, saturated px: {mask.sum()}")
    finish_figure(results_folder + f'{star_number}_max_stamp.png')

stamps = np.array(stamps)
n_stars, n_frames = len(stars), len(aligned_stack)
fluxes = np.zeros((n_stars, n_frames))
errors = np.zeros((n_stars, n_frames))

######################################################################################
# PHOTOMETRY
######################################################################################


def optimize(star_number):
    median_stamp = np.median(stamps[star_number], axis=0)
    res = photometry.optimize_parameters(median_stamp, step=p['fwhm_step'],
                                         fwhm_min=p['fwhm_min'], fwhm_max=p['fwhm_max'])
    photometry.plot_optimization(res, results_folder + f'star_{star_number}_')
    plt.close('all')
    best = res['best_params']
    print(f"  star {star_number}: aperture={best['aperture']:.2f}  sky_in={best['sky_in']:.2f}  "
          f"sky_out={best['sky_out']:.2f}  SNR={res['best_snr']:.1f}")
    return best


if mode == 'isophotal':
    print("start isophotal photometry")
    for star_number in range(n_stars):
        for frame, stamp in enumerate(stamps[star_number]):
            fluxes[star_number, frame], errors[star_number, frame] = photometry.isophotal_phot(stamp)
    err_label = 'pixel_std_star_'

else:
    print(f"start aperture photometry ({mode})")
    if mode == 'aperture_common':
        apertures = [optimize(0)] * n_stars
    else:
        apertures = [optimize(i) for i in range(n_stars)]

    for star_number in range(n_stars):
        ap = apertures[star_number]
        median_stamp = np.median(stamps[star_number], axis=0)
        median_centre = photometry.centroid(median_stamp)
        n_fallback = 0
        for frame, stamp in enumerate(stamps[star_number]):
            try:
                centre = photometry.centroid(stamp)
            except Exception:
                centre = (np.nan, np.nan)
            # a non-converged fit can land far away; frames are aligned, so the median centroid is valid
            if np.hypot(centre[0] - median_centre[0], centre[1] - median_centre[1]) > ap['aperture'] \
                    or not np.isfinite(centre).all():
                centre = median_centre
                n_fallback += 1
            fluxes[star_number, frame], errors[star_number, frame] = photometry.aperture_phot(
                stamp, centre, ap['aperture'], ap['sky_in'], ap['sky_out'])
        if n_fallback:
            print(f"  star {star_number}: median centroid used in {n_fallback}/{n_frames} frames")

        centre = median_centre
        fig, ax = plt.subplots()
        im = ax.imshow(median_stamp, origin="lower", interpolation="nearest")
        plt.colorbar(im, ax=ax, label="Flux")
        ax.plot(*centre, marker="+", color="cyan", markersize=10, label="Centroid")
        ax.add_patch(Circle(centre, ap['aperture'], edgecolor="lime", facecolor="none",
                            linewidth=1.5, label="Aperture"))
        ax.add_patch(Circle(centre, ap['sky_in'], edgecolor="orange", facecolor="none",
                            linestyle="--", label="Sky in"))
        ax.add_patch(Circle(centre, ap['sky_out'], edgecolor="red", facecolor="none",
                            linestyle="--", label="Sky out"))
        ax.set_title(f"Star {star_number} – Aperture check")
        ax.legend(loc="upper right")
        finish_figure(results_folder + f'{star_number}_aperture_check.png')
    err_label = 'err_flux_star_'

print("photometry done.")

######################################################################################
# SAVING
######################################################################################

df = pd.DataFrame({'time': np.array(times), 'exptime': np.array(exptimes),
                   'n_sat_target': (stamps[0] > p['saturation']).sum(axis=(1, 2))})
for i in range(n_stars):
    df[f'flux_star_{i}'] = fluxes[i]
    df[f'{err_label}{i}'] = errors[i]
df.to_csv(results_folder + 'fluxes.csv', index=False)

plt.figure(figsize=(10, 4))
plt.plot(times, exptimes, 'o-', drawstyle='steps-post')
plt.xlabel('Time (JD)')
plt.ylabel('EXPTIME (s)')
plt.title(f'Exposure time ({n_frames} frames)')
plt.grid(alpha=0.3)
plt.tight_layout()
finish_figure(results_folder + 'exptime_vs_jd.png')

pd.DataFrame(stars, columns=['x', 'y']).rename_axis('star').to_csv(results_folder + 'stars.csv')
print(f"results saved in {results_folder}")
