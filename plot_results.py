"""
Relative light curves (ensemble method, ../brown_dwarfs/metodo.md) from pipeline_variability.py outputs.

Usage, from the data folder:
    python3 /path/to/star_pipeline/plot_results.py pipeline_params.py

Processes every <results_dir>/<photometry_mode>/fluxes.csv found, writes the plots, relative_phot.csv
and star_selection.txt into that same folder, and a mode comparison into <results_dir>/modes_comparison.png.
"""
import os
import sys
import runpy
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import relative

MODES = ('aperture_common', 'aperture_per_star', 'isophotal')
POINTS_PER_BIN = 4

if len(sys.argv) != 2:
    sys.exit(__doc__)

params_path = os.path.abspath(sys.argv[1])
p = runpy.run_path(params_path)
results_root = os.path.join(os.path.dirname(params_path), p.get('results_dir', 'results'))
show_plots = p.get('show_plots', False)
clip_sigma = p.get('clip_sigma', 3)
title = p['target_name'] if isinstance(p.get('target_name'), str) else 'target'
if not show_plots:
    matplotlib.use('Agg')


def finish_figure(path):
    plt.gca().invert_yaxis()
    plt.savefig(path)
    if show_plots:
        plt.show()
    plt.close('all')


def bin_curve(time, y):
    n_bins = max(1, len(time) // POINTS_PER_BIN)
    b = pd.DataFrame({'time': time, 'y': y})
    b['bin'] = pd.cut(b['time'], bins=n_bins)
    s = b.groupby('bin', observed=True).agg(time=('time', 'mean'), y=('y', 'mean'), std=('y', 'std'))
    return s.dropna(subset=['y'])


def rms(x):
    x = np.asarray(x)
    return np.sqrt(np.nanmean((x - np.nanmedian(x))**2))


def light_curve_plot(time, dm, err, label, path, figsize=(10, 4)):
    r = rms(dm)
    plt.figure(figsize=figsize)
    plt.errorbar(time, dm, yerr=err, fmt='o', capsize=2, color='darksalmon')
    plt.axhspan(np.nanmedian(dm) - r, np.nanmedian(dm) + r, alpha=0.1, label=f'RMS = {r:.4f} mag')
    plt.xlabel('Time (JD)')
    plt.ylabel(r'$\Delta m$ - median')
    plt.title(label)
    plt.legend()
    finish_figure(path)


def analyse(folder):
    df = pd.read_csv(folder + 'fluxes.csv')
    time = df['time'].values
    n_stars = sum(c.startswith('flux_star_') for c in df.columns)
    flux = np.column_stack([df[f'flux_star_{i}'] for i in range(n_stars)])
    has_err = 'err_flux_star_0' in df
    flux_err = np.column_stack([df[f'err_flux_star_{i}'] for i in range(n_stars)]) if has_err else None
    mag, mag_err = relative.to_mag(flux, flux_err)

    # saturated target frames are dropped from the light curve (references are already unsaturated)
    saturated = df['n_sat_target'].values > 0 if 'n_sat_target' in df else np.zeros(len(df), bool)
    mag_target_all = mag[:, 0].copy()
    mag[saturated, 0] = np.nan
    if saturated.any():
        print(f"  {saturated.sum()} frames with saturated target excluded")

    # references: every star but the target, measured in all frames
    ref_ids = np.arange(1, n_stars)
    complete = relative.complete_stars(mag[:, ref_ids])
    dropped = ref_ids[~complete]
    ref_ids = ref_ids[complete]
    ref_err = mag_err[:, ref_ids] if has_err else None

    mag_ens, mag_ens_err, rejected = relative.ensemble_star(mag[:, ref_ids], ref_err, sigma=clip_sigma)
    t_err = mag_err[:, 0] if has_err else None
    dm, dm_err = relative.relative_phot(mag[:, 0], t_err, mag_ens, mag_ens_err)
    dm_norm = dm - np.nanmedian(dm)
    # same curve keeping saturated frames, same zero point as the clean one
    dm_all, _ = relative.relative_phot(mag_target_all, t_err, mag_ens, mag_ens_err)
    dm_all_norm = dm_all - np.nanmedian(dm)

    out = pd.DataFrame({'time': time, 'mag_target': mag_target_all, 'mag_target_err': t_err,
                        'mag_ens': mag_ens, 'mag_ens_err': mag_ens_err,
                        'mag_relative': dm, 'mag_relative_err': dm_err, 'delta_mag': dm_norm,
                        'delta_mag_with_saturated': dm_all_norm,
                        'N_refs_used': (~rejected).sum(axis=1), 'target_saturated': saturated})
    out.to_csv(folder + 'relative_phot.csv', index=False)

    # star_selection.txt
    n_rej = rejected.sum(axis=0)
    rej_frames = [(f, ref_ids[j]) for f, j in zip(*np.nonzero(rejected))]
    lines = [
        "target: star 0",
        f"{saturated.sum()} frames excluded (saturated target): {np.nonzero(saturated)[0].tolist()}",
        f"sigma clip on the ensemble: per frame, |delta_i - mean| > {clip_sigma} std, iterative",
        f"{len(ref_ids) + len(dropped)} candidate references; {len(dropped)} dropped (not in every frame): "
        f"{dropped.tolist()}",
        f"{len(ref_ids)} references: {ref_ids.tolist()}",
        "" if len(ref_ids) > 10 else
        f"NOTE: with N = {len(ref_ids)} <= 10 references the {clip_sigma}-sigma clip cannot reject anything "
        "(max deviation is sqrt(N-1) std); raise n_refs to make it effective.",
        f"{rejected.sum()} rejections (star, frame) out of {rejected.size}; per star:",
        pd.Series(n_rej, index=ref_ids, name='n_rejected').rename_axis('star').to_string(),
        "",
        "rejections (frame, star): " + (str(rej_frames) if rej_frames else "none"),
    ]
    with open(folder + 'star_selection.txt', 'w') as f:
        f.write("\n".join(lines) + "\n")
    print(f"  {len(ref_ids)} refs, {len(dropped)} dropped, {rejected.sum()} clipped points")

    # plots
    light_curve_plot(time, dm_norm, dm_err, f'{title} – {os.path.basename(folder[:-1])}',
                     folder + 'rel_lightcurve.png')
    b = bin_curve(time, dm_norm)
    light_curve_plot(b['time'], b['y'], b['std'], f'{title} – binned', folder + 'rel_lightcurve_binned.png')

    plt.figure(figsize=(10, 4))
    plt.errorbar(time[~saturated], dm_all_norm[~saturated],
                 yerr=None if dm_err is None else dm_err[~saturated],
                 fmt='o', capsize=2, color='darksalmon', label='unsaturated')
    plt.plot(time[saturated], dm_all_norm[saturated], 'x', color='crimson', ms=7,
             label=f'saturated target ({saturated.sum()})')
    plt.xlabel('Time (JD)')
    plt.ylabel(r'$\Delta m$ - median')
    plt.title(f'{title} – {os.path.basename(folder[:-1])} – with saturated frames')
    plt.legend()
    finish_figure(folder + 'rel_lightcurve_with_saturated.png')

    plt.figure(figsize=(10, 4))
    plt.errorbar(time, mag[:, 0] - np.nanmedian(mag[:, 0]), yerr=t_err, fmt='--o', capsize=2, label='target')
    plt.errorbar(time, mag_ens, yerr=mag_ens_err, fmt='--o', capsize=2, label='ensemble')
    plt.xlabel('Time (JD)')
    plt.ylabel('m - median(m)')
    plt.title(f'{title}: target vs ensemble')
    plt.legend()
    finish_figure(folder + 'rel_target_vs_ensemble.png')

    # each reference against an ensemble that excludes it (variable references stand out)
    plt.figure(figsize=(10, 6))
    for k, star in enumerate(ref_ids):
        others = np.delete(np.arange(len(ref_ids)), k)
        e, ee, _ = relative.ensemble_star(mag[:, ref_ids[others]],
                                          None if ref_err is None else ref_err[:, others], sigma=clip_sigma)
        d, _ = relative.relative_phot(mag[:, star], None, e, ee)
        d -= np.nanmedian(d)
        plt.plot(time, d, '--o', ms=3, label=f'{star} (rms {rms(d):.3f})')
    plt.xlabel('Time (JD)')
    plt.ylabel('m_ref - mag_ens - median')
    plt.title(f'{title}: reference stars')
    plt.legend(ncol=2, fontsize='small')
    finish_figure(folder + 'rel_lightcurves_refstars.png')

    return time, dm_norm


curves = {}
for mode in MODES:
    folder = os.path.join(results_root, mode, '')
    if os.path.exists(folder + 'fluxes.csv'):
        print(f"{mode}:")
        curves[mode] = analyse(folder)

if not curves:
    sys.exit(f"no fluxes.csv found under {results_root}/<mode>/")

fig, axes = plt.subplots(len(curves), 1, sharex=True, sharey=True,
                         figsize=(10, 2.8 * len(curves)), squeeze=False)
for ax, (mode, (time, dm)) in zip(axes[:, 0], curves.items()):
    b = bin_curve(time, dm)
    ax.plot(time, dm, '.', alpha=0.4, label=f'RMS = {rms(dm):.4f}')
    ax.plot(b['time'], b['y'], 'o-', label=f'binned RMS = {rms(b["y"]):.4f}')
    ax.set_title(mode, fontsize=10)
    ax.set_ylabel(r'$\Delta m$')
    ax.legend(fontsize=8, loc='lower left')
axes[-1, 0].set_xlabel('Time (JD)')
fig.suptitle(f'{title} – photometry modes')
fig.tight_layout()
finish_figure(os.path.join(results_root, 'modes_comparison.png'))
print(f"comparison saved in {results_root}/modes_comparison.png")
