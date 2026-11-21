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
import exofop
import plate_solve
import relative

MODES = ('aperture_common', 'aperture_per_star', 'isophotal')
POINTS_PER_BIN = 3

if len(sys.argv) != 2:
    sys.exit(__doc__)

params_path = os.path.abspath(sys.argv[1])
p = runpy.run_path(params_path)
results_root = os.path.join(os.path.dirname(params_path), p.get('results_dir', 'results'))
show_plots = p.get('show_plots', False)
clip_sigma = p.get('clip_sigma', 3)
title = p.get('target_label') or (p['target_name'] if isinstance(p.get('target_name'), str) else 'target')
if not show_plots:
    matplotlib.use('Agg')

# Only with `exofop = True` (or a manual `ephemeris`): target coordinates -> BJD_TDB time axis,
# and predicted transit marked. ExoFOP answers are cached next to the params file;
# `ephemeris = (T0 [BJD_TDB], P [d], duration [h])` overrides the ExoFOP one.
use_exofop = p.get('exofop', False)
toi = None
target_coord = None
if use_exofop or p.get('ephemeris'):
    try:
        if use_exofop:
            if not exofop.parse_toi(p.get('target_name')):
                sys.exit(f"exofop = True needs a 'TOI-<n>[.<pp>]' target_name, got {p.get('target_name')!r}")
            toi = exofop.toi_info(p['target_name'], cache_dir=os.path.dirname(params_path))
        target_coord = plate_solve.to_skycoord(p.get('target_name'), cache_dir=os.path.dirname(params_path),
                                               use_exofop=use_exofop)
    except Exception as e:
        print(f"WARNING: could not resolve the target ({e}); time axis stays in JD_UTC")
if p.get('ephemeris'):
    T0, P, dur = p['ephemeris']
    toi = {**(toi or {}), 'epoch': T0, 'period': P, 'duration': dur, 'epoch_err': 0, 'period_err': 0}
TIME_LABEL = 'BJD$_{TDB}$ (mid-exposure)' if target_coord is not None else 'Time (JD$_{UTC}$, exposure start)'
transits = []


def mark_transits(ax=None):
    """Predicted ingress/egress (dashed) and mid-transit (dotted, shaded ±1σ) of `transits`."""
    ax = ax or plt.gca()
    for k, tr in enumerate(transits):
        for t in (tr['t1'], tr['t4']):
            ax.axvline(t, color='k', ls='--', lw=1.2, label='predicted ingress/egress' if k == 0 and t == tr['t1'] else None)
        ax.axvline(tr['tc'], color='k', ls=':', lw=1)
        if tr['tc_err'] > 0:
            ax.axvspan(tr['tc'] - tr['tc_err'], tr['tc'] + tr['tc_err'], color='k', alpha=0.08, lw=0)


def finish_figure(path):
    plt.gca().invert_yaxis()
    plt.savefig(path)
    if show_plots:
        plt.show()
    plt.close('all')


def bin_curve(time, y, gap_factor=5):
    """Bins of POINTS_PER_BIN consecutive valid points; a bin never spans a gap
    (> gap_factor x median cadence). A segment's leftover points go into its last bin."""
    b = pd.DataFrame({'time': time, 'y': y}).dropna().sort_values('time').reset_index(drop=True)
    dt = np.diff(b['time'].values)
    segment = np.concatenate([[0], np.cumsum(dt > gap_factor * np.median(dt))]) if len(dt) else np.zeros(len(b), int)
    b['bin'] = 0
    next_bin = 0
    for seg in np.unique(segment):
        idx = np.nonzero(segment == seg)[0]
        k = np.arange(len(idx)) // POINTS_PER_BIN
        n_full = max(1, len(idx) // POINTS_PER_BIN)
        b.loc[idx, 'bin'] = next_bin + np.minimum(k, n_full - 1)
        next_bin += n_full
    s = b.groupby('bin').agg(time=('time', 'mean'), y=('y', 'mean'), std=('y', 'std'),
                                                 n=('y', 'count'))
    s['sem'] = s['std'] / np.sqrt(s['n'])
    return s.dropna(subset=['y'])


def rms(x):
    x = np.asarray(x)
    return np.sqrt(np.nanmean((x - np.nanmedian(x))**2))


def light_curve_plot(time, dm, err, label, path, figsize=(10, 4), binned=None):
    """binned: output of bin_curve, drawn on top of the (faded) unbinned points."""
    plt.figure(figsize=figsize)
    if binned is None:
        plt.errorbar(time, dm, yerr=err, fmt='o', capsize=2, color='darksalmon')
    else:
        plt.errorbar(time, dm, yerr=err, fmt='o', capsize=0, color='darksalmon', alpha=0.4,
                     label='data')
        plt.errorbar(binned['time'], binned['y'], yerr=binned['sem'], fmt='o', ms=7, capsize=2,
                     color='tab:blue', zorder=3,
                     label=f'binned ({POINTS_PER_BIN} pts)')
    mark_transits()
    plt.xlabel(TIME_LABEL)
    plt.ylabel(r'$\Delta m$ - median')
    plt.title(label)
    if binned is not None or transits:
        plt.legend()
    finish_figure(path)


def analyse(folder):
    global transits
    df = pd.read_csv(folder + 'fluxes.csv')
    jd_utc = df['time'].values                   # DATE-OBS = exposure start
    if target_coord is not None:
        mid = jd_utc + df['exptime'].values / 2 / 86400 if 'exptime' in df else jd_utc
        time = exofop.jd_utc_to_bjd_tdb(mid, target_coord)
    else:
        time = jd_utc
    transits = exofop.transits_in(toi, time.min(), time.max()) if toi else []
    for tr in transits:
        print(f"  predicted transit (epoch {tr['epoch']}): Tc = {tr['tc']:.5f} +/- {tr['tc_err']:.5f} BJD_TDB, "
              f"ingress {tr['t1']:.5f}, egress {tr['t4']:.5f}")
    if toi and not transits:
        print("  no predicted transit inside this observation")
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

    out = pd.DataFrame({'time': jd_utc, 'bjd_tdb': time if target_coord is not None else np.nan,
                        'mag_target': mag_target_all, 'mag_target_err': t_err,
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
    b = bin_curve(time, dm_norm)
    light_curve_plot(time, dm_norm, dm_err, f'{title} – {os.path.basename(folder[:-1])}',
                     folder + 'rel_lightcurve.png', binned=b)
    light_curve_plot(b['time'], b['y'], b['std'], f'{title} – binned', folder + 'rel_lightcurve_binned.png')

    plt.figure(figsize=(10, 4))
    plt.errorbar(time[~saturated], dm_all_norm[~saturated],
                 yerr=None if dm_err is None else dm_err[~saturated],
                 fmt='o', capsize=2, color='darksalmon', label='unsaturated')
    plt.plot(time[saturated], dm_all_norm[saturated], 'x', color='crimson', ms=7,
             label=f'saturated target ({saturated.sum()})')
    mark_transits()
    plt.xlabel(TIME_LABEL)
    plt.ylabel(r'$\Delta m$ - median')
    plt.title(f'{title} – {os.path.basename(folder[:-1])} – with saturated frames')
    plt.legend()
    finish_figure(folder + 'rel_lightcurve_with_saturated.png')

    plt.figure(figsize=(10, 4))
    plt.errorbar(time, mag[:, 0] - np.nanmedian(mag[:, 0]), yerr=t_err, fmt='--o', capsize=2, label='target')
    plt.errorbar(time, mag_ens, yerr=mag_ens_err, fmt='--o', capsize=2, label='ensemble')
    plt.xlabel(TIME_LABEL)
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
    plt.xlabel(TIME_LABEL)
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
    ax.plot(time, dm, '.', alpha=0.4, label='data')
    ax.plot(b['time'], b['y'], 'o-', label=f'binned ({POINTS_PER_BIN} pts)')
    mark_transits(ax)
    ax.set_title(mode, fontsize=10)
    ax.set_ylabel(r'$\Delta m$')
    ax.legend(fontsize=8, loc='lower left')
axes[-1, 0].set_xlabel(TIME_LABEL)
fig.suptitle(f'{title} – photometry modes')
fig.tight_layout()
finish_figure(os.path.join(results_root, 'modes_comparison.png'))
print(f"comparison saved in {results_root}/modes_comparison.png")
