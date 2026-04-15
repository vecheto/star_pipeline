# Star Photometry Pipeline — CLAUDE.md

## Project Overview

Astronomical aperture photometry pipeline for processing stellar FITS images and producing light curves. Supports variable star monitoring, exoplanet transit detection, and multi-filter photometry. All observations are stored as multi-frame FITS files; the pipeline calibrates, aligns, and measures stellar flux across frames.

---

## Pipeline Entry Points

| Script | Target | Data Dir | Output |
|--------|--------|----------|--------|
| `pipeline_variable_star.py` | RRPic (variable star) | `estrellavariable/light/` | `RRPic/fluxes.csv` |
| `pipeline_variable_star_copy.py` | RRPic_2 | `croppedvariablestars/` | `RRPic_2/fluxes.csv` |
| `pipeline_transit.py` | CoRoT-7 (exoplanet transit) | `COROT7/` | `corot-7/fluxes.csv` |
| `plot_transit.py` | Post-processing | `RRPic_2/fluxes.csv` | `relative_flux*.png` |
| `fits_to_table.py` | J0559-1404 (astrometry) | `J0559-1404_NTT_sequence.fits` | Light curve plots |
| `hr_diagram_example.py` | HR diagram template | `omega/B/` | (incomplete — TODO) |

**Typical run order for a variable star:**
```
python pipeline_variable_star.py   # calibrate → align → photometry → fluxes.csv
python plot_transit.py             # differential photometry → light curve plots
```

---

## Module Architecture

```
calibrations.py     ─┐
alignement.py        ├─→  pipeline_*.py  →  plot_transit.py
photometry.py       ─┘                   →  fits_to_table.py
```

### `calibrations.py`
Standard CCD reduction.

- `create_master_dark(dark_dir)` — median-combines dark frames; checks uniform exposure time
- `create_master_flat(flat_dir, master_dark=None)` — normalises and median-combines flats; rejects saturated frames (>60 000 counts)
- `calibrate(data_fits, master_flat, master_dark, subframe, isbias, use_subframe)` — applies `(data − dark) / flat`; supports subframe extraction

### `alignement.py`  *(note: typo in filename is intentional)*
Sub-pixel image registration and stacking.

- `phase_correlation_alignment(images, upsamplefactor=10)` — Fourier-space cross-correlation; first frame is reference; applies fractional-pixel shift in frequency domain
- `stack(images, method)` — combine aligned frames: `'mean'`, `'median'`, `'max'`, or `'sum'`

### `photometry.py`
Core measurement engine.

- `cut_stamp(frames, center, radius, recenter=True)` — extracts N×(2r+1)×(2r+1) postage stamps; optionally recenters using 2D Gaussian centroid
- `centroid(stamp)` — calls `photutils.centroids.centroid_2dg`
- `aperture_phot(stamp, centro, radio, skyradio_int, skyradio_ext)` — circular aperture + annular sky; sky estimated with sigma-clipped median; returns `(net_flux, net_error)` with full error propagation (Poisson + sky noise)
- `optimize_parameters(stamp, fwhm_min, fwhm_max, step, min_sky_width)` — grid search over aperture and sky radii scaled by FWHM; selects combination maximising SNR; returns `best_params`, `best_snr`, all tested values
- `snr(net_flux, net_error)` — `flux / error`
- `estimate_background(image, box_size, filter_size, sigma)` — 2D background map + RMS using median with sigma clipping
- `extract_sources(image, background, background_rms, threshold_sigma, npixels)` — connected-component source detection; returns segmentation map + catalog
- `plot_optimization(res, path)` — visualises SNR vs aperture/sky parameters; marks optimum
- `plot_image(image, title, nsigma=1)` — auto-scaled display (median ± n·σ)

---

## Data Flow

```
Raw FITS (multi-frame)
  │
  ├─ calibrations.calibrate()          ← dark subtraction + flat division
  │
  ├─ alignement.phase_correlation_alignment()   ← sub-pixel registration
  │
  ├─ photometry.cut_stamp()            ← extract stamps per star
  │
  ├─ photometry.optimize_parameters()  ← grid-search best aperture (runs once on median stamp)
  │
  ├─ photometry.aperture_phot()        ← measure flux + error per frame per star
  │
  └─ pandas DataFrame → fluxes.csv
       columns: time (JD), flux_star_N, err_flux_star_N

fluxes.csv
  └─ plot_transit.py
       ├─ differential photometry (target / median of comparisons)
       └─ binned light curve → *.png plots
```

---

## Key Parameters

| Parameter | Typical Value | Where Used |
|-----------|--------------|------------|
| Stamp radius | 15 px (RRPic), 18 px (transit/RRPic_2) | `cut_stamp` |
| FWHM search range | 0.5–5 or 0.5–6 | `optimize_parameters` |
| Saturation limit | 60 000 counts | flat rejection, stamp check |
| Background box | 50×50 px, filter 3×3 | `estimate_background` |
| Background σ-clip | 3σ | background, sky estimation |
| Source threshold | 5–10σ | `extract_sources` |
| Min source pixels | 10–40 | `extract_sources` |
| Phase corr. upsample | ×10 | `phase_correlation_alignment` |

---

## Output Files

| File | Description |
|------|-------------|
| `{target}/fluxes.csv` | Time-series flux table (JD + flux + error per star) |
| `{target}/optimization.png` | SNR vs aperture/sky radii plot |
| `{target}/stamps/*.fits` | Per-star stamp FITS cubes (for CoRoT-7) |
| `relative_flux.png` | Differential light curve |
| `relative_flux_binned.png` | Binned light curve |
| `flux_star_N.png` | Individual star raw light curves |
| `df_sigma.csv` | Processed magnitude data (J0559-1404) |

---

## Data Directories

```
estrellavariable/light/    ← RRPic raw FITS
croppedvariablestars/      ← RRPic_2 raw FITS
COROT7/                    ← CoRoT-7 transit FITS
J0559-1404/                ← astrometric multi-epoch catalog
  ├─ J0559-1404_NTT_sequence.fits   (2.1 MB, multi-HDU star catalog)
  └─ J0559-1404_NTT_results.fits    (14.4 KB)
omega/B/, omega/G/, omega/R/        ← multi-filter OMEGA observations (B/G/R)
```

---

## Dependencies (`requirements.txt`)

```
astropy==7.0.1       # FITS I/O, WCS, time
photutils==2.2.0     # aperture photometry, centroids, background
numpy==1.26.0
scipy==1.15.2        # Fourier transforms, image shifts
scikit-image==0.25.2 # phase cross-correlation
pandas==2.2.3        # flux tables
matplotlib==3.10.1   # plots
astroquery==0.4.10   # (available, not yet used in pipeline scripts)
```

Virtual environment: `pipe/` (Python 3.12)

---

## Known Incomplete / TODO

- **HR diagram** (`hr_diagram_example.py`): calibration + alignment + source detection done; magnitude colour-colour plot not yet implemented
- **PSF photometry**: stub exists in photometry.py but not implemented
- **Master calibration creation**: `create_master_dark` / `create_master_flat` exist but pipeline scripts pass `None` for dark/flat (calibration currently skipped in RRPic pipelines)
- `alignement.py` — filename typo is inherited from the original code; do not rename without updating all imports

---

## `fits_to_table.py` — Special Case

This script reads a pre-existing astrometric catalog (multi-epoch FITS, not raw images). Each HDU is one epoch with per-star photometry already done externally. The script:
1. Extracts `EPOCH`, star `ID`, `X/Y`, `RA/DEC` (rad→deg), `mag`, `errmag`
2. Builds a long-format pandas DataFrame
3. Sigma-clips per-epoch magnitudes
4. Plots per-star and all-star light curves (magnitude deviation from epoch median, inverted Y axis per astronomical convention)

This is the entry point for **relative photometry from an external catalog** — different from the image-based pipelines above.
