# Star Photometry Pipeline — CLAUDE.md

## Descripción

Pipeline de fotometría de apertura para imágenes FITS crudas: calibra (dark/flat),
alinea con precisión sub-píxel, detecta fuentes, optimiza la apertura por SNR y produce
curvas de luz. Objetivo actual: el tránsito de **HATS-26**.

> El análisis de **J0559-1404** (fotometría relativa a partir de un catálogo astrométrico
> externo, sin reducción de imágenes) vive desde 2026-07-23 en su propio repo:
> `~/J0559-1404_phot`. Los dos proyectos son independientes y no comparten código.

---

## Entry points

| Script | Fotometría | Salida |
|--------|-----------|--------|
| `pipeline_variability.py` | Apertura, optimizada **una vez** sobre el stamp mediano de la primera estrella | `results/fluxes.csv` |
| `pipeline_variability_full.py` | Apertura, optimización **por estrella** (`fwhm_max=7`) | `results/fluxes.csv` |
| `pipeline_variability_isophot.py` | Isofotal (`isophotal_phot`), sin apertura fija | `results_isophot/isofluxes.csv` |
| `plot_results.py` | Post-proceso: fotometría diferencial + binning | `relative_flux*.png`, `stars_flux*.png` |
| `hr_diagram_example.py` | Plantilla de diagrama HR sobre `omega/B/` | incompleto — ver TODO |

Orden típico:

```bash
python pipeline_variability_full.py   # calibrar → alinear → fotometría → fluxes.csv
python plot_results.py                # fotometría diferencial → curvas de luz
```

---

## Rutas de datos (fuera del repo)

Los pipelines apuntan a directorios **hermanos** de este repo, no versionados:

```
../HATS-26_transit/          ← FITS crudos del tránsito; las salidas van a su results/
../RH200Calibs/              ← master_dark_100s.fits, master_dark_120s.fits, master_flat_L.fits
omega/B/, omega/G/, omega/R/ ← observaciones multi-filtro (6.9 GB, gitignored)
```

Ninguno de los dos primeros existe actualmente en `~`; hay que montarlos antes de correr
los pipelines de HATS-26.

---

## Arquitectura

```
calibrations.py ─┐
alignement.py    ├─→ pipeline_variability*.py → plot_results.py
photometry.py    │
sources.py      ─┘
```

### `calibrations.py`
- `sort_by_date(fits_paths)` — ordena por fecha de cabecera
- `create_master_dark(dark_dir, ...)` — combina darks por mediana; verifica tiempo de exposición uniforme
- `create_master_flat(flat_dir, master_dark=None, ...)` — normaliza y combina flats; descarta saturados (>60 000 cuentas)
- `extract_subframe(image, center, size, save=None)` / `coords_to_subframe(coords, center, size)`
- `calibrate(data_fits, master_flat, master_dark, isbias=False, isfit=True)` — aplica `(data − dark) / flat`

### `alignement.py` *(el typo del nombre es heredado — no renombrar sin actualizar los imports)*
- `phase_correlation_alignment(images, upsamplefactor=10)` — correlación cruzada en Fourier; la primera imagen es la referencia
- `stack(images, method)` — `'mean'`, `'median'`, `'max'` o `'sum'`

### `sources.py`
- `extract_sources(image, threshold_sigma=5.0, npixels=10, saturation_threshold=60000, ...)` — detección por componentes conexas; devuelve mapa de segmentación + catálogo
- `sources_from_fits(fits_path, hdu=0, **kwargs)`
- `brightest_xy(catalog, n, sort_by='peak')` — las `n` fuentes más brillantes

### `photometry.py`
- `cut_stamp(frames, center, radius, recenter=True)` — stamps N×(2r+1)×(2r+1); recentra con centroide gaussiano 2D
- `centroid(stamp)` — `photutils.centroids.centroid_2dg`
- `aperture_phot(stamp, centro, radio, skyradio_int, skyradio_ext)` — apertura circular + anillo de cielo (mediana con sigma-clip); devuelve `(net_flux, net_error)` con propagación completa (Poisson + ruido de cielo)
- `isophotal_phot(stamp)` — flujo isofotal; devuelve `(flux, pixel_std)`
- `optimize_parameters(stamp, fwhm_min=1, fwhm_max=10, step=0.1, min_sky_width=1)` — grid search de radios escalados por FWHM; maximiza SNR
- `plot_optimization(res, path)` — SNR vs radios, marca el óptimo
- `snr(net_flux, net_error)` — `flux / error`
- `plot_image(image, title, nsigma=1)` — escalado automático (mediana ± n·σ)
- `psf_phot()` — **stub, sin implementar**

---

## Flujo de datos

```
FITS crudos (multi-frame)
  ├─ calibrations.calibrate()                    ← dark + flat
  ├─ alignement.phase_correlation_alignment()    ← registro sub-píxel
  ├─ sources.extract_sources()                   ← detección de estrellas
  ├─ photometry.cut_stamp()                      ← stamp por estrella
  ├─ photometry.optimize_parameters()            ← mejor apertura por SNR
  ├─ photometry.aperture_phot()                  ← flujo + error por frame y estrella
  └─ pandas → fluxes.csv   (time [JD], flux_star_N, err_flux_star_N)

fluxes.csv → plot_results.py
  ├─ flujo relativo = objetivo / mediana(comparación), normalizado
  ├─ error propagado: σ_m = √(π/2)·√(Σσᵢ²)/n
  └─ curvas binned + RMS
```

---

## Parámetros clave

| Parámetro | Valor | Dónde |
|-----------|-------|-------|
| `target_star` (HATS-26) | `[2242, 1792]` (X, Y) | pipelines |
| Rango FWHM | 0.5–5 (`variability`), 0.5–7 (`_full`) | `optimize_parameters` |
| Límite de saturación | 60 000 cuentas | flats, `extract_sources` |
| σ-clip | 3σ | fondo y cielo |
| Umbral de detección | 5σ | `extract_sources` |
| Píxeles mínimos por fuente | 10 | `extract_sources` |
| Upsample correlación de fase | ×10 | `phase_correlation_alignment` |
| Binning | `len(fluxes)/4` bins | `plot_results.py` |

---

## Salidas (en `../HATS-26_transit/results/`)

| Archivo | Contenido |
|---------|-----------|
| `fluxes.csv` | Serie temporal: JD + flujo + error por estrella |
| `isofluxes.csv` | Ídem, fotometría isofotal (`results_isophot/`) |
| `star_N_optimization.png` | SNR vs radios de apertura/cielo |
| `N_aperture_check.png` | Verificación visual de la apertura |
| `stamps/stamp*.fits` | Cubos de stamps por estrella |
| `relative_flux.png`, `relative_flux_binned.png` | Curva de luz diferencial |
| `stars_flux.png`, `stars_flux_miuns_median.png` | Flujos crudos y residuos |

---

## Dependencias

```
astropy==7.0.1       photutils==2.2.0     numpy==1.26.0
scipy==1.15.2        scikit-image==0.25.2 pandas==2.2.3
matplotlib==3.10.1   astroquery==0.4.10
```

Entorno virtual: `pipe/` (Python 3.12). `requirements.txt` es un `pip freeze` completo.

---

## Pendiente / TODO

- **Diagrama HR** (`hr_diagram_example.py`): calibración, alineamiento y detección listos;
  falta el gráfico color-color. Además apunta a rutas (`calibration/darks/…`) que no existen.
- **Fotometría PSF**: `psf_phot()` es un stub vacío.
- Los tres `pipeline_variability*.py` comparten ~80 % del código (inicialización, calibración,
  alineamiento, detección) y solo divergen en el método fotométrico. Si se tocan, conviene
  factorizar la parte común antes de seguir duplicando.
- `create_master_dark` / `create_master_flat` existen, pero los pipelines cargan masters
  ya hechos desde `../RH200Calibs/`.
