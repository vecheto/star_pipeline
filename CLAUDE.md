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

| Script | Qué hace |
|--------|----------|
| `pipeline_variability.py` | Pipeline único: calibra → alinea → plate solve → fotometría. Modo vía `photometry_mode` |
| `pipeline_params.py` | Plantilla de parámetros; se copia a la carpeta de datos y se edita |
| `plot_results.py` | Post-proceso: fotometría diferencial + binning (ruta a `fluxes.csv` aún hardcodeada) |
| `hr_diagram_example.py` | Plantilla de diagrama HR sobre `omega/B/` (incompleta) |
| `legacy/` | Los tres pipelines anteriores (`variability`, `_full`, `_isophot`), solo referencia |

`photometry_mode`:
- `'aperture_common'` — apertura optimizada **una vez** sobre el objetivo y aplicada a todas las estrellas (misma fracción de PSF → robusto a cambios de seeing en el cociente)
- `'aperture_per_star'` — apertura optimizada por estrella
- `'isophotal'` — `isophotal_phot`; la columna de error es `pixel_std_star_N`, no un error de flujo

Uso:

```bash
cd ../HATS-26_transit                 # carpeta con los FITS + pipeline_params.py editado
~/VICENTE/star_pipeline/pipe/bin/python ~/VICENTE/star_pipeline/pipeline_variability.py pipeline_params.py
```

Después, `.../plot_results.py pipeline_params.py` hace la fotometría relativa con el método
ensemble de `../brown_dwarfs/metodo.md` ([relative.py](relative.py)): referencias medidas en
todos los frames; δᵢ(t) = mᵢ − mediana_t(mᵢ); sigma clip por frame (|δᵢ − media| > `clip_sigma`·std,
iterativo, punto a punto, como `brown_dwarfs/calibrated.zero_point`); ensemble = media ponderada
1/σᵢ²; Δm = m★ − δ_ens. Con N ≤ 10 referencias el clip a 3σ no rechaza nada → `n_refs` ≥ ~20.
Excluye los frames con el objetivo saturado (`n_sat_target` > 0 en `fluxes.csv`).
Escribe en cada `results/<modo>/`: `relative_phot.csv`, `star_selection.txt`, `rel_lightcurve*.png`
(`rel_lightcurve_with_saturated.png` muestra además los frames saturados con ×),
`rel_target_vs_ensemble.png`, `rel_lightcurves_refstars.png` (cada referencia contra un ensemble
sin ella: identifica referencias variables); y `results/modes_comparison.png`. Gráficos en
magnitudes con eje Y invertido.

`flat = None` → flat = 1; `darks = None` → dark = 0.
Rutas relativas en `pipeline_params.py` se resuelven desde su carpeta. Salidas en
`<carpeta>/results/<photometry_mode>/`. Con `show_plots = False` las figuras solo se guardan.
Referencias: no saturadas (pico crudo), aisladas (sin otra detección a < 1.5·`stamp_size`) y lejos del objetivo.
Se detectan sobre el **máximo** del stack alineado: se descartan estrellas que
saturan en cualquier frame. `stamp_size` debe ser ≳ `fwhm_max·FWHM`; si no,
`optimize_parameters` recorta los radios al borde del stamp y avisa.

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

### `plate_solve.py` *(no llamarlo `astrometry.py`: tapa el paquete Python de astrometry.net y `solve-field` falla)*
- `header_pointing(header)` — RA/Dec del header (`RAJ2000`/`DECJ2000`, `RA`/`DEC`, `OBJCTRA`/`OBJCTDEC`); `None` si no hay
- `header_pixel_scale(header)` — ″/px desde `FOCALLEN`, `XPIXSZ`, `XBINNING`
- `solve_wcs(image, center, radius, scale, timeout)` — `solve-field` local (índices Tycho-2 07–19 en `/usr/share/astrometry`); timeout por reloj matando el grupo de procesos (`--cpulimit` no se respeta)
- `target_pixel(wcs, target)` — (x, y) 0-based de un nombre resoluble por Sesame o un `SkyCoord`

Orden en el pipeline: `target_xy` manual → plate solve centrado en el **objetivo** → centrado
en el apuntado del header → fuente más cercana al centro. El objetivo va primero porque el
apuntado del header puede estar mal (VW Hyi, nov 2026: montura reportaba 4.8° fuera).
Siempre revisar `results/<modo>/finder.png`.

### `sources.py`
- `extract_sources(image, threshold_sigma=5.0, npixels=10, saturation_threshold=60000, ...)` — detección por componentes conexas; devuelve mapa de segmentación + catálogo
- `sources_from_fits(fits_path, hdu=0, **kwargs)`
- `brightest_xy(catalog, n, sort_by='peak')` — las `n` fuentes más brillantes

### `photometry.py`
- `cut_stamp(frames, center, radius, recenter=True)` — stamps N×(2r+1)×(2r+1); recentra con centroide gaussiano 2D
- `centroid(stamp)` — `centroid_2dg` sobre el stamp con cielo restado (photutils 2.x no ajusta fondo); si cae fuera, `centroid_com`
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
| Posición del objetivo | por plate solve desde `target_name` | `pipeline_params.py` |
| Rango FWHM | `fwhm_min`–`fwhm_max` (0.5–7) | `pipeline_params.py` |
| Límite de saturación | 60 000 cuentas | flats, `extract_sources` |
| σ-clip | 3σ | fondo y cielo |
| Umbral de detección | 5σ | `extract_sources` |
| Píxeles mínimos por fuente | 10 | `extract_sources` |
| Upsample correlación de fase | ×10 | `phase_correlation_alignment` |
| Binning | `len(fluxes)/4` bins | `plot_results.py` |

---

## Salidas (en `<carpeta de datos>/results/<photometry_mode>/`)

| Archivo | Contenido |
|---------|-----------|
| `fluxes.csv` | JD + `exptime` + `flux_star_N` + `err_flux_star_N` (o `pixel_std_star_N` en isofotal) |
| `exptime_vs_jd.png` | Tiempo de exposición vs JD |
| `stars.csv` | Posición (x, y) de cada estrella en el frame de referencia; 0 = objetivo |
| `star_N_optimization.png` | SNR vs radios de apertura/cielo |
| `N_aperture_check.png` | Verificación visual de la apertura |
| `N_max_stamp.png` | Máximo del stamp en el tiempo, con píxeles saturados |
| `stamps/star_N.fits` | Cubo de stamps por estrella |
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
- `plot_results.py` aún lee `../HATS-26_transit/results/fluxes.csv`; adaptarlo a `results/<modo>/` y al archivo de parámetros.
- `create_master_dark` / `create_master_flat` existen, pero los pipelines cargan masters
  ya hechos desde `../RH200Calibs/`.
