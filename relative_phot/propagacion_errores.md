# Propagación de errores en fotometría diferencial

## Contexto

En fotometría diferencial se calcula el flujo relativo de una estrella objetivo respecto
a la mediana de un conjunto de estrellas de comparación. El objetivo es cancelar
variaciones atmosféricas comunes a todas las estrellas.

---

## 1. Flujo relativo

El flujo relativo normalizado es:

```
r = f_objetivo / (N · m)
```

donde:
- `f_objetivo` — flujo medido de la estrella objetivo en cada época
- `m` — mediana de los flujos de las estrellas de comparación en cada época
- `N` — constante de normalización: `N = median(f_objetivo / m)` (hace que r ≈ 1 fuera del tránsito)

`N` es una constante derivada de los datos, por lo que no aporta incertidumbre.

---

## 2. Regla de propagación para un cociente

Para una cantidad `r = A / B`, la propagación de errores estándar (primer orden) da:

```
σ_r / r = √( (σ_A / A)² + (σ_B / B)² )
```

Aplicado aquí:

```
σ_r / r = √( (σ_f / f_objetivo)² + (σ_m / m)² )
```

```
σ_r = r · √( (σ_f / f_objetivo)² + (σ_m / m)² )
```

---

## 3. Incertidumbre de la mediana de comparación (σ_m)

La mediana de `n` valores independientes con errores distintos `σ_i` tiene una
incertidumbre aproximada de:

```
σ_m = √(π/2) · √(Σ σ_i²) / n
```

El factor `√(π/2) ≈ 1.253` viene de la eficiencia estadística de la mediana:
para distribuciones gaussianas, la varianza del estimador mediana es π/2 veces
la varianza del estimador media con el mismo número de muestras. Es decir, la
mediana es ~25% menos eficiente que la media.

En código:

```python
def epoch_comp_stats(group):
    n = len(group)
    return pd.Series({
        'mag_median':        group['mag'].median(),
        'flux_median':       group['flux'].median(),
        'sigma_mag_median':  np.sqrt(np.pi / 2) * np.sqrt((group['errmag']**2).sum()) / n,
        'sigma_flux_median': np.sqrt(np.pi / 2) * np.sqrt((group['flux_err']**2).sum()) / n,
    })
```

---

## 4. Propagación para la diferencia de magnitudes

Para la diferencia `Δm = m_estrella - m_mediana`, la propagación es simplemente
la suma en cuadratura (la resta de variables independientes):

```
σ_Δm = √( σ_estrella² + σ_mediana² )
```

En código:

```python
rel['errmag'] = np.sqrt(ref_star['errmag']**2 + medians['sigma_mag_median']**2)
```

---

## 5. Resumen de fórmulas

| Cantidad            | Fórmula de error                                                              |
|---------------------|-------------------------------------------------------------------------------|
| Mediana de comp.    | `σ_m = √(π/2) · √(Σ σᵢ²) / n`                                               |
| Diferencia de mag   | `σ_Δm = √( σ_estrella² + σ_mediana² )`                                        |
| Flujo relativo      | `σ_r = r · √( (σ_f/f)² + (σ_m/m)² )`                                         |

---

## 6. Implementación en `fits_to_table.py`

```python
# mediana de comparación por época (excluyendo el objetivo)
epoch_stats = df_comp.groupby('epoch').apply(epoch_comp_stats)

# diferencia de magnitud con error propagado
rel['mag'] = ref_star['mag'] - epoch_stats['mag_median']
rel['errmag'] = np.sqrt(ref_star['errmag']**2 + epoch_stats['sigma_mag_median']**2)

# flujo relativo con error propagado
dif['relative_flux'] = df_star['flux'] / epoch_stats['flux_median']
norm = np.median(dif['relative_flux'])
dif['relative_flux'] /= norm

dif['relative_err'] = dif['relative_flux'] * np.sqrt(
    (df_star['flux_err'] / df_star['flux'])**2 +
    (epoch_stats['sigma_flux_median'] / epoch_stats['flux_median'])**2
)
```

---

## 7. Notas importantes

- La mediana de comparación **no debe incluir la estrella objetivo**. Si se incluye,
  se introduce una correlación espuria entre numerador y denominador.
- `σ_f` y `σ_m` se asumen independientes (lo cual es válido si las estrellas
  están suficientemente separadas y el ruido dominante es fotónico o de lectura).
- La fórmula de `σ_m` con el factor `√(π/2)` es una aproximación válida para
  errores similares entre las estrellas de comparación. Si los errores son muy
  heterogéneos, convendría usar una mediana ponderada.
