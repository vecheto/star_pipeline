import photutils
import numpy as np
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from photutils.aperture import CircularAperture, CircularAnnulus
from photutils.aperture import aperture_photometry
from astropy.stats import sigma_clipped_stats
from photutils.segmentation import detect_sources
from photutils.segmentation import SourceCatalog
import pandas as pd
import matplotlib.pyplot as plt
from photutils.psf import fit_fwhm
from photutils.centroids import centroid_2dg
import numpy as np

def cut_stamp(frames, center, radius, recenter=True):
    """
    frames: (N, Y, X)
    center: (x, y) aproximado en coordenadas del frame
    radius: radio de la stamp
    """

    x0, y0 = center

    # Corte inicial (grande)
    stamps = frames[
        :,
        y0-radius:y0+radius,
        x0-radius:x0+radius
    ]

    if not recenter:
        return stamps

    # Usar la mediana para encontrar mejor centroide
    median_stamp = np.median(stamps, axis=0)

    try:
        cx, cy = centroid_2dg(median_stamp)
    except Exception:
        # fallback si falla el ajuste
        cx, cy = np.array(median_stamp.shape[::-1]) / 2

    # Centro refinado en coordenadas del frame
    x_ref = int(x0 - radius + cx)
    y_ref = int(y0 - radius + cy)

    # Re-corte centrado correctamente
    stamps_refined = frames[
        :,
        y_ref-radius:y_ref+radius,
        x_ref-radius:x_ref+radius
    ]

    return stamps_refined


def centroid(stamp):
    cx, cy = centroid_2dg(stamp)
    return cx, cy



def aperture_phot(stamp, centro, radio, skyradio_int, skyradio_ext):
    aperture = CircularAperture(centro, r=radio)
    annulus = CircularAnnulus(centro, r_in=skyradio_int, r_out=skyradio_ext)
    
    # APERTURE ON SOURCE
    phot_table = aperture_photometry(stamp, [aperture, annulus])
    
    # SKY FLUX EXTIMATION
    annulus_mask = annulus.to_mask(method='center')
    annulus_data = annulus_mask.multiply(stamp)
    sky_pixels = annulus_data[annulus_data != 0]
    sky_mean, sky_median, sky_std = sigma_clipped_stats(sky_pixels, sigma=3.0)
    sky_sum = sky_median * aperture.area
    
    # FLUX - BACKGROUND
    flux_source = phot_table['aperture_sum_0'][0]
    net_flux = flux_source - sky_sum
    
    # ERROR
    # Poisson in source (sqrt(N_counts))
    error_source = np.sqrt(flux_source)
    
    # Error in sky (ruido por píxel * sqrt(n_pix en apertura))
    error_sky = sky_std * np.sqrt(aperture.area)
    
    # Error total (suma en cuadratura)
    net_error = np.sqrt(error_source**2 + error_sky**2)
    
    return net_flux, net_error


def estimate_background(image, box_size=(50, 50), filter_size=(3, 3), sigma=3.0):
    sigma_clip = SigmaClip(sigma=sigma)
    bkg_estimator = MedianBackground()
    bkg = Background2D(
        image,
        box_size=box_size,
        filter_size=filter_size,
        sigma_clip=sigma_clip,
        bkg_estimator=bkg_estimator
    )
    return bkg.background, bkg.background_rms


def extract_sources(image, background, background_rms, threshold_sigma=5.0, npixels=10):
    image_sub = image - background
    threshold = threshold_sigma * background_rms
    segm = detect_sources(
        image_sub,
        threshold=threshold,
        npixels=npixels
    )
    
    # Crear catálogo
    cat = SourceCatalog(image_sub, segm)
    #catalog = cat.to_table().to_pandas() 
    return segm, cat

# TODO
def psf_phot():
    return

def optimize_parameters(stamp, fwhm_min=1, fwhm_max=10, step=0.1, min_sky_width=1):

    # --- asegurar stamp impar ---
    ny, nx = stamp.shape
    if ny % 2 == 0:
        stamp = stamp[:-1, :]
    if nx % 2 == 0:
        stamp = stamp[:, :-1]

    # --- centro y FWHM ---
    centro = centroid(stamp)
    fwhm = fit_fwhm(stamp)

    # --- grilla de búsqueda (escalada con FWHM) ---
    ap = np.arange(fwhm_min*fwhm, fwhm_max*fwhm, step)
    sky_in = np.arange(fwhm_min*fwhm + step, fwhm_max*fwhm, step)
    sky_out = np.arange(fwhm_min*fwhm + 2*step, fwhm_max*fwhm, step)

    AP, IN, OUT = np.meshgrid(ap, sky_in, sky_out, indexing="ij")
    params = np.column_stack([AP.ravel(), IN.ravel(), OUT.ravel()])

    min_sky_width_fwhm = min_sky_width * fwhm

    mask = (
        (params[:, 0] < params[:, 1]) &
        (params[:, 1] < params[:, 2]) &
        ((params[:, 2] - params[:, 1]) >= min_sky_width_fwhm)
    )

    params = params[mask]

    print(f"Iniciando optimización con {len(params)} combinaciones")

    # --- evaluar SNR ---
    snr_values = np.empty(len(params))
    for i, (ap_, inn_, out_) in enumerate(params):
        print(f"Probando combinacion {i} : {ap_, inn_, out_}")
        try:
            net_flux, net_error = aperture_phot(stamp, centro, ap_, inn_, out_)
            snr_values[i] = snr(net_flux, net_error)
        except Exception:
            snr_values[i] = np.nan

    best_idx = np.nanargmax(snr_values)

    best_params = params[best_idx]
    best_snr = snr_values[best_idx]

    # --- diccionario de salida ---
    result = {
        "best_params": {
            "aperture": best_params[0],
            "sky_in": best_params[1],
            "sky_out": best_params[2],
        },
        "best_snr": best_snr,
        "fwhm": fwhm,
        "centroid": centro,
        "params": params,
        "snr_values": snr_values,
    }

    return result   

def plot_optimization(res, path):
    params = res["params"]
    snr = res["snr_values"]
    best = res["best_params"]
    fwhm = res['fwhm']

    # Aperture
    unique_ap = np.unique(params[:, 0])
    snr_max_ap = [
        np.nanmax(snr[params[:, 0] == a])
        for a in unique_ap
        ]

    # Sky in
    unique_si = np.unique(params[:, 1])
    snr_max_si = [
        np.nanmax(snr[params[:, 1] == a])
        for a in unique_si
        ]
    
    # Sky out
    unique_so = np.unique(params[:, 2])
    snr_max_so = [
        np.nanmax(snr[params[:, 2] == a])
        for a in unique_so
        ]
    
    plt.plot(unique_ap/fwhm, snr_max_ap, label='Aperture', color='green')
    plt.plot(unique_si/fwhm, snr_max_si, label='Sky In', color='blue')
    plt.plot(unique_so/fwhm, snr_max_so, label='Sky Out', color='red')
    
    plt.axvline(best["aperture"]/fwhm, color="green", ls="--")
    plt.axvline(best["sky_in"]/fwhm, color="blue", ls="--")
    plt.axvline(best["sky_out"]/fwhm, color="red", ls="--")
    plt.axhline(np.nanmax(snr), ls="--", label='Max SNR')

    plt.legend()
    plt.xlabel("FHWM")
    plt.ylabel("Max SNR")
    plt.savefig(path+'optimization.png')
    plt.show()


def snr(net_flux, net_error):
    if net_error>0:
        snr_value = net_flux/net_error
    else:
        snr_value = np.nan       
    return snr_value


def plot_image(image,title, nsigma=1):
    median = np.median(image)
    std = np.std(image)
    plt.imshow(image, 
           vmin=median-nsigma*std,
           vmax=median+nsigma*std,
           cmap='inferno',
           origin='lower')
    plt.title(title)
    plt.show()