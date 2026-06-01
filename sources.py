import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import detect_sources, SourceCatalog


def extract_sources(image, threshold_sigma=5.0, npixels=10, saturation_threshold=60000,
                    box_size=(50, 50), filter_size=(3, 3), sigma=3.0,
                    sort_by='peak', edge_margin=0, output_file=None):
    """
    Extract sources from a 2D image array using photutils segmentation.

    Parameters
    ----------
    image : np.ndarray
        2D calibrated image.
    threshold_sigma : float
        Detection threshold in background RMS units.
    npixels : int
        Minimum connected pixels to count as a source.
    saturation_threshold : float or None
        Discard sources whose peak pixel >= this value. None keeps all.
    box_size : tuple
        Box size for 2D background estimation.
    filter_size : tuple
        Median filter size for background smoothing.
    sigma : float
        Sigma for background sigma-clipping.
    sort_by : str
        'peak'  → sort by peak pixel value.
        'flux'  → sort by total segment flux.
    edge_margin : int
        Discard sources whose centroid is within this many pixels of any border.
        0 (default) keeps all sources.
    output_file : str or None
        If given, save the catalog as CSV at this path.

    Returns
    -------
    pd.DataFrame
        Columns: id, x, y, peak, flux, area, ellipticity, fwhm
        Sorted strongest → weakest. Empty DataFrame if no sources found.
    """
    sigma_clip = SigmaClip(sigma=sigma)
    bkg = Background2D(
        image,
        box_size=box_size,
        filter_size=filter_size,
        sigma_clip=sigma_clip,
        bkg_estimator=MedianBackground(),
    )
    image_sub = image - bkg.background
    threshold = threshold_sigma * bkg.background_rms

    segm = detect_sources(image_sub, threshold=threshold, npixels=npixels)
    if segm is None:
        print("no sources detected")
        return pd.DataFrame()

    cat = SourceCatalog(image_sub, segm)

    df = pd.DataFrame({
        'id':          np.asarray(cat.label),
        'x':           np.asarray(cat.xcentroid, dtype=float),
        'y':           np.asarray(cat.ycentroid, dtype=float),
        'peak':        np.asarray(cat.max_value, dtype=float),
        'flux':        np.asarray(cat.segment_flux, dtype=float),
        'area':        np.asarray(cat.area, dtype=int),
        'ellipticity': np.asarray(cat.ellipticity, dtype=float),
        'fwhm':        np.asarray(cat.fwhm, dtype=float),
    })

    if saturation_threshold is not None:
        n_before = len(df)
        df = df[df['peak'] < saturation_threshold].reset_index(drop=True)
        removed = n_before - len(df)
        if removed:
            print(f"{removed} saturated sources removed (peak >= {saturation_threshold})")

    if edge_margin > 0:
        nrows, ncols = image.shape
        mask = (
            (df['x'] >= edge_margin) &
            (df['x'] <= ncols - 1 - edge_margin) &
            (df['y'] >= edge_margin) &
            (df['y'] <= nrows - 1 - edge_margin)
        )
        n_before = len(df)
        df = df[mask].reset_index(drop=True)
        removed = n_before - len(df)
        if removed:
            print(f"{removed} edge sources removed (margin={edge_margin} px)")

    sort_col = 'peak' if sort_by == 'peak' else 'flux'
    df = df.sort_values(sort_col, ascending=False).reset_index(drop=True)
    print(f"{len(df)} sources extracted")

    if output_file is not None:
        df.to_csv(output_file, index=False)
        print(f"catalog saved to {output_file}")

    return df


def sources_from_fits(fits_path, hdu=0, **kwargs):
    """
    Load a FITS file and extract sources.

    Parameters
    ----------
    fits_path : str
        Path to the FITS file.
    hdu : int
        HDU index to read. Default 0 (primary).
    **kwargs
        Forwarded to extract_sources.

    Returns
    -------
    pd.DataFrame  (same as extract_sources)
    """
    with fits.open(fits_path) as hdul:
        image = hdul[hdu].data.astype(float)
    print(f"loaded {fits_path}  shape={image.shape}")
    return extract_sources(image, **kwargs)


def brightest_xy(catalog, n, sort_by='peak'):
    """
    Return XY positions of the N brightest sources in a catalog.

    Parameters
    ----------
    catalog : pd.DataFrame
        Output of extract_sources — must contain 'x', 'y', and sort_by columns.
    n : int
        Number of sources to return.
    sort_by : str
        Column used to rank brightness: 'peak' or 'flux'. Default 'peak'.

    Returns
    -------
    list of [float, float]
        [[x0, y0], [x1, y1], ...] sorted brightest → faintest.
    """
    col = 'peak' if sort_by == 'peak' else 'flux'
    top = catalog.nlargest(n, col)
    return top[['x', 'y']].values.tolist()
