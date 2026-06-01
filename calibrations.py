import os
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
from glob import glob


def sort_by_date(fits_paths):
    date_map = {}
    for path in fits_paths:
        with fits.open(path) as hdul:
            date_map[path] = Time(hdul[0].header['DATE-OBS'], format='isot').jd

    return sorted(fits_paths, key=lambda p: date_map[p])


def create_master_dark(dark_dir, output_file='master_dark.fits', extension='.fits'):
    # dark directory
    dark_files = glob(os.path.join(dark_dir, '*'+extension))
    
    if not dark_files:
        raise ValueError("no files in the dark directory")
    
    print(f"processing {len(dark_files)} darks...")
    
    # initialization
    with fits.open(dark_files[0]) as hdul:
        header_first = hdul[0].header
        exptime_first = header_first['EXPTIME']
    
    # append darks
    darks = []
    for dark_file in dark_files:
        print(f"processing dark: {dark_file}")
        with fits.open(dark_file) as hdul:
            exptime = hdul[0].header['EXPTIME']
            if exptime != exptime_first:
                raise ValueError("not all darks have the same exptime.")
            else:
                darks.append(hdul[0].data)

    # masterize
    darks = np.array(darks)
    master_dark = np.median(darks, axis=0).astype(np.uint16)
    
    # save
    hdu = fits.PrimaryHDU(master_dark, header=header_first)
    hdu.writeto(output_file, overwrite=True)
    print(f"master dark saved as: {output_file}")


def create_master_flat(flat_dir, master_dark=None, output_file='master_flat.fits', extension='.fit',
                      saturation_threshold=60000):
    
    # flat directory
    flat_files = glob(os.path.join(flat_dir, '*'+extension))
    
    if not flat_files:
        raise ValueError("no files in the flat directory")
    
    print(f"processing {len(flat_files)} flats...")

    # initialization
    with fits.open(flat_files[0]) as hdul:
        data_shape = hdul[0].data.shape
        header_first = hdul[0].header
        exptime_first = header_first['EXPTIME'] 
        filter_first = header_first['FILTER']
    
    print(f" processing filter {filter_first}...")

    # dark
    if master_dark is not None:
        with fits.open(master_dark) as hdul:
            master_dark = hdul[0].data
            master_dark_header = hdul[0].header
            exptime_dark = master_dark_header['EXPTIME']

        if exptime_dark != exptime_first:
            raise ValueError('master dark have not the same exptime with flats')    
    else: 
        master_dark = 0
        print('No Master Dark')

    # append flats
    flats = []
    
    for flat_file in flat_files:
        print(f"processing flat: {flat_file}")
        with fits.open(flat_file) as hdul:
            flat_header = hdul[0].header
            if flat_header['FILTER'] != filter_first:
                print(f"flat {flat_file} have not the same filter")
            else:
                flat_data = hdul[0].data - master_dark
                mean, median, std = sigma_clipped_stats(flat_data)
                if median > saturation_threshold:
                    print(f"flat {flat_file} discarded due to saturation")
                else: 
                    flats.append(flat_data / mean)
    
    # masterize
    flats = np.array(flats)
    master_flat = np.median(flats, axis=0)
    master_flat = master_flat / np.mean(master_flat)
    
    # save
    hdu = fits.PrimaryHDU(master_flat, header=header_first)
    hdu.writeto(output_file, overwrite=True)
    print(f"master flat saved as: {output_file}")


def extract_subframe(image, center, size, save=None):
    cx, cy = center
    sub = image[cy - size:cy + size, cx - size:cx + size]
    if save is not None:
        fits.writeto(save, sub, overwrite=True)
        print(f"subframe saved as: {save}")
    return sub


def coords_to_subframe(coords, center, size):
    """Transform a list of (x, y) coords from the original image to the subframe."""
    cx, cy = center
    x0, y0 = cx - size, cy - size
    return [(x - x0, y - y0) for x, y in coords]


def calibrate(data_fits, master_flat, master_dark, isbias=False, isfit=True):
    if isfit:
        data = data_fits.data

        if master_flat is None:
            if master_dark is None:
                return data

        header_data = data_fits.header
        exptime_data = header_data['EXPTIME']
        filter_data = header_data['FILTER']

        # dark
        header_dark = master_dark.header
        exptime_dark = header_dark['EXPTIME']
        dark = master_dark.data

        if exptime_dark != exptime_data:
            if not isbias:
                raise ValueError('dark and data have not the same exptime')
            else:
                print('using bias instead of dark...')

        # flat
        header_flat = master_flat.header
        filter_flat = header_flat['FILTER']
        flat = master_flat.data

        if filter_flat != filter_data:
            raise ValueError('flat and data have not the same filter')
        
    else:
        data = data_fits
        dark = master_dark
        flat = master_flat

    return (data - dark) / flat