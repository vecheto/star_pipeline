import os
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from glob import glob

def create_master_dark(dark_dir, output_file='master_dark.fits', extension='.fit'):
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
    master_dark = np.median(darks, axis=0)
    
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


def calibrate(data_fits, master_flat, master_dark, subframe=[[0,1000],[0,1000]], isbias=False, use_subframe=True):
    # data
    if use_subframe:
        data = data_fits.data[subframe[1][0]:subframe[1][1], subframe[0][0]:subframe[0][1]]
    else:
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
    if use_subframe:
        dark = master_dark.data[subframe[1][0]:subframe[1][1], subframe[0][0]:subframe[0][1]]
    else:
        dark = master_dark.data

    if exptime_dark != exptime_data:
        if not isbias:
            raise ValueError('dark and data have not the same exptime')
        else:
            print('using bias instead of dark...')
            pass

    # flat
    header_flat = master_flat.header
    filter_flat = header_flat['FILTER']
    if use_subframe:
        flat = master_flat.data[subframe[1][0]:subframe[1][1], subframe[0][0]:subframe[0][1]]
    else:
        flat = master_flat.data

    if filter_flat != filter_data:
        raise ValueError('flat and data have not the same filter')
    
    # calibrate
    data_calib = (data - dark) / (flat)

    return data_calib