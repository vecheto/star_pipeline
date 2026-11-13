import os
import signal
import subprocess
import tempfile
import warnings
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning


def header_pointing(header, keys=(('RAJ2000', 'DECJ2000'), ('RA', 'DEC'), ('OBJCTRA', 'OBJCTDEC'))):
    """Approximate pointing from header; accepts degrees or sexagesimal strings. None if absent."""
    for ra_key, dec_key in keys:
        if ra_key in header and dec_key in header:
            ra, dec = header[ra_key], header[dec_key]
            if isinstance(ra, str):
                return SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            return SkyCoord(ra, dec, unit=u.deg)
    return None


def header_pixel_scale(header):
    """arcsec/px from FOCALLEN [mm], XPIXSZ [um, unbinned] and XBINNING. None if absent."""
    try:
        return 206.265 * header['XPIXSZ'] * header.get('XBINNING', 1) / header['FOCALLEN']
    except (KeyError, TypeError, ZeroDivisionError):
        return None


def to_skycoord(target):
    return SkyCoord.from_name(target) if isinstance(target, str) else target


def solve_wcs(image, center, radius=2.0, scale=None, timeout=120):
    """
    Plate-solve an image with a local astrometry.net install (solve-field).

    center  : SkyCoord, approximate field center.
    radius  : search radius around center [deg].
    scale   : approximate pixel scale [arcsec/px]; narrows the search (±15 %).
    timeout : wall-clock limit [s]; solve-field's --cpulimit is not always honored.
    """
    downsample = 2 if min(image.shape) > 1500 else 1
    with tempfile.TemporaryDirectory() as tmp:
        img_path = os.path.join(tmp, 'frame.fits')
        fits.writeto(img_path, np.asarray(image, dtype=np.float32))
        cmd = [
            'solve-field', img_path,
            '--ra', str(center.ra.deg), '--dec', str(center.dec.deg),
            '--radius', str(radius),
            '--downsample', str(downsample),
            '--no-plots', '--overwrite', '--no-verify',
            '--dir', tmp,
        ]
        if scale is not None:
            cmd += ['--scale-units', 'arcsecperpix',
                    '--scale-low', str(0.85 * scale), '--scale-high', str(1.15 * scale)]

        # own process group, so the astrometry-engine child is killed too on timeout
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                text=True, start_new_session=True)
        try:
            out, err = proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.communicate()
            raise RuntimeError(f'solve-field timed out after {timeout} s')

        wcs_path = os.path.join(tmp, 'frame.wcs')
        if proc.returncode != 0 or not os.path.exists(wcs_path):
            raise RuntimeError(f'solve-field failed (exit {proc.returncode}):\n{out[-1500:]}\n{err[-1500:]}')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', FITSFixedWarning)  # .wcs file has no data array
            return WCS(fits.getheader(wcs_path))


def target_pixel(wcs, target):
    """Pixel (x, y), 0-based, of a target given as SkyCoord or a resolvable name."""
    x, y = wcs.world_to_pixel(to_skycoord(target))
    return float(x), float(y)
