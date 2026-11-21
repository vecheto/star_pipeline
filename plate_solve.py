import os
import shutil
import signal
import subprocess
import tempfile
import warnings
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning

import exofop


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


def to_skycoord(target, cache_dir=None, use_exofop=False):
    """SkyCoord of a target: with use_exofop, 'TOI-…' names come from ExoFOP (cached in
    cache_dir); otherwise names go to Sesame. A SkyCoord is returned as is."""
    if not isinstance(target, str):
        return target
    if use_exofop and exofop.parse_toi(target):
        return exofop.toi_skycoord(exofop.toi_info(target, cache_dir))
    return SkyCoord.from_name(target)


ASTAP_CANDIDATES = (
    'astap_cli', 'astap',
    r'C:\Program Files\astap\astap_cli.exe',
    '/opt/astap/astap_cli', '/opt/astap/astap',
    '/Applications/ASTAP.app/Contents/MacOS/astap',
)


def find_astap():
    """Path to the ASTAP command-line solver, or None."""
    for cand in ASTAP_CANDIDATES:
        path = shutil.which(cand) or (cand if os.path.isfile(cand) else None)
        if path:
            return path
    return None


def _run(cmd, timeout, name):
    """Run cmd in its own process group; kill the whole group on timeout. Returns (rc, out, err)."""
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, start_new_session=True)
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        if hasattr(os, 'killpg'):
            os.killpg(proc.pid, signal.SIGKILL)
        else:
            proc.kill()
        proc.communicate()
        raise RuntimeError(f'{name} timed out after {timeout} s')
    return proc.returncode, out, err


def _read_wcs(path):
    header = fits.getheader(path)
    header['RADESYS'] = 'ICRS'   # ASTAP copies the input RADESYS (e.g. GAPPT); solutions are ICRS
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FITSFixedWarning)  # .wcs file has no data array
        return WCS(header)


def solve_wcs(image, center, radius=2.0, scale=None, timeout=120, solver='auto', downsample=None):
    """
    Plate-solve an image with a local solver: ASTAP (astap_cli, needs a star database such
    as D50) or astrometry.net (solve-field).

    center  : SkyCoord, approximate field center.
    radius  : search radius around center [deg].
    scale   : approximate pixel scale [arcsec/px]; narrows the search.
    timeout : wall-clock limit [s] (solve-field's --cpulimit is not always honored).
    solver  : 'auto' (ASTAP if installed, else solve-field), 'astap' or 'astrometry'.
    downsample : int >= 1 binning applied by the solver before detecting stars (2-4 helps
                 with defocused/donut stars); None = solver default. The WCS is always
                 returned in full-resolution pixels.
    """
    if downsample is not None and (int(downsample) != downsample or downsample < 1):
        raise ValueError(f'downsample must be None or an integer >= 1, got {downsample}')
    if solver not in ('auto', 'astap', 'astrometry'):
        raise ValueError(f"solver must be 'auto', 'astap' or 'astrometry', got '{solver}'")
    astap = find_astap() if solver in ('auto', 'astap') else None
    if astap:
        return _solve_astap(astap, image, center, radius, scale, timeout, downsample)
    if solver == 'astap':
        raise RuntimeError('ASTAP not found (install it from hnsky.org, with a star database such as D50)')
    if shutil.which('solve-field') is None:
        raise RuntimeError('no plate solver found: install ASTAP (Windows/Linux/macOS) or astrometry.net')
    return _solve_astrometry(image, center, radius, scale, timeout, downsample)


def _solve_astap(astap, image, center, radius, scale, timeout, downsample=None):
    fov = 0 if scale is None else scale * image.shape[0] / 3600   # image height [deg]; 0 = auto
    with tempfile.TemporaryDirectory() as tmp:
        img_path = os.path.join(tmp, 'frame.fits')
        fits.writeto(img_path, np.asarray(image, dtype=np.float32))
        cmd = [astap, '-f', img_path,
               '-ra', str(center.ra.hour), '-spd', str(center.dec.deg + 90),
               '-r', str(radius), '-fov', str(fov), '-z', str(downsample or 0), '-wcs']   # -z 0 = auto
        rc, out, err = _run(cmd, timeout, 'ASTAP')
        wcs_path = os.path.join(tmp, 'frame.wcs')
        if rc != 0 or not os.path.exists(wcs_path):
            raise RuntimeError(f'ASTAP failed (exit {rc}):\n{out[-1500:]}\n{err[-1500:]}')
        return _read_wcs(wcs_path)


def _solve_astrometry(image, center, radius, scale, timeout, downsample=None):
    if downsample is None:
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
        rc, out, err = _run(cmd, timeout, 'solve-field')
        wcs_path = os.path.join(tmp, 'frame.wcs')
        if rc != 0 or not os.path.exists(wcs_path):
            raise RuntimeError(f'solve-field failed (exit {rc}):\n{out[-1500:]}\n{err[-1500:]}')
        return _read_wcs(wcs_path)


def target_pixel(wcs, target):
    """Pixel (x, y), 0-based, of a target given as SkyCoord or a resolvable name."""
    x, y = wcs.world_to_pixel(to_skycoord(target))
    return float(x), float(y)


def locate_target(image, header, catalog, target_name, target_xy=None, search_radius=2.0,
                  timeout=120, solver='auto', cache_dir=None, use_exofop=False, downsample=None):
    """
    Target pixel on `image`: manual target_xy, else plate solve centered on the target,
    then on the header pointing, else the catalog source closest to the image center.

    catalog : DataFrame with 'x', 'y' columns (detections on image), for the fallback.
    Returns ([x, y], method).
    """
    if target_xy is not None:
        print("using target_xy from params")
        return list(map(float, target_xy)), 'manual'

    ny, nx = image.shape
    target = to_skycoord(target_name, cache_dir=cache_dir, use_exofop=use_exofop)
    scale = header_pixel_scale(header)
    centers = [('target coordinates', target),
               ('header pointing', header_pointing(header))]
    for label, center in centers:
        if center is None:
            continue
        print(f"plate solving around {label} ({center.to_string('hmsdms')})")
        try:
            wcs = solve_wcs(image, center, radius=search_radius, scale=scale,
                            timeout=timeout, solver=solver, downsample=downsample)
        except RuntimeError as e:
            print(f"  failed: {str(e).splitlines()[0]}")
            continue
        x, y = target_pixel(wcs, target)
        if 0 <= x < nx and 0 <= y < ny:
            return [x, y], f'plate solve ({label})'
        print(f"  solved, but target falls outside the image ({x:.0f}, {y:.0f})")

    print("WARNING: astrometry failed; assuming the target is the source closest to the image "
          "center. Check finder.png, or set target_xy in the params file.")
    i = np.argmin(np.hypot(catalog['x'] - nx / 2, catalog['y'] - ny / 2))
    return [float(catalog['x'].iloc[i]), float(catalog['y'].iloc[i])], 'closest to image center'
