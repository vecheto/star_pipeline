"""TESS Objects of Interest from ExoFOP: coordinates and transit ephemeris.

A target name like 'TOI-7886', 'TOI-7886.01' or 'toi 270.02' is looked up in
https://exofop.ipac.caltech.edu/tess/ ; the answer is cached as JSON (exofop_TOI-<n>.json)
so a reduction can be re-run offline.
"""
import json
import os
import re
import urllib.parse
import urllib.request

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

URL = 'https://exofop.ipac.caltech.edu/tess/download_toi.php?toi={}&output=pipe'
TOI_RE = re.compile(r'^\s*TOI[\s_-]*(\d+)(?:\.(\d+))?\s*$', re.IGNORECASE)

FIELDS = {  # ExoFOP column -> key
    'TIC ID': 'tic', 'TOI': 'toi', 'RA': 'ra', 'Dec': 'dec',
    'PM RA (mas/yr)': 'pmra', 'PM Dec (mas/yr)': 'pmdec',
    'Epoch (BJD)': 'epoch', 'Epoch (BJD) err': 'epoch_err',
    'Period (days)': 'period', 'Period (days) err': 'period_err',
    'Duration (hours)': 'duration', 'Duration (hours) err': 'duration_err',
    'Depth (ppm)': 'depth_ppm', 'TESS Mag': 'tmag',
}


def parse_toi(name):
    """'TOI-7886.01' -> '7886.01', 'TOI-7886' -> '7886'; None if not a TOI name."""
    if not isinstance(name, str):
        return None
    m = TOI_RE.match(name)
    if not m:
        return None
    return m.group(1) + (f'.{int(m.group(2)):02d}' if m.group(2) else '')


def _query(toi, timeout=30):
    req = urllib.request.Request(URL.format(urllib.parse.quote(toi)),
                                 headers={'User-Agent': 'star-pipeline'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        lines = [l for l in r.read().decode('utf-8', 'replace').splitlines() if l.strip()]
    if len(lines) < 2:
        raise LookupError(f"TOI-{toi} not found in ExoFOP")
    header = lines[0].split('|')
    rows = [dict(zip(header, l.split('|'))) for l in lines[1:]]
    out = []
    for row in rows:
        d = {}
        for col, key in FIELDS.items():
            v = row.get(col, '').strip()
            if key in ('tic', 'toi', 'ra', 'dec'):
                d[key] = v
            else:
                try:
                    d[key] = float(v)
                except ValueError:
                    d[key] = None
        out.append(d)
    return out


def toi_info(name, cache_dir=None):
    """Dict with tic, toi, ra, dec (sexagesimal), pm, epoch [BJD_TDB], period [d],
    duration [h], depth_ppm, tmag (+ errors). Without a planet suffix the lowest one is used."""
    toi = parse_toi(name)
    if toi is None:
        raise ValueError(f"not a TOI name: {name!r}")
    cache = os.path.join(cache_dir, f'exofop_TOI-{toi}.json') if cache_dir else None
    if cache and os.path.exists(cache):
        with open(cache) as f:
            return json.load(f)

    rows = sorted(_query(toi), key=lambda d: float(d['toi']))
    if '.' not in toi and len(rows) > 1:
        print(f"TOI-{toi} has {len(rows)} planet candidates ({', '.join(d['toi'] for d in rows)}); "
              f"using TOI-{rows[0]['toi']} (write e.g. 'TOI-{rows[1]['toi']}' for another)")
    info = rows[0]
    print(f"ExoFOP: TOI-{info['toi']} = TIC {info['tic']}, {info['ra']} {info['dec']}, "
          f"P = {info['period']} d, T0 = {info['epoch']} BJD, dur = {info['duration']:.2f} h, "
          f"depth = {info['depth_ppm']:.0f} ppm")
    if cache:
        with open(cache, 'w') as f:
            json.dump(info, f, indent=1)
    return info


def toi_skycoord(info):
    return SkyCoord(info['ra'], info['dec'], unit=(u.hourangle, u.deg))


def jd_utc_to_bjd_tdb(jd_utc, coord):
    """BJD_TDB of JD_UTC times for a target (geocentric observer: error < 0.03 s)."""
    t = Time(jd_utc, format='jd', scale='utc',
             location=EarthLocation.from_geocentric(0, 0, 0, unit=u.m))
    return (t.tdb + t.light_travel_time(coord, 'barycentric')).jd


def transits_in(info, bjd_start, bjd_end):
    """Predicted transits overlapping [bjd_start, bjd_end]:
    list of dicts tc, t1 (ingress), t4 (egress), tc_err [BJD_TDB, days]."""
    T0, P, dur = info['epoch'], info['period'], info['duration']
    if None in (T0, P, dur):
        return []
    half = dur / 24 / 2
    n0 = int(np.floor((bjd_start - half - T0) / P))
    n1 = int(np.ceil((bjd_end + half - T0) / P))
    out = []
    for n in range(n0, n1 + 1):
        tc = T0 + n * P
        if tc + half < bjd_start or tc - half > bjd_end:
            continue
        err = np.hypot(info.get('epoch_err') or 0, n * (info.get('period_err') or 0))
        out.append({'epoch': n, 'tc': tc, 't1': tc - half, 't4': tc + half, 'tc_err': err})
    return out
