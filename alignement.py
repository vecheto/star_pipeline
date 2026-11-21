from skimage.registration import phase_cross_correlation
from scipy.ndimage import fourier_shift
import numpy as np
import matplotlib.pyplot as plt

def phase_correlation_alignment(images, upsamplefactor=10):
    # La primera imagen es la referencia
    ref_image = images[0]
    aligned = [ref_image]
    
    for img in images[1:]:
        # Calcular desplazamiento
        shift, error, diffphase = phase_cross_correlation(np.array(ref_image), img, upsample_factor=upsamplefactor,
                                                          normalization=None)
        
        # Aplicar corrección de fase (Fourier)
        offset_image = fourier_shift(np.fft.fftn(img), shift)
        offset_image = np.fft.ifftn(offset_image).real
        
        aligned.append(offset_image)
    
    return aligned

def _registration_image(img):
    """Image used only to measure the registration: hot pixels removed (3x3 median),
    background removed (31x31 median), bright cores clipped so stars, not the sensor's
    fixed pattern, dominate the correlation."""
    from scipy.ndimage import median_filter
    img = median_filter(np.asarray(img, dtype=float), 3)
    img = img - median_filter(img, 31)
    return np.clip(img, 0, np.percentile(img, 99.9))


def _rotation_match(ref, img, angles, upsamplefactor):
    """Best (corr, angle, shift) of img rotated by each angle and shifted onto ref."""
    from scipy.ndimage import rotate, shift as nd_shift
    best = (-np.inf, None, None)
    for angle in angles:
        rot = rotate(img, angle, reshape=False, order=1)
        s, _, _ = phase_cross_correlation(ref, rot, upsample_factor=upsamplefactor, normalization=None)
        corr = np.corrcoef(ref.ravel(), nd_shift(rot, s, order=1).ravel())[0, 1]
        if corr > best[0]:
            best = (corr, angle, s)
    return best


def rotation_alignment(images, upsamplefactor=10, local_range=2.0, local_step=0.1,
                       coarse_step=2.0, min_corr=0.5, verbose=True):
    """
    Align on the first image correcting rotation + translation (meridian flips, field
    rotation near the pole). The angle is searched around the previous frame's angle
    (+-local_range deg); if the correlation is below min_corr, a full 0-360 deg search is done
    and refined. Uncovered pixels are NaN. Returns (aligned, angles, shifts, corrs).
    """
    from scipy.ndimage import rotate, shift as nd_shift
    ref = _registration_image(images[0])
    aligned, angles, shifts, corrs = [], [], [], []
    prev = 0.0
    for i, img in enumerate(images):
        reg = _registration_image(img)
        local = np.arange(prev - local_range, prev + local_range + 1e-9, local_step)
        corr, angle, s = _rotation_match(ref, reg, local, upsamplefactor)
        if corr < min_corr:
            _, coarse, _ = _rotation_match(ref, reg, np.arange(0, 360, coarse_step), 1)
            fine = np.arange(coarse - coarse_step, coarse + coarse_step + 1e-9, local_step)
            corr, angle, s = _rotation_match(ref, reg, fine, upsamplefactor)
        angle = float(angle % 360)
        prev = angle
        if verbose:
            flag = '' if corr >= min_corr else '  <-- WARNING: poor match'
            print(f"  frame {i}: rotation {angle:7.2f} deg, shift ({s[1]:6.1f}, {s[0]:6.1f}) px, "
                  f"corr {corr:.2f}{flag}")
        img = np.asarray(img, dtype=float)
        rot = rotate(img, angle, reshape=False, order=1, cval=np.nan)
        aligned.append(nd_shift(rot, s, order=1, cval=np.nan))
        angles.append(angle)
        shifts.append(s)
        corrs.append(corr)
    return aligned, np.array(angles), np.array(shifts), np.array(corrs)


def stack(images, method):
    if method=='mean':
        return np.mean(images, axis=0)
    if method=='median':
        return np.median(images, axis=0)
    if method=='max':
        return np.max(images, axis=0)
    if method=='sum':
        return np.sum(images, axis=0)