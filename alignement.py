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

def stack(images, method):
    if method=='mean':
        return np.mean(images, axis=0)
    if method=='median':
        return np.median(images, axis=0)
    if method=='max':
        return np.max(images, axis=0)
    if method=='sum':
        return np.sum(images, axis=0)