#from synphot import SourceSpectrum, units, Observation, Empirical1D, SpectralElement
from specutils import Spectrum1D
import os
import numpy as np
import astropy.units as u


curdir = os.path.dirname(__file__)
datadir = os.path.join(curdir, 'zodi_data')

from . import filters
from . import config
UVEX = config.UVEX()

def wavelength_to_energy(wave):
    """
    Convert wavelength to energy. Wavelength must have unit.

    Examples
    --------
    >>> from astropy import units as u
    >>> energy = wavelength_to_energy(1239.84193 * u.nm)
    >>> np.isclose(energy.to(u.eV).value, 1)
    True
    """
    from astropy import units as u
    from astropy import constants as c
    return (c.h * c.c / wave).to(u.Joule)



def lyman_spec(kr = 30, infile=UVEX.ZODI_SPEC):
    """
    Generate Lyman-Alpha emission at 1216 Angstroms.
    
    Default is for 1 kR
    R == 1e6 / (4pi) ph / cm2 / s / sr
    == 3.15 × 10−17 erg cm-2 s-1 arcsec-2 per COS handbook

    
    Optional Parameters
    -------------------    
    kr : float
        See. Default is kr=30, or 30e3 R of LyAlpha

    Returns
    -------
    Spectrum1D
    
    """
    
    # Load Zodi spec just to have the same wavelength coverage
    data = np.genfromtxt(os.path.join(datadir, infile))
    wave = data[:,0] *u.AA
    # Scale term:
    
    mean = 1216 * u.AA
    
    ph_ergs = wavelength_to_energy(mean).to(u.erg)
    
    
    fluence = kr*1e3 * 3.15e-17 * UVEX.PIXEL.value / ph_ergs.value # ph / cm2 / sec
    dl = wave[1] - wave[0] # delta lambda, in 
    phflux = fluence / dl.to(u.AA).value # ph / cm2 / sec / A

    flux = np.zeros_like(wave.value)
    flux[np.argmin(np.abs(wave-mean))]+= phflux
    # flux[[0]] += flux
    flux *= u.ph /(u.cm**2 * u.Angstrom * u.s)

    return Spectrum1D(flux=flux, spectral_axis= wave)
    
    
def lyman_spec_plot(spec):
    '''
    Plots the LymanAlpha light spectrum contribution per pixel.
    
    Parameters
    -------
    spec
        Spectrum1D : LymanAlpha light spectrum produced by zodi_spec()
    
    Returns
    -------
    Plot of the zodiacal light spectrum per pixel.
    '''
    from matplotlib import pyplot as plt
    ax = plt.subplots()[1]  
    ax.plot(spec.spectral_axis, spec.flux)  
    ax.set_xlabel(spec.spectral_axis.unit)  
    ax.set_ylabel(spec.flux.unit)
    ax.set_title('LymanAlpha light per pixel')
    ax.set_xlim(1000, 5e3)
    plt.show()
    return