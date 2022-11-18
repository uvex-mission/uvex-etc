import astropy.units as u
import numpy as np
import os
from specutils import Spectrum1D

curdir = os.path.dirname(__file__)
datadir = os.path.join(curdir, 'zodi_data')

from . import config
UVEX = config.UVEX()

def diffuse_galactic_fuv_norm(latitude):
    """   
    Generates FUV flux from Table 4 of Murthy et al (2014).
    
    NOTE: Not applicable for observations in the plane
   
    Parameters
    -----------
    latitude : Astropy quantity or array of Astropy quantities
    
    Returns
    --------
    
    fuv_flux : Astropy quantity with PHLAM units
   
    """
    from numpy import sin, abs

    
    south = (latitude<-0*u.deg)
    north = (latitude>=0*u.deg)

    fuv_flux = np.zeros_like(latitude.value)

    fuv_flux[north] = 93.4 + 133.2 / sin(abs(latitude[north]))
    fuv_flux[south] = -205.5 + 401.8 / sin(abs(latitude[south]))

    fuv_flux *= u.ph / (u.cm**2 * u.sr * u.s * u.AA)
    
    return fuv_flux
    
def diffuse_galactic_nuv_norm(latitude):
    """   
    Generates NUV flux from Table 4 of Murthy et al (2014).
    
    NOTE: Not applicable for observations in the plane
   
    Parameters
    -----------
    latitude : Astropy quantity or array of Astropy quantities
    
    Returns
    --------
    
    fuv_flux : Astropy quantity with PHLAM units
   
    """
    from numpy import sin, abs

    
    south = (latitude<-0*u.deg)
    north = (latitude>=0*u.deg)

    fuv_flux = np.zeros_like(latitude.value)

    fuv_flux[north] = 257.5 + 185.1 / sin(abs(latitude[north]))
    fuv_flux[south] = 66.7 + 356.3 / sin(abs(latitude[south]))

    fuv_flux *= u.ph / (u.cm**2 * u.sr * u.s * u.AA)
    
    return fuv_flux
    
def galactic_nuv_spec(latitude, infile=UVEX.ZODI_SPEC):
    """
    Takes as input the galactic latitude and generates a flat (in photon units)
    spectrum with flux scaled from Murthy (2014)

    Parameters
    -----------
    
    latitude : Astropy quantity (degrees)
        Latitude of the observation
    
    Optional Parameters
    -------------------
    infile : str
        File to use to define the input wavelength range. Defaults to Zodi.
    
    Returns
    -------
    Spectrum1D object
    
    """
    
    #### PUT A CHECK HERE TO MAKE SURE OUTISDE OF 15 DEGREES
    
    # Use the Zodi 
    data = np.genfromtxt(os.path.join(datadir, infile))
    
    
    # Get the scale term:
    scale = diffuse_galactic_nuv_norm(latitude)
    

    wave = data[:, 0]*u.AA
    flux = np.zeros_like(wave.value)
    flux += scale.value
    flux *= scale.unit
    
    # Convert to per-pixel units
    flux = flux.to(u.ph /(u.cm**2 * u.Angstrom * u.arcsec**2 * u.s))
    flux *= UVEX.PIXEL

    return Spectrum1D(flux=flux, spectral_axis= wave)

def galactic_fuv_spec(latitude, infile=UVEX.ZODI_SPEC):
    """
    Takes as input the galactic latitude and generates a flat (in photon units)
    spectrum with flux scaled from Murthy (2014)

    Parameters
    -----------
    
    latitude : Astropy quantity (degrees)
        Latitude of the observation
    
    Optional Parameters
    -------------------
    infile : str
        File to use to define the input wavelength range. Defaults to Zodi.
    
    Returns
    -------
    Spectrum1D object
    
    """
    
    #### PUT A CHECK HERE TO MAKE SURE OUTISDE OF 15 DEGREES
    
    # Use the Zodi 
    data = np.genfromtxt(os.path.join(datadir, infile))
    
    
    # Get the scale term:
    scale = diffuse_galactic_fuv_norm(latitude)
    

    wave = data[:, 0]*u.AA
    flux = np.zeros_like(wave.value)
    flux += scale.value
    flux *= scale.unit
    
    # Convert to per-pixel units
    flux = flux.to(u.ph /(u.cm**2 * u.Angstrom * u.arcsec**2 * u.s))
    flux *= UVEX.PIXEL

    return Spectrum1D(flux=flux, spectral_axis= wave)
    
