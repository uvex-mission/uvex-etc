from synphot import SourceSpectrum, units, Observation, Empirical1D, SpectralElement
import astropy.units as u
import numpy as np
from specutils import Spectrum1D
import matplotlib.pyplot as plt
from importlib import reload
import os

curdir = os.path.dirname(__file__)
datadir = os.path.join(curdir, 'zodi_data')

from . import filters
from . import config

UVEX = config.UVEX()

def zodi_spec(scale = 77, infile=UVEX.ZODI_SPEC):
    """
    From here
    https://cads.iiap.res.in/tools/zodiacalCalc/Documentation

    They provied the tabulated Zodi spectrum:
    https://cads.iiap.res.in/download/simulations/scaled_zodiacal_spec.txt

    From here
    Colina et al
    http://adsabs.harvard.edu/abs/1996AJ....112..307C

    According to this, the flux has been scaled so that at 5000 Ang
    the flux has units of 252 W / m2 / sr / micon

    This takes as input the flux density as read off from Table 17 of
    https://aas.aanda.org/articles/aas/pdf/1998/01/ds1449.pdf

    scale = scale in units of [1e-8 W / m2 / sr / micron at 500 nm]
    
    Default is for polar zodiacal emission, which is 77 in the above units.

    Toward the ecliptic plane this number can grow to be >1000

    For a Sun avoidance of 45 degrees this looks like a value of 200 - 900
    based strongly on the heliocentric longitdue. However, if you try
    72, 300, and 1000 it looks like you'll probably span this space.
    
    Optional Parameters
    -------------------
    scale : float
        See above for definition. Default is 77 (suitabled for NEP)
    

    Returns
    -------
    Spectrum1D
    
    """
    
    data = np.genfromtxt(os.path.join(datadir, infile))
    # Scale term:
    scale = scale * u.ph /(u.cm**2 * u.Angstrom * u.sr * u.s)

    wave = data[:, 0]*u.AA
    flux = data[:, 1] * scale

    # Convert to per-pixel units
    flux = flux.to(u.ph /(u.cm**2 * u.Angstrom * u.arcsec**2 * u.s))
    flux *= UVEX.PIXEL

    return Spectrum1D(flux=flux, spectral_axis= wave)

def zodi_spec_coords(coord, time, diag=False):
    """
    Wrapper to return the Zodi appropriately scaled for the given coordinates and time
    """

    # Convert to ecliptic coords:
    lat, lon = get_ecliptic_coords(coord, time)
    if lat > 75:
        scale = 72
    else:
        model = load_zodi_spatial()
        scale = model(lat, lon)

    if diag:
        print(f'zodi_spec_coords:')
        print(f'Ecliptic lat/lon: {lat:8.2f} / {lon:8.2f}')
        print(f'Zodi scale: {scale}')
    return zodi_spec(scale=scale)

def zodi_spec_plot(spec):
    '''
    Plots the zodiacal light spectrum contribution per pixel.
    
    Parameters
    -------
    spec
        Spectrum1D : Zodiacal light spectrum produced by zodi_spec()
    
    Returns
    -------
    Plot of the zodiacal light spectrum per pixel.
    '''
    ax = plt.subplots()[1]  
    ax.plot(spec.spectral_axis, spec.flux)  
    ax.set_xlabel(spec.spectral_axis.unit)  
    ax.set_ylabel(spec.flux.unit)
    ax.set_title('Zodiacal light per pixel')
    ax.set_xlim(1000, 5e3)
    plt.show()
    return

def zodi_nuv(spec):
    '''
    Calculates the UVEX NUV-bandpassed zodiacal light spectrum
    
    Parameters
    -------
    spec
        Spectrum1D : Zodiacal light spectrum produced by zodi_spec()
    
    Returns
    -------
    Spectrum1D
    '''
    nuv_band = filters.nuv_bandpass()
    zodi_nuv = Observation(spec, nuv_band, force='extrap')
    return zodi_nuv

def zodi_fuv(spec):
    '''
    Calculates the UVEX FUV-bandpassed zodiacal light spectrum
    
    Parameters
    -------
    spec
        Spectrum1D : Zodiacal light spectrum produced by zodi_spec()
    
    Returns
    -------
    Spectrum1D
    '''
    fuv_band = filters.fuv_bandpass()
    zodi_fuv = Observation(spec, fuv_band, force='extrap')
    return zodi_fuv

def zodi_bandpass_plot(spec,filt='both'):
    '''
    Plots the zodiacal light spectrum convolved with the UVEX FUV and/or NUV bandpasses.
    
    Parameters
    -------
    spec
        Spectrum1D : Zodiacal light spectrum produced by zodi_spec()
    filt
        str : 'nuv', 'fuv', or 'both'
    
    Returns
    -------
    Plot of the bandpassed zodiacal light spectrum.
    '''
    if filt == 'both':
        spec_nuv = zodi_nuv(spec)
        spec_fuv = zodi_fuv(spec)
        
        fig, axes = plt.subplots(1,2,figsize=(12,4))
        axes[0].plot(spec.spectral_axis, spec_nuv(spec.spectral_axis))
        axes[0].set_xlim(1100, 5e3)
        axes[0].set_xlabel(spec.spectral_axis.unit)
        axes[0].set_ylabel(spec.flux.unit)
        axes[0].set_title('Zodi in NUV')
        
        axes[1].plot(spec.spectral_axis, spec_fuv(spec.spectral_axis))
        axes[1].set_xlim(1000, 5e3)
        axes[1].set_xlabel(spec.spectral_axis.unit)
        axes[1].set_ylabel(spec.flux.unit)
        axes[1].set_yscale('Log')
        axes[1].set_title('Zodi in FUV')
        
        
        
        pix_rate_nuv = spec_nuv.countrate(area=UVEX.AREA)
        pix_rate_fuv = spec_fuv.countrate(area=UVEX.AREA)
        print(f'Estimated counts per pixel (NUV): {pix_rate_nuv:8.2e}')
        print(f'Estimated counts per pixel (FUV): {pix_rate_fuv:8.2e}')
        
        plt.show()
        
    elif filt == 'nuv':
        spec_nuv = zodi_nuv(spec)
        
        ax = plt.figure().subplots()
        ax.plot(spec.spectral_axis, spec_nuv(spec.spectral_axis))
        ax.set_xlim(1100, 5e3)
        ax.set_xlabel(spec.spectral_axis.unit)
        ax.set_ylabel(spec.flux.unit)
        ax.set_title('Zodi in NUV')
        
        pix_rate_nuv = spec_nuv.countrate(area=UVEX.AREA)
        print(f'Estimated counts per pixel (NUV): {pix_rate_nuv:8.2e}')
        
        plt.show()
    
    elif filt == 'fuv':
        spec_fuv = zodi_fuv(spec)
        
        ax = plt.figure().subplots()
        ax.plot(spec.spectral_axis, spec_fuv(spec.spectral_axis))
        ax.set_xlim(1000, 5e3)
        ax.set_xlabel(spec.spectral_axis.unit)
        ax.set_ylabel(spec.flux.unit)
        ax.set_title('Zodi in FUV')
        
        pix_rate_fuv = spec_fuv.countrate(area=UVEX.AREA)
        print(f'Estimated counts per pixel (FUV): {pix_rate_fuv:8.2e}')
        
        plt.show()
    else:
        raise Exception('Filter must be nuv, fuv, or both.')
    return
    
    
def load_zodi_spatial(infile=UVEX.ZODI_SPATIAL):
    """
    Loads the 2D spatial map of the Zodi
    
    Data from Leinert 1997, table 17
    Units are compatible with the zodi_spec.scale parameter
    
    """
    from scipy.interpolate import interp2d
    data = np.genfromtxt(os.path.join(datadir, infile))
    lon = data[1:, 0]
    lat = data[0, 1:]
    zodi = data[1:, 1:]
    return interp2d(lat, lon, zodi)

def get_zodi_norm(coord, time, diag=False):
    """
    Compute the ecliptic positions relative to the sun at the given time
    
    Parameters
    ----------
    coord - SkyCoord object
        Target location
    
    time - Astropy Time object
        Observation epoch
    
        
    """
    from astropy.coordinates import SkyCoord, get_sun, GeocentricTrueEcliptic
    
    sun = get_sun(time).transform_to(GeocentricTrueEcliptic(equinox=time))
    target = coord.transform_to(GeocentricTrueEcliptic(equinox=time))
    
    # Latitude is already ecliptic lat:
    lat = np.abs(target.lat.deg)
    
    # If lat > 70 degrees just leave this at 72
    if lat > 75:
        return 72
    
    # Compute longitude wrt the sun:
    lon = np.abs((target.lon - sun.lon).wrap_at(180 * u.deg).deg)
    
    if diag:
        print(coord)
        print(f'{lat} {lon}')
    model = load_zodi_spatial()
    return model(lat, lon)



def get_ecliptic_coords(coord, time):
    """
    Compute the ecliptic positions relative to the sun at the given time
    
    Parameters
    ----------
    coord : SkyCoord object
        Target location
    
    time : Astropy Time object
        Observation epoch
    
    Retuns
    -------
    lat, lon : float
    
    """
    from astropy.coordinates import SkyCoord, get_sun, GeocentricTrueEcliptic
    
    sun = get_sun(time).transform_to(GeocentricTrueEcliptic(equinox=time))
    target = coord.transform_to(GeocentricTrueEcliptic(equinox=time))
    
    # Latitude is already ecliptic lat:
    lat = np.abs(target.lat.deg)
    # Compute longitude wrt the sun:
    lon = np.abs((target.lon - sun.lon).wrap_at(180 * u.deg).deg)
    
    return lat, lon
    
    