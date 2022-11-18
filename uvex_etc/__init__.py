
from astropy.stats import signal_to_noise_oir_ccd
import numpy as np
import astropy.units as u
import os


from . import filters
from . import config
from . import zodi
from . import lyman
from . import galactic
#from . import spectrometer
#from . import requirements

from pkg_resources import get_distribution

__version__ = get_distribution('uvex_etc').version



#UVEX = config.UVEX()

def get_snr(exposure, source, sky, npix=False):
    """Calculate the SNR of an observation of a point source with UVEX.

    Parameters
    ----------
    exposure : astropy.units.Quantity
        The exposure time
    source : astropy.units.Quantity
        Source count rate
    sky : astropy.units.Quantity
        Sky count rate
        
    Optional Parameters
    -------------------
    
    npix : float
        Number of effective pixels. Defaults to UVEX.NPIX

    Returns
    -------
    float
        The signal to noise ratio
    """
    UVEX=config.UVEX()

    if npix is False:
        npix = UVEX.NPIX
    
    return signal_to_noise_oir_ccd(exposure.to(u.s).value,
        source.value,
        sky.value,
        dark_eps=UVEX.DARK_CURRENT.value,
        rd=UVEX.READ_NOISE,
        npix=npix,
        gain=1.    
    )

    
def get_exposure(source, sky, snr = 5., neff=False):
    """
    Wrapper to compute the required exposure to get a SNR of at least 5.
    
    Parameters
    ----------
    source :  astropy.units.Quantity
        Source count rate
    sky : astropy.units.Quantity
        Sky count rate
        
    Optional Parameters
    -------------------
    
    snr : float
        Desired SNR
        
    neff : float
        Number of effective pixels. Defaults to UVEX.NPIX

    Returns
    -------
    float
        The time required to reach the given SNR
    """
    UVEX=config.UVEX()
    if neff is False:
        neff = UVEX.NPIX

    return calc_exposure(snr,
                        source.value, 
                        sky.value + UVEX.DARK_CURRENT.value,
                        UVEX.READ_NOISE,
                        neff)
                        

def get_required_rate(exposure, sky, snr=5, neff=False):
    
    """
    Wrapper to compute the required exposure to get a SNR of at least 5.
    
    Parameters
    ----------
    exposure :  astropy.units.Quantity
        Source count rate
    sky : astropy.units.Quantity
        Sky count rate
        
    Optional Parameters
    -------------------
    
    snr : float
        Desired SNR
        
    neff : float
        Number of effective pixels. Defaults to UVEX.NPIX

    Returns
    -------
    float
        The time required to reach the given SNR
    """
    UVEX=config.UVEX()
    if neff is False:
        neff = UVEX.NPIX
    return _req_source(snr,
                        exposure.to(u.s).value,
                        sky.value + UVEX.DARK_CURRENT.value,
                        UVEX.READ_NOISE,
                        neff)
                        


def calc_exposure(k, src_rate, bgd_rate, read_noise, neff):
    """
    Compute the time to get to a given significance (k) given the source rate,
    the background rate, the read noise, and the number
    of effective background pixels. Inversion of the standard CCD SNR equation.

    -----

    time = calc_exposure(k, src_rate, bgd_rate, read_noise, neff)

    """
    denom = 2 * src_rate**2

    nom1 = (k**2) * (src_rate + neff*bgd_rate)

    nom2 = ( k**4 *(src_rate + neff*bgd_rate)**2 +
                    4 * k**2 * src_rate**2 * neff * read_noise**2)**(0.5)
    exposure = (nom1 + nom2)/ denom
    return exposure
    

def _req_source(k, exposure, bgd_rate, read_noise, neff):
    """
    Isolate source flux to get at least SNR of k in exposure seconds
    
    Parameters
    -----------
    
    k : float
        Desired SNR
    exposure: float
        Exposure in seconds
    bgd_rate : float
        Combined sky and dark current
    read_noise : float
        Read noise per pixel
    neff : float
        Effective number of pixels
    """

    c = neff * k**2 * (read_noise**2 + exposure*(bgd_rate))
    source =  (k**2 + np.sqrt(k**4 + 4*c))/ (2*exposure)
    return source * u.ct / u.s

def limiting_mag(coord, obstime, exposure=300*u.s, 
                band='FUV',diag=False, snr=5, frames = 1, kr=2, **kwargs):
    """
    Returns the limiting magnitude for a given coordinate at a given time.
    
    Assumes a flat AB spectrum.
    
    Parameters
    -----------
    coord : SkyCoord object
        The location you want the limiting magnitude for
    obstime : Astropy Time object
        The time of the observation
        
    Optional Parameters
    -------------------
    band : 'str'
        Which band to use. 'FUV' or 'NUV'. Default is 'FUV'.
    
    exposure : Astropy Units Quantity
        Exposure time to use. Default is 300*u.s
        
    diag : boolean
        Print out diagostics
        
    frames : int
        Number of frames to stack. Assumes that SNR increases as sqrt(frames)
    
    
    """
    UVEX = config.UVEX()

    from synphot import Observation, SourceSpectrum
    from synphot.models import ConstFlux1D
    
    if band == 'FUV':
        bandpass = filters.fuv_bandpass(**kwargs)
        if diag:
            print('get_maglimit: Using FUV bandpass')
    else:
        bandpass = filters.nuv_bandpass()
        if diag:
            print('get_maglimit: Using NUV bandpass')

    # Make sure you convert to the correct coordinates
    radec = coord.icrs
    gal = coord.galactic
    
    # Get the Zodi and backgrounds at this locatio

    zodi_spec = zodi.zodi_spec_coords(radec, obstime)
    galactic_spec = galactic.galactic_nuv_spec(gal.b)
    ly_spec = lyman.lyman_spec(kr=kr)

    band_zodi = Observation(zodi_spec, bandpass, force='extrap')
    band_galactic = Observation(galactic_spec, bandpass, force = 'extrap')
    band_lyman = Observation(ly_spec, bandpass, force='extrap')

    zodi_rate = band_zodi.countrate(area=UVEX.AREA)
    galactic_rate = band_galactic.countrate(area=UVEX.AREA)
    
    if band == 'FUV':
        ly_rate = band_lyman.countrate(area=UVEX.AREA)
        sky = zodi_rate + galactic_rate + ly_rate
    else:
        sky = zodi_rate + galactic_rate

    
    # Get reference count rate
    m_ref = 22*u.ABmag
    sp = SourceSpectrum(ConstFlux1D, amplitude=m_ref)
    obs_band = Observation(sp, bandpass)
    ref_rate = obs_band.countrate(area=UVEX.AREA)
    
    # Get the required rate:
    req_rate = get_required_rate(exposure, sky, snr=snr/np.sqrt(frames))
    ratio = req_rate / ref_rate
    m_limit = m_ref - (2.5*np.log10(ratio))*u.mag
    
    
    if diag:
        print(f'Diagnostics for get_maglimit')
        print(f'---')

        print(f'Sky coordinates: {radec}')
        print(f'Galactic Latitude: {gal.b:8.2f}')
        zlat, zlon = zodi.get_ecliptic_coords(radec, obstime)
        print(f'Ecliptic latitude: {zlat:8.2f}')
        print(f'Ecliptic longitude: {zlon:8.2f}')
        print(f'---')
        tot_bgd = zodi_rate + galactic_rate
        print(f'Zodi rate: {zodi_rate:8.4e}')
        print(f'Galactic Diffuse Rate: {galactic_rate:8.4e}')
        if band == 'FUV':
            print(f'LymanAlpha rate: {ly_rate:8.2e}')
            tot_bgd += ly_rate
        print(f'Total background: {tot_bgd:8.2e}')
        print(f'Required source rate: {req_rate:8.2e}')
        print(f'---')
        print(f'Frame exposure time {exposure:8.2f}')
        print(f'Stacking {frames} frames')
        print(f'Limiting magnitude for flat AB spectrum: {m_limit:8.2f}')
        print(f'---')
        
    return m_limit
    


    


    
    

