import os
import numpy as np
import astropy.units as u
from synphot import Empirical1D, SpectralElement

from . import config
UVEX = config.UVEX()


localpath = os.path.dirname(__file__)
filter_dir = os.path.join(localpath, 'filter_data')


def nuv_bandpass(**kwargs):
    """
    Return the combined NUV bandpass
    """
    data = np.genfromtxt(os.path.join(filter_dir, 'nuv_bandpass.txt'))
    
    wave = data[:,0]*u.AA
    throughput = data[:, 1]
    
    return SpectralElement(Empirical1D,
       points = wave,
       lookup_table = throughput)

def fuv_bandpass(**kwargs):
    """
    Return the combined NUV bandpass
    """
    data = np.genfromtxt(os.path.join(filter_dir, 'fuv_bandpass.txt'))
    
    wave = data[:,0]*u.AA
    throughput = data[:, 1]
    
    return SpectralElement(Empirical1D,
       points = wave,
       lookup_table = throughput)
       
       
def lss_bandpass(**kwargs):
    """
    Return the combined NUV bandpass
    """
    data = np.genfromtxt(os.path.join(filter_dir, 'lss_bandpass.txt'))
    
    wave = data[:,0]*u.AA
    throughput = data[:, 1]
    
    return SpectralElement(Empirical1D,
       points = wave,
       lookup_table = throughput)

def lss_wave(spec_file=UVEX.SPEC_BINS):
    """
    Load the spectral bins for the UVEX-LSS
    
   
    """
    
    import numpy as np
    buffer = 0
    data = np.genfromtxt(os.path.join(filter_dir, spec_file), skip_header=2)
    
    # 0 is wavelength
    # 1 is local position (ignore this)
    # 2 is the local dispersion in nm per pix
    
    wave = data[:, 0]
    dispersion = data[:,2]
    
    min_wave = wave[0]
   
    last_wave = min_wave
    wave_bins = np.zeros(4096-buffer*2)
    
    for b in np.arange(4096-buffer*2):
        bin_disper = np.interp(last_wave, wave, dispersion)    
        bin_wave = last_wave + bin_disper
        wave_bins[b] = bin_wave
        last_wave = bin_wave
        
    return (wave_bins*10)*u.AA 
    
def make_lss_spectrum(binrate, exposure=3600*u.s):
    """
    Produce a "realistic" spectrum including dark current and read noise.

    Parameters
    -----------
    binrate : arr
        Rate of electrons per second in each pixel or array of floats (no units)
    
    nrows : int
        Nummber of rows in the *spatial* direction that the source is spread across
        
    Optional Parameters
    -------------------
    exposure : Astropy unit
        How long you want to expose for. Default is 3600-s
        
    """
    
    dark_counts = UVEX.DARK_CURRENT * exposure
    spec = np.random.poisson(dark_counts + (binrate * exposure).value).astype(float)
    # Add on read noise:
    
    spec += np.random.normal(loc=0, scale=UVEX.READ_NOISE,size=spec.shape)
    spec[spec<0] = 0
    
    return np.floor(spec)

