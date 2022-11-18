import astropy.units as u
from numpy import pi as PI
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord
import os

localpath = os.path.dirname(__file__)

    

class UVEX():
    """
    Container class to hold the CBE values
    """

    def __init__(self, default=False):
    
    

        # Characteristics
        self.EPD = 75*u.cm
        self.FOV = (3.5*u.deg)**2

        # Pixel size:
        self.PIX_UM = 10*u.micron

        # Plate scale for imager and spectrometer
        self.PLATE_SCALE = 1.03*u.arcsec / self.PIX_UM
        self.SPEC_PLATE_SCALE=0.8 * u.arcsec / self.PIX_UM

        # Pixel area
        self.PIXEL=(self.PIX_UM * self.PLATE_SCALE)**2

        # Read noise
        self.READ_NOISE_CBE=2.0 # e-

        # Dark current
        self.DARK_CURRENT_CBE=1e-3 / u.s


        # Number of effective pixels for 2-arcsec FWHM for imagers
        self.NPIX_CBE=10.15

        # Number of effective pixels for spectrometer
        # for 2-arcsec FWHM, estimate is
        self.NPIX_SPEC=3



        # Spectrometer input files
        
        
        # Spectral bins
        self.SPEC_BINS='UVEXS_Spectral_Resolution_R1600.txt'
        
        # Background input files
        
        # Zodiacal Light Spectrum files
        self.ZODI_SPEC='scaled_zodiacal_spec.txt'
        self.ZODI_SPATIAL='Leinert97_table17.txt'
        
#        self.AREA = PI * (self.EPD*0.5)**2


            
    @property
    def AREA(self):
        return 0.85*PI * (self.EPD*0.5)**2

    @property
    def DARK_CURRENT(self):
        return self.DARK_CURRENT_CBE

    @property
    def NPIX(self):
        return self.NPIX_CBE

    @property
    def READ_NOISE(self):
        return self.READ_NOISE_CBE

    @property
    def READ_NOISE(self):
        return self.READ_NOISE_CBE


    
    
## Class for holding Survey information
class Survey():
    """
    Container class to hold the ETC default depth information
    """

    def __init__(self):

        self.set_defaults()
        
    
    def set_defaults(self):
        self.nuv_exposure = 300*u.s
        self.fuv_exposure = 900*u.s

        self.dwell_time = 900*u.s
        self.allsky_dwells = 10
    
    
        # Set the calendar time for this pointing
        self.obstime = Time('2021-02-18 09:00:00')

        # NUV Performacne estiamtes
        # Define location for deep field (ECDFN)
        self.nuv_deep_skycoord = SkyCoord(269.7372, 66.025,
            unit = u.deg, obstime=self.obstime)
        # Use something near the Ecliptic Plane as "average"    
        self.nuv_ave_skycoord = SkyCoord(90., 5, unit = u.deg,
            obstime=self.obstime, frame = 'geocentrictrueecliptic')

        # Define locations for FUV
        # Use ECDF-N again
        self.fuv_deep_skycoord = SkyCoord(269.7372, 66.025,
            unit = u.deg, obstime=self.obstime)
        # Use 15-degrees out of galactic plane
        self.fuv_ave_skycoord = SkyCoord(90., 15, unit = u.deg,
            obstime=self.obstime, frame = 'galactic')

        # Ly-alpha background level in kR
        self.lya_kr = 20 
        
        
    def info(self):
        """
        Report out current values
        """
        v = vars(self)
        for i in v:
            print(i, v[i])
    
    @property
    def allsky_nuv_survey_frames(self):
        '''
        Returns the number of NUV frames
        '''
        return (self.allsky_dwells * self.dwell_time / self.nuv_exposure).value

    @property
    def allsky_fuv_survey_frames(self):
        '''
        Returns the number of NUV frames
        '''
        return (self.allsky_dwells * self.dwell_time / self.fuv_exposure).value

