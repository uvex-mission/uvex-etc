import astropy.units as u
import numpy as np

import uvex_etc
from . import config

from synphot import SourceSpectrum, Observation, SpectralElement
from synphot.models import ConstFlux1D, BlackBodyNorm1D



def nuv_imaging_performance(survey, deep=False, **kwargs):
    """
    Wrapper script to run the NUV imaging requirements.    
    """
    UVEX = config.UVEX()

    # Imports

    dwell = 10
    if deep:
        radec = survey.nuv_deep_skycoord
    else:
        radec = survey.nuv_ave_skycoord
    m = uvex_etc.limiting_mag(radec.icrs,
                                      survey.obstime,
                                      exposure=survey.nuv_exposure,
                                      frames=survey.allsky_nuv_survey_frames,
                                      band = 'NUV',
                                      kr = survey.lya_kr,
                                      **kwargs)
    return m    


    
def fuv_imaging_performance(survey, deep=False, **kwargs):
    """
    Wrapper script to run the NUV imaging requirements.    
    """
    UVEX = config.UVEX()

    # Imports

    dwell = 10
    if deep:
        radec = survey.fuv_deep_skycoord
    else:
        radec = survey.fuv_ave_skycoord
    m = uvex_etc.limiting_mag(radec.icrs,
                                      survey.obstime,
                                      exposure=survey.nuv_exposure,
                                      frames=survey.allsky_nuv_survey_frames,
                                      band = 'FUV',
                                      kr = survey.lya_kr,
                                      **kwargs)
    return m 
    

# def spec_oiii_performance(spec_req, diag=False, **kwargs):
#     """
#     Performance estimate at OIII, allows redshift
# 
#     """
#     UVEX = config.UVEX()
# 
#     from astropy.constants import c, h
#     import uvex.sensitivity.filters as filters
#     spec_band = filters.spec_bandpass(peak_blaze=spec_req.peak_blaze)
# 
#     wave = spec_req.eff_oiii_line()
#     en = (h * c / wave).to(u.erg)
#     phden = (spec_req.oiii_line_flux  / en)
#     
#     eff_snr = spec_req.oiii_snr / np.sqrt(spec_req.oiii_frames())
# 
#     sky = 0. *u.ct / u.s
#     rate = uvex.sensitivity.get_required_rate(spec_req.lss_exposure,
#         sky, snr=eff_snr, neff=spec_req.oiii_neff)
# 
#     counts = rate * spec_req.lss_exposure
#     req_eff_area = counts / (phden * spec_req.lss_exposure)
# 
# 
#     eff_area = UVEX.AREA * spec_band.__call__(spec_req.eff_oiii_line())
#     margin = (eff_area.value - req_eff_area.value) / req_eff_area.value
#     
# 
#     if diag:
#         print(f'Wavelength {wave}')
#         print(f'Line flux: {spec_req.oiii_line_flux:8.2e}')
#         print(f'Photon fluence rate: {phden:8.2e}')
#         print('')
#         print(f'Exposure time: {spec_req.lss_exposure}')
#         print(f'Allocation: {spec_req.oiii_allocation}')
#         print(f'Frames: {spec_req.oiii_frames()}')
#         print(f'Photon fluence {phden*spec_req.lss_exposure:8.2e}')
#         print()
#         print(f'Required SNR: {spec_req.oiii_snr}')
#         print(f'Effective SNR per frame: {eff_snr:8.2e}')
#         print()
#         print(f'Required rate: {rate:8.2e}')
#         print(f'Required counts {counts:8.2f}')
#         print()
#         print(f'Required effective area {req_eff_area:8.1f}, Actual {eff_area:8.1f}')
#         print(f'Margin {100*margin:8.2f}%')
# 
#     return margin
# 
# 
# def spec_nv_performance(spec_req, diag=False, **kwargs):
#     """
#     Performance estimate at OIII, allows redshift
# 
#     """
#     UVEX = config.UVEX()
# 
#     from astropy.constants import c, h
#     import uvex.sensitivity.filters as filters
#     from uvex.sensitivity.spectrometer import load_spectral_bins
# 
#     spec_bins = load_spectral_bins()
#     spec_band = filters.spec_bandpass(peak_blaze=spec_req.peak_blaze)
#     
# 
#     wave = spec_req.eff_nv_line()
# 
#     in_range = spec_bins[np.abs(spec_bins - wave) < 2*u.AA]
#     dlambda = (in_range[1] - in_range[0])
# 
# 
# 
#     en = (h * c / wave).to(u.erg)
#     phden = (spec_req.nv_cont_flux  / en) * 2*dlambda
#     
#     eff_snr = spec_req.nv_snr / np.sqrt(spec_req.nv_frames())
# 
#     sky = 0. *u.ct / u.s
#     rate = uvex.sensitivity.get_required_rate(spec_req.lss_exposure,
#         sky, snr=eff_snr, neff=spec_req.nv_neff)
# 
#     counts = rate * spec_req.lss_exposure
#     req_eff_area = counts / (phden * spec_req.lss_exposure)
# 
# 
#     eff_area = UVEX.AREA * spec_band.__call__(wave)
#     margin = (eff_area.value - req_eff_area.value) / req_eff_area.value
#     
# 
#     if diag:
#         print(f'Wavelength {wave}')
#         print(f'Local dispersion {dlambda}')
#         print(f'Line flux: {spec_req.nv_cont_flux:8.2e}')
#         print(f'Photon fluence rate: {phden:8.2e}')
#         print('')
#         print(f'Exposure time: {spec_req.lss_exposure}')
#         print(f'Allocation: {spec_req.nv_allocation}')
#         print(f'Frames: {spec_req.nv_frames()}')
#         print(f'Photon fluence {phden*spec_req.lss_exposure:8.2e}')
#         print()
#         print(f'Required SNR: {spec_req.nv_snr}')
#         print(f'Effective SNR per frame: {eff_snr:8.2e}')
#         print()
#         print(f'Required rate: {rate:8.2e}')
#         print(f'Required counts {counts:8.2f}')
#         print()
#         print(f'Required effective area {req_eff_area:8.1f}, Actual {eff_area:8.1f}')
#         print(f'Margin {100*margin:8.2f}%')
# 
#     return margin
# 
