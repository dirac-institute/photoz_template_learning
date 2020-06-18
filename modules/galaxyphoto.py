
# Classes and functions to support photometry 
# (i.e. galaxies, filters, and SED's)

import numpy as np


class Galaxy:
    '''
    Galaxy class with photometry and redshift.
    
    source is the name of the data set it came from.
    template is the template it is matched to.
    magToflux/fluxTomag converts between magnitudes and F_lambda, 
    using reference magnitude 25, and reference wavelenth 5000 angstroms
    '''
    def __init__(self, wavelen=None, mags=None, mag_err=None, fluxes=None, 
                flux_err=None, filters=None, redshift=None, source=None,
                template=None, m0=None):
        self.wavelen = wavelen
        self.mags = mags
        self.mag_err = mag_err
        self.fluxes = fluxes
        self.flux_err = flux_err
        self.filters = filters
        self.redshift = redshift
        self.source = source
        self.template = template
        self.m0 = m0

    _mag_ref = 25
    _lambda_ref = 5000 # Angstrom
    def magToflux(self):
        self.fluxes = (self._lambda_ref/self.wavelen)**2 * 10**((self._mag_ref-self.mags)/2.5)
        self.flux_err = self.fluxes/2.5 * np.log(10) * self.mag_err
    def fluxTomag(self):
        self.mags = self._mag_ref - 2.5*np.log10((self.wavelen/self._lambda_ref)**2 * self.fluxes)
        self.mag_err = 2.5/np.log(10) * self.flux_err/self.fluxes


class Bandpass:
    '''
    Class defining a bandpass filter
    
    T is the system respond function
    R is the normalized response function, calculated assuming the
    detector is photon-counting
    '''
    def __init__(self, filename, dlambda=20):
        # load from file
        wavelen,T = np.loadtxt(filename,unpack=True)
        # resample wavelen and calculate R
        self.wavelen = np.arange(min(wavelen),max(wavelen),dlambda)
        self.T = np.interp(self.wavelen,wavelen,T)
        self.R = self.T * self.wavelen
        self.R /= (self.R * dlambda).sum()
        del wavelen,T
        self.mean_wavelen = (self.wavelen * self.R).sum()/self.R.sum()
        self.eff_width = (self.R * dlambda).sum()/max(self.R)
        

def get_bandpass_dict(filter_loc='filters/',dlambda=20):
    """
    Return dictionary of bandpasses using the filter list at filter_loc
    """
    names, files = np.loadtxt(filter_loc+'filters.list', unpack=True, dtype=str)
    bandpass_dict = dict()
    for i,filename in enumerate(files):
        bandpass = Bandpass('filters/'+filename,dlambda)
        bandpass_dict[names[i]] = bandpass  
    return bandpass_dict

def get_mean_wavelen(bandpass_dict, filters=None):
    """
    Return the mean wavelengths for a list of filters
    """
    if filters is None:
        filters = bandpass_dict.keys()
    mean_wavelens = np.array([])
    for name in filters:
        band = bandpass_dict[name]
        mean_wavelens = np.append(mean_wavelens,band.mean_wavelen)
    return mean_wavelens


class Sed:
    """
    An SED defined by F_lambda. Can redshift, or calculate fluxes in bandpasses.

    redshift will redshift the SED to redshift z.
    flux will calculate the flux in the provided bandpass.
    fluxlist will return list of fluxes for the filters in the provided
    bandpass_dict. If you provide the list filters, it will only provide
    fluxes for those filters.
    """
    '''Class defining an SED'''
    def __init__(self, wavelen=None, flambda=None):
        self.wavelen = wavelen
        self.flambda = flambda
        
    def redshift(self, z):
        self.wavelen *= (1 + z)
        
    def flux(self, bandpass):
        y = np.interp(bandpass.wavelen,self.wavelen,self.flambda)
        flux = (y*bandpass.R).sum() * (bandpass.wavelen[1] - bandpass.wavelen[0])
        return flux
    
    def fluxlist(self, bandpass_dict, filters=None):
        if filters is None:
            filters = bandpass_dict.keys()
        fluxes = []
        for name in filters:
            bandpass = bandpass_dict[name]
            fluxes.append(self.flux(bandpass))
        return np.array(fluxes)