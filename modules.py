import numpy as np
from scipy.interpolate import interp1d


class Bandpass:
    '''Class defining a bandpass filter'''
    def __init__(self, filename):
        # load from file
        wavelen,sb = np.loadtxt(filename,unpack=True)
        # resample wavelen and calculate phi
        self.wavelen = np.arange(min(wavelen),max(wavelen),1)
        sb = np.interp(self.wavelen,wavelen,sb)
        self.phi = sb/self.wavelen
        self.phi /= self.phi.sum() * (self.wavelen[1] - self.wavelen[0])
        del wavelen,sb
        # calculate effective wavelen
        self.eff_wavelen = (self.wavelen*self.phi).sum()/self.phi.sum()
        
def get_bandpass_dict(filter_loc):
    names, files = np.loadtxt(filter_loc+'filters.list', unpack=True, dtype=str)
    bandpass_dict = dict()
    for i,filename in enumerate(files):
        bandpass = Bandpass('filters/'+filename)
        bandpass_dict[names[i]] = bandpass  
    return bandpass_dict

def get_eff_wavelen(bandpass_dict):
    eff_wavelen = []
    for band in bandpass_dict.values():
        eff_wavelen.append(band.eff_wavelen)
    return np.array(eff_wavelen)
    
def get_bandpass_functions(bandpass_dict):
    bandpass_functions = dict()
    for key in bandpass_dict.keys():
        bandpass = bandpass_dict[key]
        f = interp1d(bandpass.wavelen,bandpass.phi,bounds_error=False,fill_value=0)
        bandpass_functions[key] = f
    return bandpass_functions
    
class Sed:
    '''Class defining an SED'''
    def __init__(self, wavelen=None, flambda=None):
        self.wavelen = wavelen
        self.flambda = flambda
        
    def redshift(self,z):
        self.wavelen *= (1 + z)
        
    def flux(self,bandpass):
        y = np.interp(bandpass.wavelen,self.wavelen,self.flambda)
        flux = (y*bandpass.phi).sum() * (bandpass.wavelen[1] - bandpass.wavelen[0])
        return flux
    
    def fluxlist(self,bandpass_dict):
        fluxes = []
        for bandpass in bandpass_dict.values():
            fluxes.append(self.flux(bandpass))
        return np.array(fluxes)