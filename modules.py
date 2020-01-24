import numpy as np
import copy
from scipy.interpolate import interp1d
from lsst.sims.photUtils import Bandpass, BandpassDict


# load filters
names, files = np.loadtxt('filters/filters.list', unpack=True, dtype=str)
bandpassList = []
for filename in files:
    bandpass = Bandpass()
    bandpass.readThroughput('filters/'+filename)
    bandpassList.append(bandpass)
bandpass_dict = BandpassDict(bandpassList,names)

def get_bandpass_dict():
    return copy.deepcopy(bandpass_dict)

# effective wavelengths of filters
eff_wavelen = np.array([])
for Filter in bandpass_dict.values():
    eff_wavelen = np.append(eff_wavelen,Filter.calcEffWavelen()[0])
    
def get_eff_wavelen():
    return copy.deepcopy(eff_wavelen)

# dictionary to hold functions representing filters for use in perturbation alg.
filter_functions = dict()
for key in bandpass_dict.keys():
    Filter = bandpass_dict[key]
    f = interp1d(Filter.wavelen,Filter.phi,kind='cubic',bounds_error=False,fill_value=0)
    filter_functions[key] = f
    
def get_filter_functions():
    return copy.deepcopy(filter_functions)