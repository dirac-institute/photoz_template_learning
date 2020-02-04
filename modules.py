import numpy as np
import copy
from scipy.interpolate import interp1d
from sklearn.ensemble import IsolationForest



# Classes and function to support the galaxies, filters and Sed's
# --------------------------------------------------------------------------
class Galaxy:
    '''Class defining a galaxy'''
    def __init__(self, wavelen=None, mags=None, mag_err=None, fluxes=None, flux_err=None, 
                                                filters=None, redshift=None, source=None):
        self.wavelen = wavelen
        self.mags = mags
        self.mag_err = mag_err
        self.fluxes = fluxes
        self.flux_err = flux_err
        self.filters = filters
        self.redshift = redshift
        self.source = source

    _mag_ref = 25
    _lambda_ref = 5000 # Angstrom
    def magToflux(self):
        self.fluxes = (self._lambda_ref/self.wavelen)**2 * 10**((self._mag_ref-self.mags)/2.5)
        self.flux_err = self.fluxes/2.5 * np.log(10) * self.mag_err
    def fluxTomag(self):
        self.mags = self._mag_ref - 2.5*np.log10((self.wavelen/self._lambda_ref)**2 * self.fluxes)
        self.mag_err = 2.5/np.log(10) * self.flux_err/self.fluxes
    

class Bandpass:
    '''Class defining a bandpass filter'''
    def __init__(self, filename):
        # load from file
        wavelen,sb = np.loadtxt(filename,unpack=True)
        # resample wavelen and calculate phi
        self.wavelen = np.arange(min(wavelen),max(wavelen),20)
        sb = np.interp(self.wavelen,wavelen,sb)
        self.phi = sb/self.wavelen
        self.phi /= self.phi.sum() * (self.wavelen[1] - self.wavelen[0])
        del wavelen,sb
        # calculate effective wavelen
        self.eff_wavelen = (self.wavelen*self.phi).sum()/self.phi.sum()
        
def get_bandpass_dict():
    names, files = np.loadtxt('filters/filters.list', unpack=True, dtype=str)
    bandpass_dict = dict()
    for i,filename in enumerate(files):
        bandpass = Bandpass('filters/'+filename)
        bandpass_dict[names[i]] = bandpass  
    return bandpass_dict

def get_eff_wavelen(bandpass_dict, filters=None):
    if filters is None:
        filters = bandpass_dict.keys()
    eff_wavelen = []
    for name in filters:
        band = bandpass_dict[name]
        eff_wavelen.append(band.eff_wavelen)
    return np.array(eff_wavelen)
    
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
    
    
    
# Functions for assembling training sets
# --------------------------------------------------------------------------
def match_photometry(template_dict,fluxes,errs,redshift,bandpass_dict):
    
    keys = np.array(list(template_dict.keys()))
    mse_list = np.array([])
    scales = np.array([])
    
    for template in template_dict.values():
        
        sed = copy.deepcopy(template)
        sed.redshift(redshift)
        template_fluxes = sed.fluxlist(bandpass_dict)
        
        scale = np.median(template_fluxes/fluxes)
        template_fluxes_ = template_fluxes/template_fluxes[3]
        fluxes_ = fluxes/fluxes[3]
        errs_ = errs/fluxes
        
        mse = np.mean(1/errs_**2*(template_fluxes_ - fluxes_)**2)
        mse_list = np.append(mse_list,mse)
        scales = np.append(scales,scale)
        
    idx = mse_list.argmin()
    return keys[idx],scales[idx]

def create_training_sets(template_dict,data,bandpass_dict):
    
    # create a dictionary to hold the set for each template
    sets = dict()
    for key in template_dict.keys():
        sets[key] = []
        
    # sort each photometry into the training sets
    for row in data:
        redshift = row[0]
        fluxes = row[1:8]
        errs = row[8:]
        match, scale = match_photometry(template_dict,fluxes,errs,redshift,bandpass_dict)
        
        wavelen = get_eff_wavelen(bandpass_dict)/(1+redshift)
        
        filter_names = list(bandpass_dict.keys())
        for i,wavelen_ in enumerate(wavelen):
            sets[match].append([wavelen_,fluxes[i]*scale,errs[i]*scale,redshift,filter_names[i]])
            
    return sets



# Functions for training templates
# --------------------------------------------------------------------------
def perturb_template(template, data, bandpass_dict, Delta=0.005, mask=0):
    """
    Function that perturbs an SED template in accordance with the
    matched photometry. Definition of terms found in paper I am writing.
    """
    
    nbins = len(template.wavelen)
    wavelen = template.wavelen
    widths = np.array([float(wavelen[1]-wavelen[0]),
                      *[(wavelen[i+1]-wavelen[i-1])/2 for i in range(1,nbins-1)],
                      float(wavelen[-1]-wavelen[-2])])
    
    # create initial M and nu
    M = np.identity(nbins)*1/Delta**2
    nu = np.zeros(nbins)
    
    # if no mask is given, use all data
    if type(mask) == int:
        mask = np.ones(len(data))
    
    # run through all the photometry
    for i,row in enumerate(data):
        
        # skip this row if it's an outlier
        if mask[i] == -1:
            continue
            
        # get observed flux
        obs_flux = row[1]
        sigma = row[2]/obs_flux
        
        # calculate template flux
        redshift = row[3]
        filter_name = row[4]
        template_ = copy.deepcopy(template)
        template_.redshift(redshift)
        template_flux = template_.flux(bandpass_dict[filter_name])
        
        # calculate r^n * Delta lambda
        rn = np.interp(wavelen * (1+redshift), bandpass_dict[filter_name].wavelen, bandpass_dict[filter_name].phi)
        dlambda = widths * (1+redshift)
        rn_dlambda = rn * dlambda
        
        # add to M
        M += 1/sigma**2 * np.outer(rn_dlambda,rn_dlambda)
        
        # add to nu
        nu += 1/sigma**2 * (obs_flux - template_flux) * rn_dlambda
    
        
    # solve the system for the perturbation    
    sol = np.linalg.solve(M,nu)
    return sol


def train_templates(template_dict, data, bandpass_dict, N_rounds=5, N_iter=4, Delta=0.005):
    
    Tdict = copy.deepcopy(template_dict)
    Plist = copy.deepcopy(data)

    for i in range(N_rounds):
        
        print("Round "+str(i+1)+"/"+str(N_rounds))
        
        training_sets = create_training_sets(Tdict,data,bandpass_dict)
        
        for key in Tdict.keys():
            template = Tdict[key]
            training_set = training_sets[key]

            # identify outliers
            x = np.array(training_set)[:,0].astype(float)
            y = np.array(training_set)[:,1].astype(float)
            clf = IsolationForest(max_samples=len(x))
            xy = np.array([x,y]).T
            mask = clf.fit_predict(xy)

            for j in range(N_iter):
                pert = perturb_template(template,training_set,bandpass_dict,Delta=Delta,mask=mask)
                template.flambda += pert
                
    training_sets = create_training_sets(Tdict,Plist,bandpass_dict)
    return Tdict, training_sets
