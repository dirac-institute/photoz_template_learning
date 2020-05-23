
# Functions for assembling template photometry sets from a list of galaxies

import numpy as np
import copy
import multiprocessing as mp
import functools


def create_training_sets(galaxies, template_dict, bandpass_dict, Ncpus=None):
    """
    Returns a dictionary of training sets for each template in template dict.

    Each galaxy is matched and renormalized to one of the templates in 
    template_dict, using fluxes calculated with bandpass_dict.
    Each training set is a list of the matched galaxies.
    The matching is performed in parallel, using the number of CPU's specified
    by Ncpus. If none specified, will use all CPU's available.
    """

    # create the cpu pool to match galaxies in parallel
    Ncpus = mp.cpu_count() if Ncpus is None else Ncpus
    pool = mp.Pool(Ncpus)

    # match the galaxies
    galaxies_ = pool.starmap(match_galaxy,
                                [(galaxy,template_dict,bandpass_dict) 
                                    for galaxy in galaxies])

    # close the cpu pool
    pool.close()

    # assemble the dictionary of matched galaxies
    training_sets = dict()
    for template in template_dict.keys():
        training_sets[template] = []
    for galaxy in galaxies_:
        galaxy.wavelen /= (1 + galaxy.redshift)
        template = galaxy.template
        training_sets[template].append(galaxy)

    return training_sets 



def match_galaxy(galaxy, template_dict, bandpass_dict):
    """
    Return galaxy with galaxy.template equal to the matching template
    in template_dict, using the filters in bandpass_dict.
    """

    galaxy_ = copy.deepcopy(galaxy)
    template,scale = match_photometry(galaxy,template_dict,bandpass_dict)
    galaxy_.fluxes *= scale
    galaxy_.flux_err *= scale
    galaxy_.fluxTomag()
    galaxy_.template = template

    return galaxy_



def match_photometry(galaxy, template_dict, bandpass_dict):
    """
    Return the template and normalization matched to the galaxy photometry.
    """

    # list of templates
    keys = np.array(list(template_dict.keys()))
    # arrays to store the mean square errors and the normalizations
    mse_list = np.array([])
    norms = np.array([])
    
    # calculate mse and norm for each template
    for template in template_dict.values():

        # make a copy of the fluxes and fractional errors
        fluxes = galaxy.fluxes.copy()
        errs = galaxy.flux_err/galaxy.fluxes
        
        # calculate redshift template fluxes
        sed = copy.deepcopy(template)
        sed.redshift(galaxy.redshift)
        template_fluxes = sed.fluxlist(bandpass_dict, galaxy.filters)
        
        # determine the median normalization
        idx = np.where( template_fluxes > 0 )
        norm = np.median(template_fluxes[idx]/fluxes[idx])
        idx = (np.fabs(template_fluxes/fluxes - norm)).argmin()

        # renormalize in median band
        template_fluxes /= template_fluxes[idx]
        fluxes /= fluxes[idx]
        
        # calculate the mse
        mse = np.mean(1/errs**2*(template_fluxes - fluxes)**2)
        mse_list = np.append(mse_list, mse)
        norms = np.append(norms, norm)
        
    # identify the template with the lowest mse, 
    # and return it with the norm
    idx = mse_list.argmin()
    return keys[idx],norms[idx]