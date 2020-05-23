
import numpy as np 
import copy
import multiprocessing as mp
from scipy.signal import medfilt
from modules.galaxyphoto import Sed
from modules.photomatching import create_training_sets


def log_norm(x, mode, sigma, norm):
    """
    Log normal distribution normalized at the wavelength "norm"
    """
    mu = np.log(mode) + sigma**2
    f = lambda x: 1/(x*sigma*np.sqrt(2*np.pi))*np.exp(-(np.log(x)-mu)**2/(2*sigma**2))
    return f(x)/f(norm)


def new_naive_templates(N, res=100, x_min=10, x_max=15000,
                                    mode_min=1000, mode_max=5500,
                                    sigma_min=0.35, sigma_max=0.9,
                                    norm=5000):
    """
    Returns a dictionary of new naive templates.

    The naive templates are log-normal distributions, with mode and sigma linearly sampled
    from the designated ranges, and normalized at the wavelength "norm"
    """
    
    modes = np.linspace(mode_max,mode_min,N)
    sigmas = np.linspace(sigma_min,sigma_max,N)

    template_dict = dict()
    x = np.arange(x_min, x_max, res, dtype=float)
    for i in range(N):
        template = Sed()
        template.wavelen = x
        template.flambda = log_norm(x, modes[i], sigmas[i], norm)
        template_dict["N"+str(N)+"-"+str(i+1)] = template
        
    return template_dict



def train_templates(galaxies, template_dict, bandpass_dict,
                    w=0.5, Delta=None, dmse_stop=0.05, 
                    maxRounds=None, maxPerts=None, renorm=5000, 
                    Ncpus=None, verbose=True):
    """
    Trains a dictionary of templates on a list of galaxies, using the algorithm
    described in the paper. Returns the trained templates, and a training history
    which is a dictionary containing the SED and wMSE for every step of the training.

    galaxies is a list of galaxy objects.
    template_dict is a dictionary of SED templates.
    bandpass_dict is a dictionary of the filters used to observe the photometry.
    See galaxyphoto.py for details on all of these objects.

    w is the training ratio described in the paper, which will be used to calculate Delta.
    You can manually set Delta instead.
    dmse_stop is the threshold for fractional change in wMSE that ends perturbations.
    You can set maxRounds and maxPerts to set a maximum number of rounds/perturbations per
    round. These are just maxima. If you want to enforce those numbers of rounds/perturbations,
    then set dmse_stop=0.
    Ncpus is the number of cpus to use in the photometry matching/template training when parallelizing
    """

    if verbose:
        print("Columns: Template, number of perturbations, initial/final wMSE")

    # new templates to be perturbed
    new_templates = copy.deepcopy(template_dict)

    # create the history and mse0 dictionaries
    history = dict()
    mse0 = dict()
    for key in new_templates:
        history[key] = dict()
        mse0[key] = 1e9

    roundN = 0
    while True:

        roundN += 1
        print("Round",roundN)

        # create the training sets for this round
        training_sets = create_training_sets(galaxies, new_templates,
                                             bandpass_dict, Ncpus)

        # create the cpu pool to perturb galaxies in parallel
        Ncpus = mp.cpu_count() if Ncpus is None else Ncpus
        with mp.Pool(Ncpus) as pool:
            # perturb each template (unless already matches photometry well)
            roundResult = pool.starmap(perturbation_round,
                                        [(key, training_sets[key],
                                        new_templates[key],
                                        bandpass_dict, mse0[key], 
                                        w, Delta, dmse_stop, maxPerts)
                                        for key in new_templates.keys()])

        # Process the results of the round
        totPerts = 0 
        for i in roundResult:
            key = i[0]
            templates = i[1]
            mse_list = i[2]

            new_templates[key] = templates[-1]
            mse0[key] = mse_list[-1]
            history[key][roundN] = [templates,mse_list]

            totPerts += len(templates)-1

            if verbose:
                print("{:<8}".format(key),
                      "{:>2}".format(len(templates)-1),
                      "{:>10.1f}".format(mse_list[0]),
                      "{:>10.1f}".format(mse_list[-1]))

        # End the training if no templates were perturbed or if we have reached
        # the max number of rounds
        if totPerts == 0 or roundN == maxRounds:
            break

    if renorm != False:
        print("Renormalizing templates at",renorm,"angstroms")
        for template in new_templates.values():
            # regrid on 100 angstrom grid and median filter
            X = np.arange(template.wavelen[0], template.wavelen[-1], 100)
            Y = np.interp(X, template.wavelen, template.flambda)
            Y = medfilt(Y, kernel_size=9)
            # divide by magnitude at the designated wavelength
            idx = np.fabs(template.wavelen - renorm).argmin()
            template.flambda /= Y[idx]

    print("Done!")
    return new_templates, history



def perturbation_round(key,training_set, template, bandpass_dict, mse0,
                        w=0.5, Delta=None, dmse_stop=0.05, maxPerts=None):
    """
    Perform a round of perturbations on a specific template, stopping according
    to dmse_stop or when maxPerts is reached. See the definition of train_templates()
    for more details.

    mse0 is the wMSE from the previous round. It will be used to determined if the 
    new photometry set is different enough to warrant further perturbations. It is
    then updated throughout this function to track when perturbations should end.
    """
    
    # fractional change in mse from photometry matching
    mse = calc_mse(training_set, template, bandpass_dict)
    dmse = (mse - mse0)/mse0 
    mse0 = mse

    mse_list = [mse]
    templates = [copy.deepcopy(template)]

    # start perturbations
    pertN = 0
    while abs(dmse) > dmse_stop:
        pertN += 1

        # perturb the template
        sol = perturb_template(training_set, template, bandpass_dict, w, Delta)
        template.flambda += sol
        templates.append(copy.deepcopy(template))
        
        # calculate the new mse and the fractional change
        mse = calc_mse(training_set, template, bandpass_dict)
        dmse = (mse - mse0)/mse0 
        mse0 = mse # update mse0
        mse_list.append(mse)

        if pertN == maxPerts:
            break
    
    result = [key,templates,mse_list]
    return result



def perturb_template(training_set, template, bandpass_dict, w=0.5, Delta=None):
    """
    Perturbs the template according to the photometry in the training set, which
    is a list of galaxy objects.

    bandpass_dict is the dictionary of filters used to observe the galaxies, and
    is used to calculate synthetic photometry for the template.
    w is the training ratio described in the paper and is used to calculate Delta.
    Values of order 1 work well.
    Delta can also be set manually.
    """
    

    wavelen = template.wavelen
    nbins = len(wavelen)
    widths = np.diff(wavelen)
    widths = np.append(widths,widths[-1])

    M = np.zeros((nbins,nbins))
    nu = np.zeros(nbins)
    sigmas = np.array([])

    for galaxy in training_set:

        template_ = copy.deepcopy(template)
        template_.redshift(galaxy.redshift)
        template_fluxes = template_.fluxlist(bandpass_dict, galaxy.filters)

        rn = np.array([np.interp(wavelen * (1+galaxy.redshift),
                            bandpass_dict[band].wavelen, 
                            bandpass_dict[band].R) for band in galaxy.filters])
        dlambda = widths * (1+galaxy.redshift)
        rn_dlambda = rn*dlambda.T

        gn = galaxy.fluxes - template_fluxes
        sigma_n = galaxy.flux_err/galaxy.fluxes
        gos2 = gn/sigma_n**2
        
        M  += np.sum( [np.outer(row,row)/sigma_n[i]**2 for i,row in enumerate(rn_dlambda)], axis=0 )
        nu += np.sum( (np.diag(gos2) @ rn_dlambda), axis=0 )

        sigmas = np.concatenate((sigmas,sigma_n))

    if Delta is None:
        Delta = np.mean(sigmas)*np.sqrt(len(template.wavelen)/(w*len(sigmas)))
        Delta = np.clip(Delta,0,0.05)
    M += np.identity(nbins)/Delta**2

    sol = np.linalg.solve(M,nu)

    return sol 


 
def calc_mse(training_set, template, bandpass_dict):
    """
    Calculates the weighted mean square error (wMSE) between a set of observed galaxy photometry
    and synthetic photometry for an SED template.
    """
    
    se = 0
    N = 0

    for galaxy in training_set:

        template_ = copy.deepcopy(template)
        template_.redshift(galaxy.redshift)
        template_fluxes = template_.fluxlist(bandpass_dict, galaxy.filters)

        N += len(galaxy.fluxes)
        se += sum( (galaxy.fluxes/galaxy.flux_err)**2 * (galaxy.fluxes - template_fluxes)**2 )

    mse = se/N
    return mse