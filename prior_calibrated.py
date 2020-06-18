# This is the prior calibrated to our training set
# the parameters are taken from calibrate_prior.ipynb
# This is written to mimic prior_hdfn_gen.py in the bpz folder

import numpy as np
from scipy.special import gamma

def function(z, m, nt):

    # nt Templates = nell Elliptical + nsp Spiral + nim Im/Starburst
    try:
        nell, nsp, nim = nt
    except:
        nell = 1 # default to 1 elliptical
        nsp = 2 # 2 spiral
        nim = nt - nell - nsp # rest Im/SB
    nt = nell, nsp, nim
    N = sum(nt)

    # functional form of p(T|m0) for El and Im
    logistic = lambda m0,L,k,M,c: L/(1+np.exp(-k*(m0-M))) + c
    #logistic = lambda m0,L,k: L*np.exp(-k*(m0-20))

    # El and Im p(T|m0) parameters from logistic fits
    L_El, k_El, M_El, c_El = 0.4483, -1.4479, 20.9492, 0.0069
    L_Im, k_Im, M_Im, c_Im = 0.8512,  1.1952, 22.5595, 0.0885
    #L_El, k_El = 0.35, 0.147
    #L_Sp, k_Sp = 0.50, 0.450

    # p(T|m0)
    pT_El = logistic(m,L_El,k_El,M_El,c_El)
    pT_Im = logistic(m,L_Im,k_Im,M_Im,c_Im)
    pT_Sp = 1 - (pT_El + pT_Im)
    #pT_El = logistic(m,L_El,k_El)
    #pT_Sp = logistic(m,L_Sp,k_Sp)
    #pT_Im = 1 - (pT_El + pT_Sp)

    # El, Sp, and Im p(z|T,m0) parameters from Likelihood maximization
    a_El, z0_El, km_El = 3.8791, 0.4844, 0.1187
    a_Sp, z0_Sp, km_Sp = 3.3966, 0.4925, 0.1245
    a_Im, z0_Im, km_Im = 2.2086, 0.3635, 0.1277
    #a_El, z0_El, km_El = 2.46, 0.431, 0.091
    #a_Sp, z0_Sp, km_Sp = 1.81, 0.390, 0.0636
    #a_Im, z0_Im, km_Im = 0.91, 0.063, 0.123

    # pz = p(z|T,m0)
    zm = lambda m0,z0,km: z0 + km*(m0 - 20)
    norm = lambda m0,a,z0,km: 1/a * gamma((a+1)/a) * zm(m0,z0,km)**(a+1)
    pz = lambda z,m0,a,z0,km: 1/norm(m0,a,z0,km) * z**a * np.exp(-(z/zm(m0,z0,km))**a)

    pz_El = lambda z: pz(z,m,a_El,z0_El,km_El)
    pz_Sp = lambda z: pz(z,m,a_Sp,z0_Sp,km_Sp)
    pz_Im = lambda z: pz(z,m,a_Im,z0_Im,km_Im)

    # prob = p(z,T|m0) = p(T|m0)*p(z|T,m0)
    prob_El = lambda z: pT_El * pz_El(z)
    prob_Sp = lambda z: pT_Sp * pz_Sp(z)
    prob_Im = lambda z: pT_Im * pz_Im(z)

    # now calculate the array p_i that bpz needs
    # p_i is an len(z) x N array
    # Each row is a different z value, each column
    # a different SED template
    # p_i(j,k) = p(z_j,T_k|m0)/degen(T_k)
    # m0 is a constant across the whole array
    # degen(T_k) is the "degeneracy" of template k
    # i.e. how many templates are in the same spectral
    # class as template k? This means we are assuming that 
    # all templates of a given class are equally probable
    # note that I actually assemble p_i.T, so at the end,
    # I have to transpose to get the real p_i
    p_i = np.tile(z,(N,1))
    for i in range(nt[0]):
        p_i[i] = prob_El(p_i[i])/nt[0]
    for i in range(nt[1]):
        j = i + nt[0]
        p_i[j] = prob_Sp(p_i[j])/nt[1]
    for i in range(nt[2]):
        j = i + nt[0] + nt[1]
        p_i[j] = prob_Sp(p_i[j])/nt[2]

    p_i = p_i.T 
    
    # now hand the prior array back to bpz!
    return p_i