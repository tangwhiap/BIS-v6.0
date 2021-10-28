#!/usr/bin/env python
# Module of sigma computing functions.

# Authors:
#   Wenhan TANG - 07/2021
#   ...


import numpy as np
from pdb import set_trace

def sigma_constant(const, LON, LAT, **kwargs):
    assert LON.shape == LAT.shape
    Sigma = LON * 0 + const
    return Sigma, LON, LAT

def sigma_emiss_ratio(ratio, time, sectorName, getFun_prior, **kwargs):

    emiss, LON, LAT = getFun_prior(time = time, sectorName = sectorName) 
    sigma = np.abs(emiss) * ratio
    return sigma, LON, LAT



mean_difference = None
LON_mean_difference = None
LAT_mean_difference = None
def sigma_mean_difference(time, sectorName, getFun_prior, getFun_truth, Start, End, dt, time_vary = False, **kwargs):
    global mean_difference, LON_mean_difference, LAT_mean_difference
    if mean_difference is None:
        differenceList = []
        current = Start
        while(current <= End):
            prior, LON_prior, LAT_prior = getFun_prior(time = time, sectorName = sectorName)    
            truth, LON_truth, LAT_truth = getFun_truth(time = time, sectorName = sectorName)
            assert np.sum(np.abs( LON_prior - LON_truth) ) < 1e-6
            assert np.sum(np.abs( LAT_prior - LAT_truth) ) < 1e-6
            sigma_const = np.sqrt(np.abs(truth **2 - prior **2))
            differenceList.append(sigma_const)
            current += dt
        mean_difference = np.array(differenceList).mean()
        LON_mean_difference = LON_prior; LAT_mean_difference = LAT_prior
    return sigma_constant(mean_difference, LON_mean_difference, LAT_mean_difference)

def sigma_difference(time, sectorName, getFun_prior, getFun_truth, Start, End, dt, time_vary = True, **kwargs):
    prior, LON_prior, LAT_prior = getFun_prior(time = time, sectorName = sectorName)
    truth, LON_truth, LAT_truth = getFun_truth(time = time, sectorName = sectorName)
    assert np.sum(np.abs( LON_prior - LON_truth) ) < 1e-6
    assert np.sum(np.abs( LAT_prior - LAT_truth) ) < 1e-6
    sigma = np.sqrt(np.abs(prior**2 - truth**2))
    return sigma, LON_prior, LAT_prior



funNameDic = {
    "constant": sigma_constant,
    "mean_difference": sigma_mean_difference,
    "difference": sigma_difference,
    "emiss_ratio": sigma_emiss_ratio,
}

def compute_sigma(time = None, funName = None, sectorName = None, getFun_prior = None, getFun_truth = None, LON = None, LAT = None, **kwargs):
    assert funName in funNameDic
    compute_sigma_kwargs = {
        "time": time, "sectorName": sectorName, "getFun_prior": getFun_prior, "getFun_truth": getFun_truth, "LON": LON, "LAT": LAT,
    }
    compute_sigma_kwargs = {**compute_sigma_kwargs, **kwargs}
    func = funNameDic[funName]
    return func(**compute_sigma_kwargs)



