#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from ...utils.netcdf_io import get_ncvar, get_nctime, nc_addvar

import numpy as np
from pdb import set_trace

def background_optimize(objCore, objExp, recepTime, optTimeList):

    """
    Optimize the background values and its sigma values using bayesian scheme:
    * formulation:
        dBck = Qbck_opt.T * (HQHt + R + Qbck_c)^-1 * (z - Hx - bck)
        dQbck = Qbck_opt.T * (HQHt + R + Qbck_c)^-1 * Qbck_opt
    """
    #-- Get (z - Hx - bck(prior)) --#
    d = get_ncvar(objCore.get_dFile_name(recepTime = recepTime), "data")

    #-- Valid observation index --#
    obsValidIndex = ~np.isnan(d.flatten())

    #-- Remove NaN values from d --#
    d = d[obsValidIndex]

    #-- Compute Qbck_c in recepTime--#
    Qbck_c, _, _ = get_Qbck(dim1Trange = [recepTime], dim2Trange = [recepTime], objCore = objCore, objExp = objExp)

    #-- Get HQHt + R --#
    HQHt_R = get_ncvar(objCore.get_HQHtRFile_name(recepTime = recepTime), "data")

    #-- Compute (HQHt + R + Qbck_c)^-1 --#
    INV = np.linalg.inv((HQHt_R + Qbck_c)[obsValidIndex][:, obsValidIndex])
    
    #-- Get Qbck_opt --#
    Qbck_opt, _, qdim2TList = get_Qbck(dim1Trange = [recepTime], dim2Trange = optTimeList, objCore = objCore, objExp = objExp)

    #-- Remove the corresponding dimension where the observation is Nan. --#
    Qbck_opt = Qbck_opt[obsValidIndex]

    #-- Compute dBck --#
    dBck = np.matmul(Qbck_opt.T, np.matmul(INV, d))

    #-- Compute dQbck --#
    dQbck = np.matmul(Qbck_opt.T, np.matmul(INV, Qbck_opt))
    dQbck_diag = np.diag(dQbck)

    #-- Modify the background files and background sigma files for each time involved (dim2Trange) --#
    for time in list(set(qdim2TList)):

        #-- Find the index of values represent things in this time --#
        index = qdim2TList == time

        #-- Modify background files, add dBck into original background values --#
        nc_addvar(objExp.get_bckProcFile_name(recepTime = time), "data", dBck[index])
       
        #-- Compute the increment of background sigma values. --#
        sbck = get_ncvar(objExp.get_sbckProcFile_name(recepTime = time), "data").flatten()
        sbck2_opt = sbck ** 2 - dQbck_diag[index]
        sbck2_opt = np.where(sbck2_opt >= 0, sbck2_opt, 0) # If the variance value is less than zero, change it to zero.
        sbck_opt = np.sqrt(sbck2_opt)
        dsbck = sbck_opt - sbck

        #-- Modify background sigma files, add dBck into original background sigma values --#
        nc_addvar(objExp.get_sbckProcFile_name(recepTime = time), "data", dsbck.reshape(len(dsbck), 1))

    #Qbck = get_Qbck(dim1Trange = [recepTime], dim2Trange = [recepTime], objCore = objCore, objExp = objExp)




def get_Qbck(dim1Trange, dim2Trange, objCore, objExp):

    """
        dim1Trange: The time range of the first dimension of background Q matrix.
        dim2Trange: The time range of the second dimensions of background Q matrix.
        objCore: Object of optimization core.
        objExp: Object of experiment core.
    """

    #-- Convert type from pydatetime to string (yyyymmddHHMMSS) --#
    dim1Trange_str = [time.strftime("%Y%m%d%H%M%S") for time in dim1Trange]
    dim2Trange_str = [time.strftime("%Y%m%d%H%M%S") for time in dim2Trange]

    #-- Get D and E file name --#
    D_fileName = objExp.get_DbckFile_name()
    E_fileName = objExp.get_EbckFile_name()

    #-- Get background D matrix value --#
    D = get_ncvar(D_fileName, "data")

    #-- Get timelist from D file and convert to string (yyyymmddHHMMSS) --#
    dim1T = [time.strftime("%Y%m%d%H%M%S") for time in get_nctime(D_fileName, "dim1TList")]
    dim2T = [time.strftime("%Y%m%d%H%M%S") for time in get_nctime(D_fileName, "dim2TList")]

    #-- Get background E matrix value --#
    E = get_ncvar(E_fileName, "data")

    #-- Search the index of timelist in D file for each time step in dim1TList --#
    dim1Index = [dim1T.index(time) for time in dim1Trange_str]

    #-- Search the index of timelist in D file for each time step in dim2TList --#
    dim2Index = [dim2T.index(time) for time in dim2Trange_str]

    #-- Select the valid section of D matrix --#
    D_sel = D[dim1Index, dim2Index]

    #-- Compute background error correlation matrix --#
    corr = np.kron(D_sel, E)

    #-- Get background error (standart deviation) for each time in dim1Trange. --#
    #-- sback1List is a flattened array background error in all time. --#
    #-- qdim1TList is the timelist corresponding to sback1List. --#
    sbck1List = []
    qdim1TList = []
    for time in dim1Trange:
        sbck1 = get_ncvar(objExp.get_sbckProcFile_name(recepTime = time), "data").flatten()
        sbck1List.append(sbck1)
        qdim1TList += [time] * len(sbck1)
    sbck1List = np.array(sbck1List).flatten()
    qdim1TList = np.array(qdim1TList)

    #-- Get background error (standart deviation) for each time in dim2Trange. --#
    #-- sback2List is a flattened array background error in all time. --#
    #-- qdim2TList is the timelist corresponding to sback1List. --#
    sbck2List = []
    qdim2TList = []
    for time in dim2Trange:
        sbck2 = get_ncvar(objExp.get_sbckProcFile_name(recepTime = time), "data").flatten()
        sbck2List.append(sbck2)
        qdim2TList += [time] * len(sbck2)
    sbck2List = np.array(sbck2List).flatten()
    qdim2TList = np.array(qdim2TList)

    #-- Compute background error covariance matrix --#
    Qbck = np.matmul(np.matmul(np.diag(sbck1List), corr), np.diag(sbck2List))

    return Qbck, qdim1TList, qdim2TList


