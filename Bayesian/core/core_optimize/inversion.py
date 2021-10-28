#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from .background import get_Qbck
from ..indicators.temp_indicator import HQ_Indicator, HQHt_Indicator, INVHQ_Indicator
from ...utils.netcdf_io import get_ncvar, nc_write
from pdb import set_trace

import numpy as np
import multiprocessing as mtp

def compute_INV(objCore, objExp, recepTime, obsValidIndex):

    #HQHt_Ind = HQHt_Indicator(sectors = objExp.sectors, get_HQHtFile_name = objCore.get_HQHtFile_name, recepTime = recepTime)

    #HQHt = HQHt_Ind.data_sum()

    #R = get_ncvar(objExp.get_RFile_offline_name(recepTime = recepTime), "data")

    HQHt_R = get_ncvar(objCore.get_HQHtRFile_name(recepTime = recepTime), "data")

    if objExp.hasBCK:
        #Qbck = get_ncvar(objExp.get_QbckFile_name(recepTime = recepTime), "data")
        Qbck, _, _ = get_Qbck(dim1Trange = [recepTime], dim2Trange = [recepTime], objCore = objCore, objExp = objExp)
        HQHt_R = HQHt_R + Qbck

    HQHt_R = HQHt_R[obsValidIndex][:, obsValidIndex]

    INV = np.linalg.inv(HQHt_R)

    nc_write(objCore.get_INVFile_name(recepTime = recepTime), INV)

def compute_HQt_INV_d(args):
    objCore, objExp, recepTime, time, INV, z_Hx, obsValidIndex = args
    sectors = objExp.sectors
    HQ = HQ_Indicator(sectors = sectors, get_HQFile_name = objCore.get_HQFile_name, recepTime = recepTime, jtime = time)
    #INV_d = np.matmul(INV[obsValidIndex][:, obsValidIndex], z_Hx[obsValidIndex])
    INV_d = np.matmul(INV, z_Hx[obsValidIndex])
    HQt_INV_d_IndOp = HQ.trans().to_operator()[:, obsValidIndex] * INV_d
    HQt_INV_d_Ind = objExp.indicators["dX"](HQt_INV_d_IndOp, time)
    HQt_INV_d_Ind.modify_file(objExp.get_procFile_name)
    #print(HQt_INV_d_Ind.data["total"].mean())



def compute_INVHQ(args):

    objCore, objExp, recepTime, tj, obsValidIndex = args
    sectors = objExp.sectors
    INV = get_ncvar(objCore.get_INVFile_name(recepTime = recepTime), "data")
    #for tj in tj_range:
    HQj = HQ_Indicator(sectors = sectors, get_HQFile_name = objCore.get_HQFile_name, recepTime = recepTime, jtime = tj)

    # Note: Detials about this problem:
    # https://stackoverflow.com/questions/38229953/array-and-rmul-operator-in-python-numpy
    #     --- Wenhan TANG - 09/2021
    #INVHQj_IndOp = INV * HQj.to_operator()
    INVHQj_IndOp = HQj.to_operator()[obsValidIndex, :].__rmul__(INV)

    INVHQj_Ind = INVHQ_Indicator(jtime = tj, recepTime = recepTime, objIndOp = INVHQj_IndOp)
    INVHQj_Ind.to_file(get_INVHQFile_name = objCore.get_INVHQFile_name)
    

def compute_dSigma(args):

    objCore, objExp, recepTime, time, obsValidIndex = args
    sectors = objExp.sectors
    HQ_Ind = HQ_Indicator(sectors = sectors, get_HQFile_name = objCore.get_HQFile_name, recepTime = recepTime, jtime = time)
    INV = get_ncvar(objCore.get_INVFile_name(recepTime = recepTime), "data")

    HQ_Ind = HQ_Ind[obsValidIndex]
    #INV = INV[obsValidIndex][:, obsValidIndex]

    resData = {}
    for sector in HQ_Ind.typeNames:
        nY, nX = HQ_Ind.dim1Dict[sector], HQ_Ind.dim2Dict[sector]
        resData[sector] = np.full( (nX, 1), np.nan)

        HQData = {}
        for sector in HQ_Ind.data:
            HQData[sector] = HQ_Ind.data[sector].toarray()

        for ix in range(nX):
            HQ_ix = HQData[sector][:, ix]
            HQt_ix = HQ_ix.T
            temp = np.matmul(HQt_ix, INV)
            temp = np.matmul(temp, HQ_ix)
            resData[sector][ix, :] = temp

            #HQ_ix = HQ_Ind.data[sector][:, ix]
            #HQt_ix = HQ_ix.T
            #temp = HQt_ix * INV
            #temp = temp * HQ_ix
            #set_trace()
            #resData[sector][ix, :] = temp


    dSigma2_Ind = objExp.indicators["dSigma"](resData, time, isSquare = True)

    dSigma2_Ind.modify_file(objExp.get_sigmaProcFile_name)

    dSigma2_Ind.to_file(objCore.get_dsigma2File_name, recepTime)

def inversion(objCore, objExp, recepTime, tj_range):
    
    z_Hx = get_ncvar(objCore.get_dFile_name(recepTime = recepTime), "data")
    obsValidIndex = ~np.isnan(z_Hx.flatten())

    compute_INV(objCore = objCore, objExp = objExp, recepTime = recepTime, obsValidIndex = obsValidIndex)
    INV = get_ncvar(objCore.get_INVFile_name(recepTime = recepTime), "data")
    
    #INV_d = np.matmul(INV, z_Hx)
    #print(obsValidIndex)

    print("Inversion ...")
    parallelArgs = []
    for tj in tj_range:
        #compute_HQt_INV_d(objCore, objExp, recepTime, tj, INV, z_Hx, obsValidIndex)
        parallelArgs.append((objCore, objExp, recepTime, tj, INV, z_Hx, obsValidIndex))
    pool = mtp.Pool(objCore.nProc)
    pool.map(compute_HQt_INV_d, parallelArgs)
    pool.close()
    pool.join()
             
    print("Computing reduced uncertainty ...")
    parallelArgs = []
    for tj in tj_range:
        #print("sigma: ", tj)
        #compute_dSigma(objCore, objExp, recepTime, tj, obsValidIndex)
        parallelArgs.append((objCore, objExp, recepTime, tj, obsValidIndex))
    pool = mtp.Pool(objCore.nProc)
    pool.map(compute_dSigma, parallelArgs)
    pool.close()
    pool.join()

    print("Computing INV*HQ ...")
    parallelArgs = []
    for tj in tj_range:
        #compute_INVHQ(objCore, objExp, recepTime, tj)
        parallelArgs.append((objCore, objExp, recepTime, tj, obsValidIndex))
    pool = mtp.Pool(objCore.nProc)
    pool.map(compute_INVHQ, parallelArgs)
    pool.close()
    pool.join()
    
    Hx = objCore.compute_Hx(recepTime, backtime_j2i = objCore.backtime_j2i, Xtype = "Proc", class_XIndicator = objExp.indicators["X"], objExp = objExp)
    Hx = objCore.Hx2obs(sumHx = Hx)
    nc_write(objCore.get_procHxFile_name(recepTime = recepTime), Hx)

    if objExp.hasBCK:
        bck = get_ncvar(objExp.get_bckProcFile_name(recepTime = recepTime), "data")
        nc_write(objCore.get_procHxBckFile_name(recepTime = recepTime), Hx + bck)
    else:
        nc_write(objCore.get_procHxBckFile_name(recepTime = recepTime), Hx)

