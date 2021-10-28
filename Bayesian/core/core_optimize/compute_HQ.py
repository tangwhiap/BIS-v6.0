#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

#from .indicators.sectors_indicator import Sigma_Indicator
from .indicators import H_Indicator, D_Indicator, E_Indicator
from ..indicators.indicator_operator import IndicatorOperator_Zero, IndicatorOperator_One
from ..indicators.temp_indicator import HQ_Indicator, HdQ_Indicator, HQHt_Indicator

import multiprocessing as mtp
from pdb import set_trace
import time

def compute_HQj(kwargs):

    """
        Used to caculate Hi multiply with Qij.

        * Formula:
            HQj = sum(dij * Hi * diag(sigma_i)) * E * diag(sigma_j)
            ( sum for i map the ti_range, j represents the argument: tj )

        * Reference:
            Yadav, V. and Michalak, A. M.: Improving computational efficiency in large linear inverse problems: an example from carbon dioxide flux estimation, Geosci. Model Dev., 6, 583–590, 2013.
    """

    st = time.time()

    assert "ti_range" in kwargs
    assert "tj" in kwargs
    assert "recepTime" in kwargs
    #assert "get_HQFile_name" in kwargs
    #assert "Sigma_Indicator" in kwargs
    assert "objCore" in kwargs
    assert "objExp" in kwargs
    assert "Stype" in kwargs

    ti_range = kwargs["ti_range"]
    tj = kwargs["tj"]
    recepTime = kwargs["recepTime"]
    objCore = kwargs["objCore"]
    objExp = kwargs["objExp"]
    Stype = kwargs["Stype"]
    assert Stype in ["Proc", "Prior"]

    sectors = objExp.sectors

    get_HQFile_name = objCore.get_HQFile_name

    Sigma_Indicator = objCore.ExpIndicators["Sigma"]

    # Get sigma_j
    sigmaInd_tj = Sigma_Indicator(time = tj, Stype = "Proc", objExp = objExp)

    # Get D matrix
    DInd = D_Indicator(typeToSector = objExp.typeToSector_D, fun_get_filename = objExp.get_DFile_name)

    # Get E matrix
    EInd = E_Indicator(typeToSector = objExp.typeToSector_E, fun_get_filename = objExp.get_EFile_name)

    # Define the inital value for sum
    dij_Hi_Si_IndOp = IndicatorOperator_Zero() 

    for ti in ti_range: # loop for sum

        # Caculation of Hi * diag(sigma_i)
        HInd = H_Indicator(time = ti, recepTime = recepTime, typeToSector = objExp.typeToSector_H, fun_get_filename = objExp.get_HFile_name)
        sigmaInd_ti = Sigma_Indicator(time = ti, Stype = Stype, objExp = objExp)
        Hi_Si_IndOp = HInd.to_operator() * sigmaInd_ti.diag_out()

        # Caculation of dij * Hi * diag(sigma_i)
        dij_IndOp = DInd[ti, tj]
        dij_Hi_Si_IndOp = dij_Hi_Si_IndOp + dij_IndOp * Hi_Si_IndOp

    # Caculation of sum(dij * Hi * diag(sigma_i)) * E * diag(sigma_j)
    HQj_IndOp = dij_Hi_Si_IndOp * EInd.to_operator() * sigmaInd_tj.diag_out()

    if objCore.myConfig["isCompute_dQ"]:
        try:
            HdQj_Ind = HdQ_Indicator(sectors = sectors, jtime = tj, recepTime = recepTime, get_HdQFile_name = objCore.get_HdQFile_name)
            HQj_IndOp = HQj_IndOp - HdQj_Ind.trans_to_sparse().to_operator()

        except FileNotFoundError:
            pass

    HQ = HQ_Indicator(jtime = tj, recepTime = recepTime, objIndOp = HQj_IndOp)

            

    HQ.to_file(get_HQFile_name = get_HQFile_name)

    ed = time.time()
    tConsume = ed - st
    print(tj, "-> time consume: " + "%4.2f" % tConsume + " s")


def compute_HQ(objCore, objExp, recepTime, tj_range, backtime_j2i, backtime_j2i_kwargs = {}, Stype = "Proc"):

    assert Stype in ["Proc", "Prior"]
    parallelArgs = []
    for tj in tj_range:
        #print("Computing HQ: ", tj)
        HQj_kwargs = {
            "ti_range": backtime_j2i(recepTime, **backtime_j2i_kwargs),
            "tj": tj, "recepTime": recepTime, "objCore": objCore,
            "objExp": objExp, "Stype": Stype,
        }
        parallelArgs.append(HQj_kwargs)
        #compute_HQj(HQj_kwargs)
    pool = mtp.Pool(objCore.nProc)
    pool.map(compute_HQj, parallelArgs)
    pool.close()
    pool.join()

def compute_HQHt(objCore, objExp, recepTime, backtime_j2i, backtime_j2i_kwargs = {}):
    ti_range = backtime_j2i(recepTime, **backtime_j2i_kwargs)
    HQHtsum_IndOp = IndicatorOperator_Zero() 
    for ti in ti_range:
        try:
            HQ_Ind = HQ_Indicator(sectors = objExp.sectors, get_HQFile_name = objCore.get_HQFile_name, jtime = ti, recepTime = recepTime)
        except FileNotFoundError:
            print("Warning! Skip ti = " + ti.strftime("%Y-%m-%d %H:%M:%S"))
            continue
        H_Ind = H_Indicator(time = ti, recepTime = recepTime, typeToSector = objExp.typeToSector_H, fun_get_filename = objExp.get_HFile_name)
        #H_Ind.trans()
        HQHtsum_IndOp = HQHtsum_IndOp + HQ_Ind.to_operator() * H_Ind.trans().to_operator()
    HQHt_Ind = HQHt_Indicator(objIndOp = HQHtsum_IndOp, recepTime = recepTime)
    HQHt_Ind.to_file(get_HQHtFile_name = objCore.get_HQHtFile_name)


