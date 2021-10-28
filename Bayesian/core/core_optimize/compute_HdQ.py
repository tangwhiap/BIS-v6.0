#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 09/2021
#   ...

from ..indicators.indicator_operator import IndicatorOperator_Zero
from ..indicators.temp_indicator import HQ_Indicator, HdQ_Indicator, INVHQ_Indicator, HiHQti_Indicator
from .indicators import H_Indicator
from ...utils.sparse_matrix import Sparse_to_nc
from ...utils.netcdf_io import nc_write

from pdb import set_trace

def compute_HdQj(kwargs):

    assert "history_tList" in kwargs
    #assert "ti_range" in kwargs
    assert "tj" in kwargs
    assert "recepTime" in kwargs
    assert "objCore" in kwargs
    assert "objExp" in kwargs

    history_tList = kwargs["history_tList"]
    #ti_range = kwargs["ti_range"]
    tj = kwargs["tj"]
    recepTime = kwargs["recepTime"]
    objCore = kwargs["objCore"]
    objExp = kwargs["objExp"]

    #print("Processing tj = " + tj.strftime("%Y-%m-%d %H:%M:%S"))
    sectors = objExp.sectors


    HdQj_IndOp = IndicatorOperator_Zero()
    for hisTime in history_tList:
        try:
            INVHQj_Ind = INVHQ_Indicator(sectors = sectors, get_INVHQFile_name = objCore.get_INVHQFile_name, jtime = tj, recepTime = recepTime)
        except FileNotFoundError:
            #print("continue for hisTime = " + hisTime.strftime("%Y-%m-%d %H:%M:%S"))
            continue

        #Hi_HQti_Ind = HiHQti_Indicator(sectors = sectors, recepTime = recepTime, hisTime = hisTime)
        Hi_HQti_Ind = HiHQti_Indicator(sectors = sectors, get_HiHQtiFile_name = objCore.get_HiHQtiFile_name, recepTime = recepTime, hisTime = hisTime)
        HdQj_IndOp = HdQj_IndOp + Hi_HQti_Ind.to_operator() * INVHQj_Ind.to_operator()

    if not isinstance(HdQj_IndOp, IndicatorOperator_Zero):
        #print("Compute HdQ on tj = " +  tj.strftime("%Y-%m-%d %H:%M:%S") + " for recepTime = " + recepTime.strftime("%Y-%m-%d %H:%M:%S"))
        HdQj_Ind = HdQ_Indicator(objIndOp = HdQj_IndOp, jtime = tj, recepTime = recepTime)
        HdQj_Ind.to_file(get_HdQFile_name = objCore.get_HdQFile_name)
        #for sectorName in sectors:
        #    nc_write(outFile = objCore.get_HdQFile_name(jtime = tj, recepTime = recepTime, sectorName = sectorName), arr = HdQj_IndOp.data[sectorName])
    else:
        pass
        #print("Warning! Skip HdQ on tj = " + tj.strftime("%Y-%m-%d %H:%M:%S") + " for recepTime = " + recepTime.strftime("%Y-%m-%d %H:%M:%S"))


def compute_HdQ(objCore, objExp, recepTime, tj_range, backtime_j2i, backtime_j2i_kwargs = {}):

    #print("Computing HdQ begin ...")
    sectors = objExp.sectors 
    ti_range = backtime_j2i(recepTime, **backtime_j2i_kwargs)

    history_tList = []

    for itime in objCore.objIter.iterTimeList:
        if itime < recepTime:
            history_tList.append(itime)
    #print("history_tList:")
    #print(history_tList)

    for hisTime in history_tList:
        sum_Hi_HQti_IndOp = IndicatorOperator_Zero()
        for ti in ti_range:
            try:
                HQi_Ind = HQ_Indicator(sectors = sectors, get_HQFile_name = objCore.get_HQFile_name, recepTime = hisTime, jtime = ti)
            except FileNotFoundError:
                continue
            Hi_Ind = H_Indicator(time = ti, recepTime = recepTime, typeToSector = objExp.typeToSector_H, fun_get_filename = objExp.get_HFile_name)

            Hi_HQti_IndOp = Hi_Ind.to_operator() * HQi_Ind.trans().to_operator()
            sum_Hi_HQti_IndOp  = sum_Hi_HQti_IndOp + Hi_HQti_IndOp

        #for sectorName in sectors:
            #nc_write(outFile = objCore.get_HiHQiFile_name(recepTime = recepTime, hisTime = hisTime, sectorName = sectorName), arr = sum_Hi_HQi_IndOp.data[sectorName])
            #Sparse_to_nc(spr = sum_Hi_HQi_IndOp.data[sectorName], FileDirName = objCore.get_HiHQiFile_name(recepTime = recepTime, hisTime = hisTime, sectorName = sectorName))

        #print("Write HiHQti, recepTime = " + recepTime.strftime("%Y-%m-%d %H:%M:%S") + ", hisTime = " + hisTime.strftime("%Y-%m-%d %H:%M:%S"))
        if not isinstance(sum_Hi_HQti_IndOp, IndicatorOperator_Zero):
            Hi_HQti_Ind = HiHQti_Indicator(objIndOp = sum_Hi_HQti_IndOp, recepTime = recepTime, hisTime = hisTime)
            Hi_HQti_Ind.to_file(get_HiHQtiFile_name = objCore.get_HiHQtiFile_name)

    parallelArgs = []
    for tj in tj_range:
        kwargs = {
            "history_tList": history_tList,
            "tj": tj, "recepTime": recepTime,
            "objCore": objCore, "objExp": objExp,
        }
        parallelArgs.append(kwargs)

    for kwargs in parallelArgs:
        compute_HdQj(kwargs)





    

