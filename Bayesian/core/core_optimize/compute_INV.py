#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from .background import get_Qbck
from ..indicators.temp_indicator import HQ_Indicator, HQHt_Indicator
from ...utils.netcdf_io import get_ncvar, nc_write

import numpy as np
from pdb import set_trace

def compute_INV(objCore, objExp, recepTime):

    #HQHt_Ind = HQHt_Indicator(sectors = objExp.sectors, get_HQHtFile_name = objCore.get_HQHtFile_name, recepTime = recepTime)

    #HQHt = HQHt_Ind.data_sum()

    #R = get_ncvar(objExp.get_RFile_offline_name(recepTime = recepTime), "data")

    HQHt_R = get_ncvar(objCore.get_HQHtRFile_name(recepTime = recepTime), "data")

    if objExp.hasBCK:
        #Qbck = get_ncvar(objExp.get_QbckFile_name(recepTime = recepTime), "data")
        Qbck, _, _ = get_Qbck(dim1Trange = [recepTime], dim2Trange = [recepTime], objCore = objCore, objExp = objExp)
        HQHt_R = HQHt_R + Qbck

    INV = np.linalg.inv(HQHt_R)

    nc_write(objCore.get_INVFile_name(recepTime = recepTime), INV)
