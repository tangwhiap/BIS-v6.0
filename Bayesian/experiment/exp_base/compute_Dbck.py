#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from ...utils.netcdf_io import D_write

import numpy as np

def compute_Dbck(objExp, dim1TList, dim2TList, method = "expLag", *args, **kwargs):

    loc = locals()
    exec("Dbck = compute_Dbck_" + method + "(objExp = objExp, dim1TList = dim1TList, dim2TList = dim2TList, *args, **kwargs)")
    Dbck = loc["Dbck"]
    D_write(objExp.get_DbckFile_name(), Dbck, dim1TList, dim2TList)

def compute_Dbck_expLag(objExp, LtBck, dim1TList, dim2TList):
    LagT = []
    for it in dim1TList:
        lag_it = [(it - time).total_seconds()/3600 for time in dim2TList]
        LagT.append( np.abs(lag_it) )
    LagT = np.array(LagT)
    Dbck = np.exp( - LagT / LtBck )
    return Dbck
