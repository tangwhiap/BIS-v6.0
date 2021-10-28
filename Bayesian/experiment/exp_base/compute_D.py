#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from ...utils.netcdf_io import D_write

import numpy as np
from pdb import set_trace

def compute_D(method = "expLag", *args, **kwargs):
    exec("compute_D_" + method + "(*args, **kwargs)")

def compute_D_expLag(objExp, typeName, Lt, dim1TList, dim2TList):
    #print("Computing D matrix ...")
    #tlist = np.arange(n_tback).astype('int')
    LagT = []
    for it in dim1TList:
        lag_it = [(it - time).total_seconds()/3600 for time in dim2TList]
        LagT.append( np.abs(lag_it) )
    LagT = np.array(LagT)
    D = np.exp( - LagT / Lt )
    D_write(outFile = objExp.get_DFile_name(typeName), array_D = D, dim1TList = dim1TList, dim2TList = dim2TList)

