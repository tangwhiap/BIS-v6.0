#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from ...utils.distance import DIS_point_to_points
from ...utils.netcdf_io import nc_write

import numpy as np
from pdb import set_trace

def compute_Ebck(objExp, method = "expdis2D", *args, **kwargs):

    loc = locals()
    exec("Ebck = compute_Ebck_" + method + "(objExp, *args, **kwargs)")
    Ebck = loc["Ebck"]
    nc_write(objExp.get_EbckFile_name(), Ebck)

def compute_Ebck_expdis2D(objExp, LsBck, lonList, latList):
    disMat = station_dismat(lonList, latList)
    Ebck = np.exp( - disMat / LsBck )
    return Ebck
    
    
def station_dismat(lonList, latList):
   
    assert len(lonList) == len(latList)
    Ns = len(lonList)
    dismat = np.full( (Ns, Ns), np.nan )
    for ista, (sta_lon, sta_lat) in enumerate(zip(lonList, latList)):
        dis_list = DIS_point_to_points(sta_lon, sta_lat, lonList, latList)
        dismat[ista, :] = dis_list
    return dismat
 
