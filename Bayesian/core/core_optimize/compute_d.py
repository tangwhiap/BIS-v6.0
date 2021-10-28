#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

#from .compute_Hx import compute_Hx, Hx2obs

from ...utils.netcdf_io import get_ncvar, nc_write

from pdb import set_trace

def compute_d(objCore, objExp, recepTime, bckPrior = True):

    classX = objExp.indicators["X"]

    Hx = objCore.compute_Hx(recepTime, backtime_j2i = objCore.backtime_j2i, Xtype = "Proc", class_XIndicator = classX, objExp = objExp)

    Hx = objCore.Hx2obs(sumHx = Hx)

    z = get_ncvar(objExp.get_obsFile_name(recepTime = recepTime), "data")


    if objExp.hasBCK:
        bck = get_ncvar(objExp.get_bckProcFile_name(recepTime = recepTime), "data")
        d = z - (Hx + bck)
    else:
        d = z - Hx

    if bckPrior:
        nc_write(objCore.get_priorHxFile_name(recepTime = recepTime), Hx)

    #else:
    #    nc_write(objCore.get_procHxFile_name(recepTime = recepTime), Hx)

    nc_write(objCore.get_dFile_name(recepTime = recepTime), d)

    if bckPrior:
        if objExp.hasBCK:
            nc_write(objCore.get_priorHxBckFile_name(recepTime = recepTime), Hx + bck)
        else:
            nc_write(objCore.get_priorHxBckFile_name(recepTime = recepTime), Hx)

