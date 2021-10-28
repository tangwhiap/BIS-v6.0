#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...
from ...utils.netcdf_io import get_ncvar, nc_write

import multiprocessing as mtp

def compute_concentration_time(args):

    objExp, recepTime, Xtype, outFile_Hx, outFile_HxBck = args

    assert Xtype.lower() in ["prior", "proc", "post"]
    
    timeList = objExp.objCore.backtime_j2i(recepTime)
    classX = objExp.indicators["X"]
    sumHx = objExp.objCore.Hx(timeList = timeList, recepTime = recepTime, Xtype = Xtype, class_XIndicator = classX, objExp = objExp)
    Hx = objExp.objCore.Hx2obs(sumHx = sumHx)
    if objExp.hasBCK:
        if Xtype.lower() == "prior":
            bck = get_ncvar(objExp.get_bckPriorFile_name(recepTime = recepTime), "data")
        else:
            bck = get_ncvar(objExp.get_bckProcFile_name(recepTime = recepTime), "data")
    else:
        bck = 0
    nc_write(outFile_Hx, Hx)
    nc_write(outFile_HxBck, Hx + bck)

def compute_concentration(objExp, Xtype, outFile_Hx_Fun, outFile_HxBck_Fun):

    assert Xtype.lower() in ["prior", "proc", "post"]
    parallelArgs = [(objExp, time, Xtype, outFile_Hx_Fun(recepTime = time), outFile_HxBck_Fun(recepTime = time)) for time in objExp.objCore.objIter.obsTimeList]
    pool = mtp.Pool(objExp.objCore.nProc)
    pool.map(compute_concentration_time, parallelArgs)
    pool.close()
    pool.join()

def compute_final_concentration(objExp):
    Xtype = "Post"
    outFile_Hx_Fun = objExp.get_finalHxFile_name
    outFile_HxBck_Fun = objExp.get_finalHxBckFile_name
    compute_concentration(objExp, Xtype, outFile_Hx_Fun, outFile_HxBck_Fun)

def compute_init_concentration(objExp):
    Xtype = "Prior"
    outFile_Hx_Fun = objExp.get_initHxFile_name
    outFile_HxBck_Fun = objExp.get_initHxBckFile_name
    compute_concentration(objExp, Xtype, outFile_Hx_Fun, outFile_HxBck_Fun)


