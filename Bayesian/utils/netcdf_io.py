#!/usr/bin/env python


from .sparse_matrix import Sparse_to_nc 

import numpy as np
import netCDF4 as nc
import datetime as dtm

# nc.date2num(time, units = "seconds since 1970-01-01 00:00:00Z", calendar = "standard")
def get_ncvar(FileName, VARS):
    ncf = nc.Dataset(FileName, "r")
    if not isinstance(VARS, list):
        assert isinstance(VARS, str)
        VARS = [VARS]
    dataList = []
    for VAR in VARS:
        assert isinstance(VAR, str)
        data = ncf.variables[VAR][:].filled(np.nan)
        dataList.append(data)
    ncf.close()
    return dataList[0] if len(dataList) == 1 else tuple(dataList)

def get_ncdim(FileName, DIMS):
    ncf = nc.Dataset(FileName, "r")
    if not isinstance(DIMS, list):
        assert isinstance(DIMS, str)
        DIMS = [DIMS]
    dimList = []
    for DIM in DIMS:
        assert isinstance(DIM, str)
        dim = ncf.dimensions[DIM].size
        dimList.append(dim)
    ncf.close()
    return dimList[0] if len(dimList) == 1 else tuple(dimList)

#def get_ncattr(FileName, VAR = None, ATTR = None):

def get_nctime(FileName, VARS):
    ncf = nc.Dataset(FileName, "r")
    if not isinstance(VARS, list):
        assert isinstance(VARS, str)
        VARS = [VARS]
    dataList = []
    for VAR in VARS:
        assert isinstance(VAR, str)
        data = ncf.variables[VAR][:].filled(np.nan)
        tcalendar = ncf.variables[VAR].calendar
        tunits = ncf.variables[VAR].units
        timeList = nc.num2date(data, units = tunits, calendar = tcalendar)
        timeList = [dtm.datetime(time.year, time.month, time.day, time.hour, time.minute, time.second) for time in timeList]
        dataList.append(timeList)
    ncf.close()
    return dataList[0] if len(dataList) == 1 else tuple(dataList)

        
def X_Sigma_Indicator_write(outFile, emiss, LON, LAT):
    ncf_XS = nc.Dataset(outFile, "w", format = "NETCDF4")
    assert emiss.shape == LON.shape and LON.shape == LAT.shape
    nY, nX = emiss.shape
    ncf_XS.createDimension("dim1", nY * nX)
    ncf_XS.createDimension("dim2", 1)
    ncf_XS.createDimension("dimY", nY)
    ncf_XS.createDimension("dimX", nX)
    emissFlatten = emiss.flatten()
    ncf_XS.createVariable("data", "f", ("dim1", "dim2"))[:] = emissFlatten.reshape(len(emissFlatten), 1)
    ncf_XS.createVariable("LAT", "f", ("dimY", "dimX"))[:] = LAT
    ncf_XS.createVariable("LON", "f", ("dimY", "dimX"))[:] = LON
    ncf_XS.createVariable("emiss", "f", ("dimY", "dimX"))[:] = emiss
    ncf_XS.close()

def output_write(outFile, emiss, LON, LAT):
    print("Writing to " + outFile)
    ncf = nc.Dataset(outFile, "w", format = "NETCDF4")
    assert emiss.shape == LON.shape and LON.shape == LAT.shape
    nY, nX = emiss.shape
    ncf.createDimension("dimY", nY)
    ncf.createDimension("dimX", nX)
    ncf.createVariable("LAT", "f", ("dimY", "dimX"))[:] = LAT
    ncf.createVariable("LON", "f", ("dimY", "dimX"))[:] = LON
    ncf.createVariable("emiss", "f", ("dimY", "dimX"))[:] = emiss
    ncf.close()

def obs_write(outFile, obs):
    assert isinstance(obs, np.ndarray)
    if len(obs.shape) == 1:
        obs = obs.reshape(len(obs), 1)
    ncf_obs = nc.Dataset(outFile, "w", format = "NETCDF4")
    nY, nX = obs.shape
    ncf_obs.createDimension("dim1", nY)
    ncf_obs.createDimension("dim2", nX)
    ncf_obs.createVariable("data", "f", ("dim1", "dim2"))[:] = obs
    ncf_obs.close()

def D_write(outFile, array_D, dim1TList, dim2TList):
    assert isinstance(array_D, np.ndarray)
    ncf_D = nc.Dataset(outFile, "w", format = "NETCDF4")
    nY, nX = array_D.shape
    assert len(dim1TList) == nY
    assert len(dim2TList) == nX
    ncf_D.createDimension("dim1", nY)
    ncf_D.createDimension("dim2", nX)
    ncf_D.createVariable("data", "f", ("dim1", "dim2"))[:] = array_D

    tunits = "seconds since 1970-01-01 00:00:00Z"
    tcalendar = "standard"

    ncvar_dim1TList = ncf_D.createVariable("dim1TList", "i8", ("dim1"))
    ncvar_dim1TList[:] = nc.date2num(dim1TList, units = tunits, calendar = tcalendar)
    ncvar_dim1TList.setncattr("units", tunits)
    ncvar_dim1TList.setncattr("calendar", tcalendar)

    ncvar_dim2TList = ncf_D.createVariable("dim2TList", "i8", ("dim2"))
    ncvar_dim2TList[:] = nc.date2num(dim2TList, units = tunits, calendar = tcalendar)
    ncvar_dim2TList.setncattr("units", tunits)
    ncvar_dim2TList.setncattr("calendar", tcalendar)
    
    ncf_D.close()

def nc_write(outFile, arr):
    assert isinstance(arr, np.ndarray)
    ncf_arr = nc.Dataset(outFile, "w", format = "NETCDF4")
    nY, nX = arr.shape
    ncf_arr.createDimension("dim1", nY)
    ncf_arr.createDimension("dim2", nX)
    ncf_arr.createVariable("data", "f", ("dim1", "dim2"))[:] = arr
    ncf_arr.close()


def nc_addvar(modFile, varName, arr):
    assert isinstance(arr, np.ndarray)
    ncf_mod = nc.Dataset(modFile, "r+")
    nY, nX = ncf_mod.variables[varName].shape
    assert (nY, nX) == arr.shape, "Error! Conflict shape."
    ncf_mod.variables[varName][:] += arr
    ncf_mod.close()

#def HQ_write(outFile, spr_HQ):
#    Sparse_to_nc(spr_HQ, outFile)
