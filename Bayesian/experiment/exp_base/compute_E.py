#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from ...utils.distance import DIS_point_to_points
from ...utils.show_progress import show_progress

import numpy as np
import netCDF4 as nc

def compute_E(method = "expdis2D", *args, **kwargs):
    exec("compute_E_" + method + "(*args, **kwargs)")

def compute_E_expdis2D(objExp, typeName, Ls, LON, LAT):

    LsPar = 3
    assert len(LON.shape) == 2 and LON.shape == LAT.shape
    myConfig = objExp.myConfig
    #Ls = myConfig["Ls_typeE"][typeName]
    nlat, nlon = LON.shape
    IX, IY = np.meshgrid(np.arange(nlon), np.arange(nlat))
    IG = np.arange(nlat * nlon).reshape(nlat, nlon)
    #LON_list = LON.flatten()
    #LAT_list = LAT.flatten()
    Ng = nlat * nlon
    ncf_E = nc.Dataset(objExp.get_EFile_name(typeName), "w", format = "NETCDF4")
    ncf_E.createDimension("DataDim", None)
    ncf_E.createDimension("PtrDim", Ng + 1)
    ncf_E.setncattr("ny", Ng)
    ncf_E.setncattr("nx", Ng)
    ncf_E.setncattr("Ls", Ls)
    ncvar_data = ncf_E.createVariable("data", "f", ("DataDim"))
    ncvar_indices = ncf_E.createVariable("indices", "i8", ("DataDim"))
    ncvar_indptr = ncf_E.createVariable("indptr", "i8", ("PtrDim"))


    Pt_DataDim = 0
    ncvar_indptr[0] = 0

    ig = 0
    Nprog = nlat * nlon
    for iy in range(nlat):
        #print(iy)
        for ix in range(nlon):
            show_progress(ig, Nprog)
            
            plon = LON[iy, ix]
            plat = LAT[iy, ix]
            #if iy == 56:
            #    print(iy, ix)
            #north_dis = DIS_point_to_points(plon, plat, LON[iy + 1:, ix], LAT[iy + 1:, ix])
            #east_dis = DIS_point_to_points(plon, plat, LON[iy, ix + 1:], LAT[iy, ix + 1:])
            lat_dis = DIS_point_to_points(plon, plat, LON[:, ix], LAT[:, ix])
            lon_dis = DIS_point_to_points(plon, plat, LON[iy, :], LAT[iy, :])
            
            #sn_r = np.sum(north_dis <= 2 * Ls)
            #we_r = np.sum(east_dis <= 2 * Ls)
            #n_lim = min(iy + sn_r, nlat)
            #s_lim = max(iy - sn_r, 0)
            #w_lim = max(ix - we_r, 0)
            #e_lim = min(ix + we_r, nlon)
            #print(iy, ix, s_lim, n_lim, w_lim, e_lim)
            y_sel = lat_dis <= LsPar * Ls
            x_sel = lon_dis <= LsPar * Ls
            LON_region = LON[y_sel][:, x_sel]
            LAT_region = LAT[y_sel][:, x_sel]
            #IX_region = IX[s_lim:n_lim, w_lim:e_lim]
            #IY_region = IY[s_lim:n_lim, w_lim:e_lim]
            IG_region = IG[y_sel][:, x_sel]
            LON_region_list = LON_region.flatten()
            LAT_region_list = LAT_region.flatten()
            IG_region_list = IG_region.flatten()
            #IX_region_list = IX_region.flatten()
            #IY_region_list = IY_region.flatten()
            #ind = np.where(IX_region_list + IY_region_list > iy + ix)[0]
            #print(len(ind) / len(LON_region_list))
            
            disList = DIS_point_to_points(plon, plat, LON_region_list[:], LAT_region_list[:])
            valid_ind = disList <= LsPar * Ls
            
            IG_region_list = IG_region_list[valid_ind]
            disList = disList[valid_ind]
            assert np.sum(disList > LsPar * Ls) == 0
            DLen = len(disList)
            assert DLen == len(IG_region_list)
            ncvar_data[Pt_DataDim : Pt_DataDim + DLen] = np.exp( -disList / Ls )
            ncvar_indices[Pt_DataDim : Pt_DataDim + DLen] = IG_region_list
            ncvar_indptr[ig + 1] = Pt_DataDim + DLen
            
            Pt_DataDim += DLen
            ig += 1
    ncf_E.close()




