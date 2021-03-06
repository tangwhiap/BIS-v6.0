#!/usr/bin/env python

import datetime as dtm
import numpy as np
import netCDF4 as nc
import xarray as xr
from pdb import set_trace
from scipy import interpolate
from ._utils import areaS, resample_points

'''
Athors: Wenhan TANG - 03/2021
'''

hasInterpolated = True
interpolate_method = "linear"

DataDir = "/home/tangwh/datasets/NDRC/Processed"

sectorToType = {
    "power": "point", "cement": "point", "indusp": "point", "other": "point", "heat": "point",
    "resid": "area", "indus": "area", "mobile": "line",
}


def interface(Time, LON, LAT, sector = None, **kwargs):

    Time = Time + dtm.timedelta(hours = 8) # NDRC using Beijing location time, convert UTC time to +8 zone's LST
    time_first = dtm.datetime(2016,1,1,0,0)
    time_last = dtm.datetime(2016,12,31,23,0)
    wrf_weekday = Time.weekday()
    ndrc_time = Time.replace(year=2016) # NDRC data in 2016    
    ndrc_weekday = ndrc_time.weekday()
    daydt = wrf_weekday - ndrc_weekday
    if daydt > 3:
       daydt = daydt - 7
    elif daydt < -3:
       daydt = daydt + 7
    #print("before: ", ndrc_time)
    ndrc_time = ndrc_time + dtm.timedelta(days=daydt)
    #print("after: ", ndrc_time)
    if ndrc_time < time_first:
        ndrc_time = ndrc_time + dtm.timedelta(days=7)
    if ndrc_time > time_last:
        ndrc_time = ndrc_time - dtm.timedelta(days=7)
    #print("after: ", ndrc_time)

    if sector is None or sector == "total":
        emiss_shape = "total"
        emiss_sector = "total"
    else:
        emiss_shape = sectorToType[sector]
        emiss_sector = sector

    if emiss_shape == "point":
        isPoint = True
    else:
        isPoint = False

    ncf = nc.Dataset(DataDir + "/hourly/" + emiss_shape + "/NDRC_" + emiss_sector + "_hourly.nc", "r")

    #ndrc_out = ds.sel(time = ndrc_time)
    #lon = ndrc_out["lon"].values
    #lat = ndrc_out["lat"].values
    ncf_time = ncf.variables["time"][:].filled(np.nan)
    tunit = ncf.variables["time"].units
    tcalendar = ncf.variables["time"].calendar
    tnum = nc.date2num(ndrc_time, units = tunit, calendar = tcalendar)
    tindex = int(np.where(ncf_time == tnum)[0])

    if isPoint:
        pLon = ncf.variables["pLon"][:].filled(np.nan)
        pLat = ncf.variables["pLat"][:].filled(np.nan)
        pData = ncf.variables["pData"][tindex].filled(np.nan)
        emiss_regrided = resample_points(pData, pLon, pLat, LON, LAT)

        # convert tCO2/cell/hour --> mol/km^2/hour
        matrixS = areaS(LON, LAT)
        emiss_regrided = emiss_regrided / matrixS
        mass_CO2 = 44.
        emiss_regrided = emiss_regrided * 1000000. / mass_CO2

    else:
        lon = ncf.variables["lon"][:].filled(np.nan)
        lat = ncf.variables["lat"][:].filled(np.nan)
        LON_orig, LAT_orig = np.meshgrid(lon, lat)
        emiss = ncf.variables["emiss"][tindex].filled(np.nan)

        #set_trace()
        # convert tCO2/cell/hour --> mol/km^2/hour
        # convert tCO2/cell/hour --> umol/m^2/s
        matrixS = areaS(LON_orig, LAT_orig)
        #dx_long = np.gradient(LON)[1] * np.cos(np.deg2rad(LAT)) * 111.320
        #dy_long = np.gradient(LAT)[0] * 110.574
        #area_div = np.absolute(dx_long * dy_long) 

        emiss = emiss / (matrixS * 1000000.)
        mass_CO2 = 44.
        emiss = emiss * 1000000. / mass_CO2 * 1000000. / 3600.

        emiss_regrided = interpolate.griddata((LON_orig.flatten(), LAT_orig.flatten()), emiss.flatten(), (LON, LAT), method = "linear")

        

        #return {"lon":LON.flatten(), "lat":LAT.flatten(), "value":emiss_regrided.flatten()}

    return emiss_regrided, LON, LAT
