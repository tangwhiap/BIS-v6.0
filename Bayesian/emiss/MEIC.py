#!/usr/bin/env python
import numpy as np
import pandas as pd
import datetime as dtm
import netCDF4 as nc

hasInterpolated = False
interpolate_method = "linear"

MEIC_Dir = "/home/tangwh/datasets/MEIC/Rescaled"

def interface(time, lon_s = -180, lon_e = 180, lat_s = -90, lat_e = 90, **kwargs):

    sectorType = kwargs["sector"] 
    ds = nc.Dataset(MEIC_Dir + "/MEIC_2016_CO2_hourly.nc", "r")
    lonlist_file = ds.variables["lon"][:].filled(np.nan)
    latlist_file = ds.variables["lat"][:].filled(np.nan)
    lonlist_ind = (lonlist_file >= lon_s - 1) & (lonlist_file <= lon_e + 1)
    latlist_ind = (latlist_file >= lat_s - 1) & (latlist_file <= lat_e + 1)

    lonlist = lonlist_file[lonlist_ind]
    latlist = latlist_file[latlist_ind]
    timelist = ds.variables["time"][:].filled(np.nan)
    nlon = len(lonlist)
    nlat = len(latlist)

    if type(time) is type(pd.Timestamp(0)):
        this_time = time.to_pydatetime()
    else:
        this_time = time
    this_weekday = this_time.weekday()
    this_time = this_time.replace(year=2016)
    time_first = dtm.datetime(2016,1,1,0,0)
    time_last = dtm.datetime(2016,12,31,23,0)
    meic_weekday = this_time.weekday()
    dtday = this_weekday - meic_weekday
    if dtday > 3:
        dtday = dtday - 7
    elif dtday < -3:
        dtday = dtday + 7
    this_time = this_time + dtm.timedelta( days = dtday )
    if this_time > time_last:
        this_time = this_time - dtm.timedelta( days = 7 )
    if this_time < time_first:
        this_time = this_time + dtm.timedelta( days = 7 )
#emiss_out = ds.emiss.sel(time = np.datetime64(this_time), lon = lonlist, lat = latlist).values
    emiss_out = ds.variables[sectorType][timelist == nc.date2num(this_time, units = ds.variables["time"].units, calendar = ds.variables["time"].calendar), latlist_ind, lonlist_ind].filled(np.nan)[0]
#set_trace()
    emiss_out = np.where(np.isnan(emiss_out), 0, emiss_out)
    # convert tCO2/cell/hour --> umol/m^2/s
    LON, LAT = np.meshgrid(lonlist, latlist)
    dx_long = np.gradient(LON)[1] * np.cos(np.deg2rad(LAT)) * 111.320
    dy_long = np.gradient(LAT)[0] * 110.574
    area_div = np.absolute(dx_long * dy_long)
    # gCO2/m^2/hour
    emiss_out = emiss_out * 1000000 / area_div / 1000000
    # umol/m^2/hour
    emiss_out = emiss_out / 44 * 1000000
    # umol/m^2/s
    emiss_out = emiss_out / 3600
    #LON_domain, LAT_domain = np.meshgrid(cfg.lon, cfg.lat)
    #emiss_out = interpolate.griddata( (LON.flatten(), LAT.flatten()), emiss_out.flatten(), (LON_domain, LAT_domain))
    ds.close()
    return emiss_out, LON, LAT
