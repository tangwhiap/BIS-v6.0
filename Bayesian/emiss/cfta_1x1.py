#!/usr/bin/env python

'''
Authors: Wenhan TANG - 12/2020
'''

import datetime as dtm
from datetime import datetime
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
from pdb import set_trace

hasInterpolated = False
interpolate_method = "linear"

DataDir = "/home/tangwh/datasets/VEGAS_Cfta"

def interface(time, lon_s = -180, lon_e = 180, lat_s = -90, lat_e = 90, **kwargs):

    #DataDir = DataDir.strip().split(":")
	
    infile = DataDir + "/VEGAS_CFta_hr_Jan2000_Oct2020_1.0.nc"
    datain = Dataset(infile,'r')
    gtimes = datain.variables['time']
    glons = datain.variables['lon'][:]
    glons = ((glons+180)%360)-180
    glats = datain.variables['lat'][:]

    tindex = nc.date2index(time,gtimes)

    valid_lon = (glons >= lon_s - 1) & (glons <= lon_e + 1)
    valid_lat = (glats >= lat_s - 1) & (glats <= lat_e + 1)
    gcfta = datain.variables['cfta'][tindex, :, :].filled(np.nan)
    cfta = gcfta[valid_lat][:, valid_lon]

    #To march the dims of gcfta with glons and glats
    #gcfta = gcfta.reshape(1,gcfta.shape[0],gcfta.shape[1]) 

    ## convert CFta units from KgC/m2/sec to mol/km2/hr for WRF
    #gcfta = gcfta * 1000./12.*1000000.*3600. ## 1000g/kg, 12g/mol, 1e6 m2/km2, 3600 sec/hr

    ## convert CFta units from KgC/m2/sec to umol/m2/sec
    cfta = cfta * 1000.0 / 12.0 * 1000000.0

    #glons,glats=np.meshgrid(glons,glats)
    #return {"lon" : glons.flatten(), "lat" : glats.flatten(), "value" : gcfta.flatten()}
    lon = glons[valid_lon]
    lat = glats[valid_lat]
    LON, LAT = np.meshgrid(lon, lat)

    return cfta, LON, LAT
