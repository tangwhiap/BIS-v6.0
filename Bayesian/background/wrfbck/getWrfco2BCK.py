#!/usr/bin/env python
# Compute Hx for WRF-CO2-letkf
# Authors:
#   Wenhan TANG - 05/2021
#   ...

"""
str_Time = "201901020000" # YYYYmmddHHMM

interval = 10 # Minuts

DomID = 1 # DonID = 1 for domain: d01

Nensemble = 20 # Number of ensemble

InDir = "/home/tangwh/modeling/WDAS_v1.1/WRF-CO2-v3.0-DA/model/ASSI/output/wrfco2" # Directory of wrfco2*.nc

OutDir = "/home/tangwh/modeling/WDAS_v1.1/Hx/output" # Direcotory of output files.

ObsDir = "/home/tangwh/modeling/WDAS_v1.1/Hx/obs" # Directory of obs files.

ObsFileName = "Picarro_CO2_2019-01-01_2019-01-02.nc"

ObsInput_dict = {"pr003CMAH1": (1, 1), "pr003CMAH2": (2, 1), "pr004CMAH1": (3, 1), "pr004CMAH2": (4, 1), "pr005CMAH1": (5, 1), "pr005CMAH2": (6, 1), "pr006CMAH1": (7, 1), "pr006CMAH2": (8, 1), "pr007CMAH1": (9, 1), "pr007CMAH2": (10, 1), "pr008CMAH1": (11, 1), "pr008CMAH2": (12, 1)} #{ #station name# : (OBS_id, obs_error) }

undef = -9999.0

BinaryOut = True # True: binary output; False: ASCII output.
"""
#############################################
#  Code for global variables

import numpy as np
import datetime as dtm
import netCDF4 as nc
from scipy.interpolate import interp1d
from cftime import date2num
from glob import glob
from pdb import set_trace

#Wrfco2_Prefix = "wrfout"

"""
Time = dtm.datetime.strptime(str_Time, "%Y%m%d%H%M")

ncf_obs = nc.Dataset(ObsDir + "/" + ObsFileName, "r")
obs_timeList = ncf_obs.variables["time"][:].filled(np.nan)
obs_timeUnits = ncf_obs.variables["time"].units
obs_timeCalendar = ncf_obs.variables["time"].calendar
obs_StaName = ncf_obs.variables["station"][:]
obs_lon = ncf_obs.variables["lon"][:].filled(np.nan)
obs_lat = ncf_obs.variables["lat"][:].filled(np.nan)
obs_agl = ncf_obs.variables["AGL"][:].filled(np.nan)

ncf_wrfco2 = nc.Dataset(InDir + "/wrfco2_d" + str(DomID).zfill(2) + "_" + Time.strftime("%Y-%m-%d_%H:%M:%S"), "r")

"""


#SampleFile = glob(InDir + "/wrfco2_d" + str(DomID).zfill(2) + "_*")[0]

############################################
# Definition of useful functions

"""
def getOBS(ObsName): 
    global Time, interval, ncf_obs, obs_timeList, obs_timeUnits, obs_tineCalendar, obs_StaName, undef
    assert ObsName in obs_StaName

    # Average between "Time - 5" to "Time + 4" (if interval == 10)
    ave_start = Time - dtm.timedelta(minutes = int(interval / 2))
    ave_end = ave_start + dtm.timedelta(minutes = interval)

    ave_start_int = date2num(ave_start, units = obs_timeUnits, calendar = obs_timeCalendar)
    ave_end_int = date2num(ave_end, units = obs_timeUnits, calendar = obs_timeCalendar)
    time_range = ((obs_timeList >= ave_start_int) & (obs_timeList <= ave_end_int))
    ObsCO2 = ncf_obs.variables["CO2"][ObsName == obs_StaName, time_range][0].filled(np.nan)
    #set_trace()
    ObsCO2_mean = np.nanmean(ObsCO2)
    return undef if np.isnan(ObsCO2_mean) else ObsCO2_mean
""" 
    

def GetDomInfo(SampleFile):
    #global nx, ny, nz
    sf = SampleFile
    lonlist = sf.variables["XLONG"][0,0].filled(np.nan)
    latlist = sf.variables["XLAT"][0,:,0].filled(np.nan)
    nx = sf.dimensions["west_east"].size
    ny = sf.dimensions["south_north"].size
    nz = sf.dimensions["bottom_top"].size
    assert nx == len(lonlist)
    assert ny == len(latlist)
    lon_s = lonlist[0]
    lat_s = latlist[0]
    lon_e = lonlist[-1]
    lat_e = latlist[-1]
    dlon = (lonlist[1:] - lonlist[:-1]).mean()
    dlat = (latlist[1:] - latlist[:-1]).mean()
    return {"lon_s": lon_s, "lon_e": lon_e, "lat_s": lat_s, "lat_e": lat_e, "dlon": dlon, "dlat": dlat}

def anchor(LocDic, wrfdom):
    #global ObsInput_dict
    LocDic_new = {}
    def isIn(lon,lat):
        if lon >= wrfdom["lon_s"] and lon <= wrfdom["lon_e"] and lat >= wrfdom["lat_s"] and lat <= wrfdom["lat_e"]:
            return True
        else:
            return False

    for location in LocDic:
        lon = LocDic[location]["lon"]
        lat = LocDic[location]["lat"]
        if not(isIn(lon, lat)):# or not(location in ObsInput_dict):
            continue
        
        wsx = int((lon - wrfdom["lon_s"]) / wrfdom["dlon"])
        wslon = wsx * wrfdom["dlon"] + wrfdom["lon_s"]
        wsy = int((lat - wrfdom["lat_s"]) / wrfdom["dlat"])
        wslat = wsy * wrfdom["dlat"] + wrfdom["lat_s"]
        wel = lon - wslon
        wer = wrfdom["dlon"]- wel
        snb = lat - wslat
        snt = wrfdom["dlat"] - snb
        assert wel >= 0 and wer >=0 and snb >=0 and snt >=0
        LocDic_new[location] = LocDic[location].copy()
        #LocDic_new[location]["ID"] = ObsInput_dict[location][0]
        #LocDic_new[location]["error"] = ObsInput_dict[location][1]
        LocDic_new[location]["wsx"] = wsx
        LocDic_new[location]["wsy"] = wsy
        LocDic_new[location]["wel"] = wel
        LocDic_new[location]["wer"] = wer
        LocDic_new[location]["snb"] = snb
        LocDic_new[location]["snt"] = snt
    return LocDic_new

def belinear(array, location):
    wsx = location["wsx"]
    wsy = location["wsy"]
    wel = location["wel"]
    wer = location["wer"]
    snb = location["snb"]
    snt = location["snt"]

    ws = array[wsy, wsx]
    es = array[wsy, wsx + 1]
    wn = array[wsy + 1, wsx]
    en = array[wsy + 1, wsx + 1]
    return (ws * (wer * snt) + es * (wel * snt) + en * (wel * snb) + wn * (wer * snb)) / ((wel + wer) * (snb + snt))

def spatial_interpolate_3D(array, location):
    # dim: (z, lat, lon)
    array = np.swapaxes(array, 0, 1) #dim: (lat, z, lon)
    array = np.swapaxes(array, 1, 2) #dim: (lat, lon, z)
    profile = belinear(array, location)
    return profile

def spatial_interpolate_2D(array, location):
    return belinear(array, location)

def profile_interpolate_1D(profile_bef, coord1, coord2):
    assert len(profile_bef.shape) == len(coord1.shape), "Conflict dimensions of variable and coordinate"
    assert len(profile_bef) == len(coord1), "Conflict length of variable and coordinate"
    f=interp1d(coord1, profile_bef, bounds_error = False, kind = "quadratic")
    return f(coord2)

def make_LocDic(Receptors, Wrfco2Dir, Wrfco2_Prefix, Domid):
    fileList = glob(Wrfco2Dir + "/" + Wrfco2_Prefix + "_d" + str(Domid).zfill(2) + "_*")
    ncf_wrfco2_sample = nc.Dataset(fileList[0], "r")
    wrfdom = GetDomInfo(ncf_wrfco2_sample)
    ncf_wrfco2_sample.close()
    LocDic = {}
    for name in Receptors:
        station_info = Receptors[name]
        lon = float(station_info[0])
        lat = float(station_info[1])
        agl = float(station_info[2])
        LocDic[name] = {"lon": lon, "lat": lat, "agl": agl}
    
    LocDic = anchor(LocDic, wrfdom)
    return LocDic

def getBCK(Time, LocDic, Wrfco2Dir, Wrfco2_Prefix, Domid):
    ncf_wrfco2 = nc.Dataset(Wrfco2Dir + "/" + Wrfco2_Prefix + "_d" + str(Domid).zfill(2) + "_" + Time.strftime("%Y-%m-%d_%H:%M:%S"), "r")
    Z_3D = ( ncf_wrfco2.variables["PH"][0].filled(np.nan) + ncf_wrfco2.variables["PHB"][0].filled(np.nan) ) / 9.8
    Z_3D = (Z_3D[1:] + Z_3D[:-1]) / 2
    #P_3D = ncf_wrfco2.variables["P"][0].filled(np.nan) + ncf_wrfco2.variables["PB"][0].filled(np.nan)
    #P_3D = P_3D / 100 # Convert unit from Pa to hPa
    #P_SFC = ncf_wrfco2.variables["PSFC"][0].filled(np.nan)
    #P_SFC = P_SFC / 100 # Convert unit from Pa to hPa
    HGT_2D = ncf_wrfco2.variables["HGT"][0].filled(np.nan)
    co2_bck = ncf_wrfco2.variables["CO2_BCK"][0, :, :, :].filled(np.nan)
    ncf_wrfco2.close()
    BCK_dict = {}
    for station in LocDic:
        Z_eta_layers = spatial_interpolate_3D(Z_3D, LocDic[station]) - spatial_interpolate_2D(HGT_2D, LocDic[station])
        co2_profile = spatial_interpolate_3D(co2_bck, LocDic[station])
        co2_station = profile_interpolate_1D(np.array([co2_profile[0]] + co2_profile.tolist()), np.array([0] + Z_eta_layers.tolist()), LocDic[station]["agl"])
        #P_profile = spatial_interpolate_3D(P_3D, LocDic[station])
        #P_surface = spatial_interpolate_2D(P_SFC, LocDic[station])
        #set_trace()
        #P_station = profile_interpolate_1D(np.array([P_surface] + P_profile.tolist()), np.array([0] + Z_eta_layers.tolist()), LocDic[station]["agl"])
        BCK_dict[station] = float(co2_station)
    return BCK_dict

if __name__ == "__main__":
    wrfdom = GetDomInfo(ncf_wrfco2)
    LocDic = anchor(LocDic, wrfdom)
    #set_trace()
    Hx_dict = Hx(LocDic, iMember)
    
#ncf_obs.close()
#ncf_wrfco2.close()
        
