#!/usr/bin/env python

#-> API code of STILT output (STILT_forBIS_v*) <-#

# Authors:
#   Wenhan TANG - 07/2021
#   ...

import numpy as np
import netCDF4 as nc
import glob
import datetime as dtm

#footDir = "/home/tangwh/modeldata/STILT_cases/BTH_d02_bt48/stilt_out"
from pdb import set_trace

class stilt(object):

    def __init__(self, receptors, footDir):
        self.footDir = footDir
        if isinstance(receptors, dict):
            self.receptors = receptors
        elif isinstance(receptors, list):
            self.receptors = {}
            for receptor in receptors:
                self.receptors[receptor] = self.get_site_location(receptor, outType = "string")
        else:
            assert False, "TypeError of receptors"

        self.recepTimeList = self.get_site_timeList()
        self.receptorsName = list(self.receptors.keys())

    def get_footprint(self, *args, **kwargs):
        return stilt_forBIS_interface(self, *args, **kwargs)

    def get_dom_config(self, *args, **kwargs):
        return get_dom_config(self, *args, **kwargs)

    def get_site_location(self, receptor, outType = "string"):
        fileList = glob.glob(self.footDir + "/" + receptor + "/by-id/*")
        assert len(fileList) >= 1
        fileSample = fileList[0]
        temp = fileSample.strip().split("/")[-1]
        temp = temp.strip().split("_")
        strLon = temp[1]
        strLat = temp[2]
        strAgl = temp[3]
        if outType.lower() == "string":
            return {"lon": strLon, "lat": strLat, "agl": strAgl}
        else:
            return {"lon": float(strLon), "lat": float(strLat), "agl": float(strAgl)}

    def get_site_timeList(self, receptorName = None):

        if receptorName is None:
            receptorName = list(self.receptors.keys())[0]

        fileList = glob.glob(self.footDir + "/" + receptorName + "/by-id/*")
        assert len(fileList) >= 1
        timeList = []
        for ifile in fileList:
            temp = ifile.strip().split("/")[-1]
            temp = temp.strip().split("_")
            strTime = temp[0]
            timeList.append(dtm.datetime.strptime(strTime, "%Y%m%d%H%M"))

        return timeList
        

def stilt_forBIS_interface(self, receptorName, time, recepTime):

    location = self.receptors[receptorName]
    strLon = location["lon"]
    strLat = location["lat"]
    strAgl = location["agl"]

    footDirFileName = self.footDir + "/" + receptorName + "/by-id/" +  recepTime.strftime("%Y%m%d%H%M") + "_" + strLon + "_" + strLat + "_" + strAgl + "/" + recepTime.strftime("%Y%m%d%H%M") + "_" + strLon + "_" + strLat + "_" + strAgl + "_foot.nc"

    ncf_Foot = nc.Dataset(footDirFileName, "r")
    nLon = ncf_Foot.dimensions["lon"].size
    nLat = ncf_Foot.dimensions["lat"].size
    lonList = ncf_Foot.variables["lon"][:].filled(np.nan)
    latList = ncf_Foot.variables["lat"][:].filled(np.nan)
    tunits = ncf_Foot.variables["time"].units
    tcalendar = ncf_Foot.variables["time"].calendar
    tv = nc.date2num(time, units = tunits, calendar = tcalendar)
    ftv = ncf_Foot.variables["time"][:].filled(np.nan)
    tid = np.where(ftv == tv)[0]
    assert len(tid) == 0 or len(tid) == 1
    if len(tid) == 0:
        footPrint = np.full( (nLat, nLon), 0.0 )
    else:
        footPrint = ncf_Foot.variables["foot"][int(tid)].filled(np.nan)

    assert footPrint.shape == (nLat, nLon)
    ncf_Foot.close()
    return footPrint, lonList, latList
        
def get_dom_config(self, receptorName = None):

    if receptorName is None:
        receptorName = list(self.receptors.keys())[0]
    footFiles = glob.glob(self.footDir + "/" + receptorName + "/by-id/*")
    footFiles = glob.glob(footFiles[0] + "/*_foot.nc")
    assert len(footFiles) >= 1, "Couln'd found " + self.footDir + "/" + receptorName + "/by-id/*/*_foot.nc"
    footFileSample = footFiles[0]
    ncf_footSample = nc.Dataset(footFileSample)
    lon = ncf_footSample.variables["lon"][:].filled(np.nan)
    lat = ncf_footSample.variables["lat"][:].filled(np.nan)
    LON, LAT = np.meshgrid(lon, lat)
    ncf_footSample.close()
    return LON, LAT



