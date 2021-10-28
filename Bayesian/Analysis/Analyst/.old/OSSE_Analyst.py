#!/usr/bin/env python
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import datetime as dtm
from .BIS_Analyst import BIS_Case_Analyst

class OSSE_Analyst(BIS_Case_Analyst):
    def __init__(self, **kwargs):
        BIS_Case_Analyst.__init__(self, **kwargs)
        self.CaseType = "OSSE"

    def get_truth_emiss(self, time):
        assert False, "Please define the function \"get_truth_emiss\"."

    def Truth_emiss_point_value(self, Time, plon, plat):
        dif_lon = np.abs(plon - self.lon)
        ix = int(np.where(dif_lon == dif_lon.min())[0][0])
        dif_lat = np.abs(plat - self.lat)
        iy = int(np.where(dif_lat == dif_lat.min())[0][0])
        getFun = self.get_truth_emiss
        emiss = getFun(Time)[iy, ix]
        return emiss

    def get_SiteLoc_Truth_emiss(self, start = None, end = None, dt = None):
        if start == None:
            start = self.Start
        if end == None:
            end = self.End
        if dt == None:
            dt = self.dt

        timeList = []
        emiss_dic = {}
        
        for station in self.site_location_dic:
            emiss_dic[station] = []
            

        Current = start
        while(Current <= end):
            timeList.append(Current)
            for station in self.site_location_dic:
                plon = self.site_location_dic[station][0]
                plat = self.site_location_dic[station][1]
                emiss_value = self.Truth_emiss_point_value(Time = Current, plon = plon, plat = plat)
                emiss_dic[station].append(emiss_value)
                
            Current += dtm.timedelta(hours = dt)

        for station in self.site_location_dic:
            emiss_dic[station] = np.array(emiss_dic[station])
            
        timeList = np.array(timeList)
        return timeList, emiss_dic

class OSSE_Analyst_zero_one(OSSE_Analyst):

    def __init__(self, **kwargs):
        OSSE_Analyst.__init__(self, **kwargs)

    """
    def get_truth_emiss(self, time):
        ncf = nc.Dataset(self.CaseDir + "/emiss/Truth.nc", "r")
        ftime = ncf.variables["time"][:].filled(np.nan)
        utime = ncf.variables["time"].units
        numtime = nc.date2num(time, units = utime)
        emiss = ncf.variables["emiss"][ftime == numtime, :, :].filled(np.nan)[0, :, :]
        assert emiss.shape[0] == self.nlat
        assert emiss.shape[1] == self.nlon
        return emiss 
    """
    def get_truth_emiss(self, time):
        ncf = nc.Dataset(self.CaseDir + "/emiss/Truth_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc", "r")
        emiss = ncf.variables["emiss"][:, :].filled(np.nan)
        assert emiss.shape[0] == self.nlat
        assert emiss.shape[1] == self.nlon
        return emiss     
    
class OSSE_Analyst_const_meic(OSSE_Analyst):

    def __init__(self, **kwargs):
        OSSE_Analyst.__init__(self, **kwargs)
     
    def get_truth_emiss(self, time):
        ncf = nc.Dataset(self.CaseDir + "/emiss/Truth_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc", "r")
        emiss = ncf.variables["emiss"][:, :].filled(np.nan)
        assert emiss.shape[0] == self.nlat
        assert emiss.shape[1] == self.nlon
        return emiss     
