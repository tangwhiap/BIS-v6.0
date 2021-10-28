#!/usr/bin/env python

#-> Picarro sites observation <-#

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from ..base.main_class import OBS_base
from .GetStationInfo import getinfo
from .getdata import getdata

import numpy as np
import xarray as xr
import datetime as dtm
from pdb import set_trace

class OBS_picarro(OBS_base):

    makeType = "offline"
    def __init__(self, directory, sites, location, errors):

        self.picarroDir = directory
        self.picarroSites = sites
        self.nSite = len(self.picarroSites)
        self.sitesError = errors

        assert self.makeType.lower() in ["offline", "online"]
        isOffline = True if self.makeType.lower() == "offline" else False
        isOnline = True if self.makeType.lower() == "online" else False
        infoDir = "/".join(__file__.split("/")[:-1])
        self.stationsInfo = getinfo(infoDir)

        self.obsName_dic = {siteName: index for index, siteName in enumerate(self.picarroSites)}
        self.obsName = self.picarroSites
        self.obsLoc = location
        

    def get_data_offline(self, start, end, dt):

        if isinstance(start, dtm.datetime):
            start = start.strftime("%Y-%m-%d_%H:%M:%S")
            end = end.strftime("%Y-%m-%d_%H:%M:%S")

        if isinstance(dt, dtm.timedelta):
            dt = dt.total_seconds() / 3600.0

        self.timeList = None
        self.obsDic = {}
        for station in self.picarroSites:

            assert station in self.picarroSites
            dirName, staLon, staLat, staAgl = self.stationsInfo[station]
            inputDir = self.picarroDir + "/" + dirName + "/K30"
            timeList, obs = getdata(start, end, dt, inputDir)

            self.obsDic[station] = xr.Dataset({"co2": (["time"], obs)}, coords = {"time": (["time"], timeList)})
            if self.timeList is None:
                self.timeList = timeList
            else:
                assert timeList == self.timeList

    def get_obs_orig(self, time, sites = None):
        if sites is None:
            sites = self.picarroSites
        if isinstance(sites, str):
            sites = [sites]

        #ind = self.timeList.index(time)
        obsList = []
        for site in sites:
            ds = self.obsDic[site]
            obsData = ds.sel(time = np.datetime64(time))["co2"].values
            obsList.append(obsData)
        obsList = np.array(obsList)#.reshape(self.nSite, 1)
        return obsList

    def get_obs_proc(self, time, sites = None):
        if sites is None:
            sites = self.picarroSites
        if isinstance(sites, str):
            sites = [sites]

        obsList = []
        for site in sites:
            ds = self.obsDic[site]
            timeList = ds.time.values
            timeStart = time
            timeEnd = time + dtm.timedelta(hours = 1)
            time_range = timeList[ ( timeList >= np.datetime64(timeStart)) & ( timeList <= np.datetime64(timeEnd)) ] 
            co2obs_list = ds.sel(time = time_range)["co2"].values
            co2obs_list = self.filter_min_max(co2obs_list, vmin = 350, vmax = 800)
            co2obs = np.nanmean(co2obs_list)
            obsList.append(float(co2obs))

        obsList = np.array(obsList)

        return obsList



    def get_errors(self, time):
        return self.sitesError







