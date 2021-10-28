#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from ...utils.netcdf_io import get_ncvar, get_nctime, get_ncdim
from ...utils.sparse_matrix import getNCSPR_spr, getNCSPR_array
from ...utils.distance import areaS
from ...utils.mask_tools import Maskout

from .local_point import LocPoint

import numpy as np
import netCDF4 as nc
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import datetime as dtm
from pdb import set_trace
import glob
import os


class Analyst(object):
    
    def __init__(self, objCase):

        self.objCase = objCase
        self.caseName = objCase.CaseConfig["CASENAME"]
        self.caseConf = objCase.CaseConfig

        self.caseStart = objCase.objExp.Start
        self.caseEnd = objCase.objExp.End
        self.caseDt = objCase.objExp.dt
        self.caseDtHrs = objCase.objExp.dtHrs

        self.isOSSE = objCase.objExp.isOSSE
        self.isREAL = objCase.objExp.isREAL

        self.iterTimeList = objCase.objExp.objCore.objIter.iterTimeList
        self.obsTimeList = objCase.objExp.objCore.objIter.obsTimeList

        self.sectors = objCase.objExp.sectors
        self.sampleSector = self.sectors[0]

        

        self.get_priorFile_name = objCase.objExp.get_priorFile_name
        self.get_postFile_name = objCase.objExp.get_postFile_name
        self.get_sigmaPriorFile_name = objCase.objExp.get_sigmaPriorFile_name
        self.get_sigmaProcFile_name = objCase.objExp.get_sigmaProcFile_name
        self.get_obsFile_name = objCase.objExp.get_obsFile_name
        self.get_bckPriorFile_name = objCase.objExp.get_bckPriorFile_name
        self.get_bckProcFile_name = objCase.objExp.get_bckProcFile_name
        self.get_sbckPriorFile_name = objCase.objExp.get_sbckPriorFile_name
        self.get_sbckProcFile_name = objCase.objExp.get_sbckProcFile_name
        self.get_initHxFile_name = objCase.objExp.get_initHxFile_name
        self.get_initHxBckFile_name = objCase.objExp.get_initHxBckFile_name
        self.get_finalHxFile_name = objCase.objExp.get_finalHxFile_name
        self.get_finalHxBckFile_name = objCase.objExp.get_finalHxBckFile_name
        self.get_bckPriorFile_name = objCase.objExp.get_bckPriorFile_name
        self.get_bckProcFile_name = objCase.objExp.get_bckProcFile_name
        self.get_sbckPriorFile_name = objCase.objExp.get_sbckPriorFile_name
        self.get_sbckProcFile_name = objCase.objExp.get_sbckProcFile_name

        

        self.get_dsigma2File_name = objCase.objExp.objCore.get_dsigma2File_name

        self.get_priorHxFile_name = objCase.objExp.objCore.get_priorHxFile_name
        self.get_priorHxBckFile_name = objCase.objExp.objCore.get_priorHxBckFile_name
        self.get_procHxFile_name = objCase.objExp.objCore.get_procHxFile_name
        self.get_procHxBckFile_name = objCase.objExp.objCore.get_procHxBckFile_name

        self.get_postHourlyFile_name = objCase.objExp.get_postHourlyFile_name
        self.get_spostHourlyFile_name = objCase.objExp.get_spostHourlyFile_name
        self.get_postDailyFile_name = objCase.objExp.get_postDailyFile_name
        self.get_spostDailyFile_name = objCase.objExp.get_spostDailyFile_name
        self.get_postWeeklyFile_name = objCase.objExp.get_postWeeklyFile_name
        self.get_spostWeeklyFile_name = objCase.objExp.get_spostWeeklyFile_name
        self.get_postMonthlyFile_name = objCase.objExp.get_postMonthlyFile_name
        self.get_spostMonthlyFile_name = objCase.objExp.get_spostMonthlyFile_name
        self.get_postAllFile_name = objCase.objExp.get_postAllFile_name
        self.get_spostAllFile_name = objCase.objExp.get_spostAllFile_name

        self.get_priorHourlyFile_name = objCase.objExp.get_priorHourlyFile_name
        self.get_spriorHourlyFile_name = objCase.objExp.get_spriorHourlyFile_name
        self.get_priorDailyFile_name = objCase.objExp.get_priorDailyFile_name
        self.get_spriorDailyFile_name = objCase.objExp.get_spriorDailyFile_name
        self.get_priorWeeklyFile_name = objCase.objExp.get_priorWeeklyFile_name
        self.get_spriorWeeklyFile_name = objCase.objExp.get_spriorWeeklyFile_name
        self.get_priorMonthlyFile_name = objCase.objExp.get_priorMonthlyFile_name
        self.get_spriorMonthlyFile_name = objCase.objExp.get_spriorMonthlyFile_name
        self.get_priorAllFile_name = objCase.objExp.get_priorAllFile_name
        self.get_spriorAllFile_name = objCase.objExp.get_spriorAllFile_name


        if self.objCase.objExp.isOSSE:
            self.get_truthFile_name = objCase.objExp.get_truthFile_name
            self.get_truthHourlyFile_name = objCase.objExp.get_truthHourlyFile_name
            self.get_truthDailyFile_name = objCase.objExp.get_truthDailyFile_name
            self.get_truthWeeklyFile_name = objCase.objExp.get_truthWeeklyFile_name
            self.get_truthMonthlyFile_name = objCase.objExp.get_truthMonthlyFile_name
            self.get_truthAllFile_name = objCase.objExp.get_truthAllFile_name

        self.hasBCK = objCase.objExp.hasBCK
        self.get_Domain_info()
        self.get_OBS_info()
        self.create_OBS_anchor()
        self.get_output_timeList()

        self.regionsMaskout = {}
    
        print("Successfully create an analyst object of case: " + self.caseName + ".")


    def get_Domain_info(self):

        self.LON_dic = {}
        self.LAT_dic = {}
        self.nDim1_dic = {}
        self.nDim2_dic = {}
        self.areaS_dic = {}

        for sector in self.sectors: 

            SampleFile = self.get_priorFile_name(time = self.caseStart, sectorName = sector)
            LON = get_ncvar(SampleFile, "LON")
            LAT = get_ncvar(SampleFile, "LAT")

            assert LON.shape == LAT.shape
            dim1, dim2 = LON.shape
            #self.nLat, self.nLon = LON.shape
            #self.nGrid = self.nLat * self.nLon
            self.LON_dic[sector] = LON
            self.LAT_dic[sector] = LAT
            self.nDim1_dic[sector] = dim1
            self.nDim2_dic[sector] = dim2
            self.areaS_dic[sector] = areaS(LON, LAT)

        self.sampleLON = self.LON_dic[self.sampleSector]
        self.sampleLAT = self.LAT_dic[self.sampleSector]
        self.sampleAreaS = self.areaS_dic[self.sampleSector]


    def get_OBS_info(self):

        if self.isOSSE:
            objFootDic = self.objCase.objExp.objFootDic
            self.sampleFootName = list(objFootDic.keys())[0]
            self.obsName = self.objCase.objExp.objFootDic[self.sampleFootName].receptorsName
            obsDic = self.objCase.objExp.objFootDic[self.sampleFootName].receptors
            obsLoc = {}
            for site in obsDic:
                obsLoc[site] = {"lon": float(obsDic[site]["lon"]), "lat": float(obsDic[site]["lat"]), "agl": float(obsDic[site]["agl"])}
        if self.isREAL:
            self.obsName = self.objCase.objExp.objOBS.obsName
            obsDic = self.objCase.objExp.objOBS.obsLoc
            obsLoc = {}
            for site in obsDic:
                obsLoc[site] = {"lon": float(obsDic[site][0]), "lat": float(obsDic[site][1]), "agl": float(obsDic[site][2])}

        self.obsLoc = obsLoc

    def get_output_timeList(self):
        
        self.timeList_dic = {}
        self.timeList_dic["obs"] = self.obsTimeList

        expConfig = self.objCase.objExp.myConfig
        time_format = "%Y-%m-%d_%H:%M:%S"
        after_time_format = ".nc"
        
        #elf.myConfig["emissPostDir"] + "/" + self.myConfig["emissPost_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"
        before_time_format = expConfig["emissPostDir"] + "/" + expConfig["emissPost_Prefix"] + "_" + self.sampleSector + "_"
        self.timeList_dic["orig"] = self.read_timeList(before_time_format, time_format, after_time_format)

        before_time_format = expConfig["outHourlyDir"] + "/" + expConfig["postHourly_Prefix"] + "_" + self.sampleSector + "_"
        self.timeList_dic["hourly"] = self.read_timeList(before_time_format, time_format, after_time_format)

        before_time_format = expConfig["outDailyDir"] + "/" + expConfig["postDaily_Prefix"] + "_" + self.sampleSector + "_"
        self.timeList_dic["daily"] = self.read_timeList(before_time_format, time_format, after_time_format)

        before_time_format = expConfig["outWeeklyDir"] + "/" + expConfig["postWeekly_Prefix"] + "_" + self.sampleSector + "_"
        self.timeList_dic["weekly"] = self.read_timeList(before_time_format, time_format, after_time_format)

        before_time_format = expConfig["outMonthlyDir"] + "/" + expConfig["postMonthly_Prefix"] + "_" + self.sampleSector + "_"
        self.timeList_dic["monthly"] = self.read_timeList(before_time_format, time_format, after_time_format)

    def read_timeList(self, before_time_format, time_format, after_time_format):
        fileList = glob.glob(before_time_format + "*" + after_time_format)
        fileList.sort()
        timeList = []
        for file_ in fileList:
            time = dtm.datetime.strptime(file_, before_time_format + time_format + after_time_format)
            timeList.append(time)
        return timeList

    def create_OBS_anchor(self):
        obsPoints = {}
        for obsName in self.obsLoc:
            obsPoints[obsName] = LocPoint(staName = obsName, staLon = self.obsLoc[obsName]["lon"], staLat = self.obsLoc[obsName]["lat"])
        self.obsPoints = obsPoints

    def build_region_maskout(self, regionName, shpName = None, LON = None, LAT = None):
        if regionName not in self.regionsMaskout:
            assert shpName is not None
            if LON is None:
                LON = self.sampleLON
            if LAT is None:
                LAT = self.sampleLAT
            self.regionsMaskout[regionName] = Maskout(LON = LON, LAT = LAT, Region_list = shpName)
            self.regionsMaskout[regionName].mask_array_out()
        return self.regionsMaskout[regionName]

        

    def get_prior_emiss(self, time, sectorName = None, timeScale = "orig"):

        if sectorName is None:
            sectorName = self.sampleSector

        if timeScale.lower() == "orig":
            fileName = self.get_priorFile_name(time, sectorName)

        elif timeScale.lower() == "hourly":
            fileName = self.get_priorHourlyFile_name(time, sectorName)

        elif timeScale.lower() == "daily":
            fileName = self.get_priorDailyFile_name(time, sectorName)

        elif timeScale.lower() == "weekly":
            fileName = self.get_priorWeeklyFile_name(time, sectorName)

        elif timeScale.lower() == "monthly":
            fileName = self.get_priorMonthlyFile_name(time, sectorName)

        elif timeScale.lower() == "all":
            fileName = self.get_priorAllFile_name(time, sectorName)

        else:
            assert False

        emiss, LON, LAT = get_ncvar(fileName, ["emiss", "LON", "LAT"])

        return emiss, LON, LAT


    def get_posterior_emiss(self, time, sectorName = None, timeScale = "orig"):

        if sectorName is None:
            sectorName = self.sampleSector

        if timeScale.lower() == "orig":
            fileName = self.get_postFile_name(time, sectorName)

        elif timeScale.lower() == "hourly":
            fileName = self.get_postHourlyFile_name(time, sectorName)

        elif timeScale.lower() == "daily":
            fileName = self.get_postDailyFile_name(time, sectorName)

        elif timeScale.lower() == "weekly":
            fileName = self.get_postWeeklyFile_name(time, sectorName)

        elif timeScale.lower() == "monthly":
            fileName = self.get_postMonthlyFile_name(time, sectorName)

        elif timeScale.lower() == "all":
            fileName = self.get_postAllFile_name(time, sectorName)

        else:
            assert False


        emiss, LON, LAT = get_ncvar(fileName, ["emiss", "LON", "LAT"])

        return emiss, LON, LAT


    def get_prior_sigma(self, time, sectorName = None, timeScale = "orig"):

        if sectorName is None:
            sectorName = self.sampleSector

        if timeScale.lower() == "orig":
            fileName = self.get_sigmaPriorFile_name(time, sectorName)

        elif timeScale.lower() == "hourly":
            fileName = self.get_spriorHourlyFile_name(time, sectorName)

        elif timeScale.lower() == "daily":
            fileName = self.get_spriorDailyFile_name(time, sectorName)

        elif timeScale.lower() == "weekly":
            fileName = self.get_spriorWeeklyFile_name(time, sectorName)

        elif timeScale.lower() == "monthly":
            fileName = self.get_spriorMonthlyFile_name(time, sectorName)

        elif timeScale.lower() == "all":
            fileName = self.get_spriorAllFile_name(time, sectorName)

        else:
            assert False

        emiss, LON, LAT = get_ncvar(fileName, ["emiss", "LON", "LAT"])

        return emiss, LON, LAT
        
    def get_posterior_sigma(self, time, sectorName = None, timeScale = "orig"):

        if sectorName is None:
            sectorName = self.sampleSector

        if timeScale.lower() == "orig":
            fileName = self.get_sigmaProcFile_name(time, sectorName)

        elif timeScale.lower() == "hourly":
            fileName = self.get_spostHourlyFile_name(time, sectorName)

        elif timeScale.lower() == "daily":
            fileName = self.get_spostDailyFile_name(time, sectorName)

        elif timeScale.lower() == "weekly":
            fileName = self.get_spostWeeklyFile_name(time, sectorName)

        elif timeScale.lower() == "monthly":
            fileName = self.get_spostMonthlyFile_name(time, sectorName)

        elif timeScale.lower() == "all":
            fileName = self.get_spostAllFile_name(time, sectorName)

        else:
            assert False

        emiss, LON, LAT = get_ncvar(fileName, ["emiss", "LON", "LAT"])

        return emiss, LON, LAT

    def get_truth_emiss(self, time, sectorName = None, timeScale = "orig"):

        assert self.isOSSE

        if sectorName is None:
            sectorName = self.sampleSector
        if timeScale.lower() == "orig":
            emiss, LON, LAT = get_ncvar(self.get_truthFile_name(time, sectorName), ["emiss", "LON", "LAT"])
        if timeScale.lower() == "hourly":
            emiss, LON, LAT = get_ncvar(self.get_truthHourlyFile_name(time, sectorName), ["emiss", "LON", "LAT"])
        if timeScale.lower() == "daily":
            emiss, LON, LAT = get_ncvar(self.get_truthDailyFile_name(time, sectorName), ["emiss", "LON", "LAT"])
        if timeScale.lower() == "weekly":
            emiss, LON, LAT = get_ncvar(self.get_truthWeeklyFile_name(time, sectorName), ["emiss", "LON", "LAT"])
        if timeScale.lower() == "monthly":
            emiss, LON, LAT = get_ncvar(self.get_truthMonthlyFile_name(time, sectorName), ["emiss", "LON", "LAT"])
        if timeScale.lower() == "all":
            emiss, LON, LAT = get_ncvar(self.get_truthAllFile_name(sectorName), ["emiss", "LON", "LAT"])

        return emiss, LON, LAT

    def get_reduced_uncertainty(self, time, recepTime, sectorName = None):

        if sectorName is None:
            sectorName = self.sampleSector

        ny, nx = self.nDim1_dic[sectorName], self.nDim2_dic[sectorName]
        LON, LAT = self.LON_dic[sectorName], self.LAT_dic[sectorName]

        fileName = self.get_dsigma2File_name(time, recepTime, sectorName)
        if os.path.exists(fileName):
            data = get_ncvar(fileName, "data")
        else:
            #print("Warning! file: " + fileName + " couldn't be found, return zeros matrix.")
            data = np.zeros(ny * nx)

        emiss = data.reshape(ny, nx)

        return emiss, LON, LAT

    def get_reduced_uncertainty_time(self, time, sectorName = None):
        if sectorName is None:
            sectorName = self.sectors[0]
        reducedUncertainty = 0
        for recepTime in self.iterTimeList:
            emiss, LON, LAT = self.get_reduced_uncertainty(time, recepTime, sectorName)
            reducedUncertainty = reducedUncertainty + emiss
        return reducedUncertainty, LON, LAT

    def get_obs(self, time):
        obs = get_ncvar(self.get_obsFile_name(recepTime = time), "data")
        obs = obs.flatten()
        return obs

    def get_priorHx(self, time):
        Hx = get_ncvar(self.get_priorHxFile_name(recepTime = time), "data")
        Hx = Hx.flatten()
        return Hx

    def get_procHx(self, time):
        Hx = get_ncvar(self.get_procHxFile_name(recepTime = time), "data")
        Hx = Hx.flatten()
        return Hx
        
    def get_priorHxBck(self, time):
        HxBck = get_ncvar(self.get_priorHxBckFile_name(recepTime = time), "data")
        HxBck = HxBck.flatten()
        return HxBck

    def get_procHxBck(self, time):
        HxBck = get_ncvar(self.get_procHxBckFile_name(recepTime = time), "data")
        HxBck = HxBck.flatten()
        return HxBck

    def get_initHx(self, time):
        Hx = get_ncvar(self.get_initHxFile_name(recepTime = time), "data")
        Hx = Hx.flatten()
        return Hx

    def get_initHxBck(self, time):
        HxBck = get_ncvar(self.get_initHxBckFile_name(recepTime = time), "data")
        HxBck = HxBck.flatten()
        return HxBck

    def get_finalHx(self, time):
        Hx = get_ncvar(self.get_finalHxFile_name(recepTime = time), "data")
        Hx = Hx.flatten()
        return Hx

    def get_finalHxBck(self, time):
        HxBck = get_ncvar(self.get_finalHxBckFile_name(recepTime = time), "data")
        HxBck = HxBck.flatten()
        return HxBck

    def get_bckPrior(self, time):
        if self.hasBCK:
            bck = get_ncvar(self.get_bckPriorFile_name(recepTime = time), "data")
            bck = bck.flatten()
        else:
            bck = np.zeros(len(self.obsLoc))
        return bck

    def get_sbckPrior(self, time):
        if self.hasBCK:
            sbck = get_ncvar(self.get_sbckPriorFile_name(recepTime = time), "data")
            sbck = sbck.flatten()
        else:
            sbck = np.zeros(len(self.obsLoc))
        return sbck

    def get_bckProc(self, time):
        if self.hasBCK:
            bck = get_ncvar(self.get_bckProcFile_name(recepTime = time), "data")
            bck = bck.flatten()
        else:
            bck = np.zeros(len(self.obsLoc))
        return bck

    def get_sbckProc(self, time):
        if self.hasBCK:
            sbck = get_ncvar(self.get_sbckProcFile_name(recepTime = time), "data")
            sbck = sbck.flatten()
        else:
            sbck = np.zeros(len(self.obsLoc))
        return sbck


    def get_recepTimelist(self, start, end, timeList_type, getFun, getFun_kwargs = {}, site_nameList = None):
        if site_nameList is None:
            site_nameList = self.obsName
        timeList = np.array(self.timeList_dic[timeList_type])
        timeList = timeList[(timeList >= start) & (timeList <= end)]
        value_timelist = []
        for time in timeList:
            value_thisTime = getFun(time, **getFun_kwargs)
            value_timelist.append(value_thisTime)
        value_timelist = np.array(value_timelist)
        value_timelist = value_timelist.T
        assert len(value_timelist.shape) == 2 and value_timelist.shape[0] == len(site_nameList)
        #value_timelist_dic = {}
        #for ind, siteName in enumerate(self.obsName):
        #    value_timelist_dic[siteName] = value_timelist[ind, :]
        A = {"data": (["site", "time"], value_timelist)}
        B = {"site": (["site"], site_nameList), "time": (["time"], timeList)}
        xrds = xr.Dataset(A, coords = B)
        return xrds

    def get_obs_timelist(self, start, end):
        getFun = self.get_obs
        timeList_type = "obs"
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, getFun = getFun)

    def get_priorHx_timelist(self, start, end):
        getFun = self.get_priorHx
        timeList_type = "obs"
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, getFun = getFun)

    def get_priorHxBck_timelist(self, start, end):
        getFun = self.get_priorHxBck
        timeList_type = "obs"
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, getFun = getFun)   

    def get_procHx_timelist(self, start, end):
        getFun = self.get_procHx
        timeList_type = "obs"
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, getFun = getFun)

    def get_procHxBck_timelist(self, start, end):
        getFun = self.get_procHxBck
        timeList_type = "obs"
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, getFun = getFun)

    def get_initHx_timelist(self, start, end):
        getFun = self.get_initHx
        timeList_type = "obs"
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, getFun = getFun)

    def get_initHxBck_timelist(self, start, end):
        getFun = self.get_initHxBck
        timeList_type = "obs"
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, getFun = getFun)

    def get_finalHx_timelist(self, start, end):
        getFun = self.get_finalHx
        timeList_type = "obs"
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, getFun = getFun)

    def get_finalHxBck_timelist(self, start, end):
        getFun = self.get_finalHxBck
        timeList_type = "obs"
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, getFun = getFun)

    def get_bckPrior_timelist(self, start, end):
        getFun = self.get_bckPrior
        timeList_type = "obs"
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, getFun = getFun)

    def get_bckProc_timelist(self, start, end):
        getFun = self.get_bckProc
        timeList_type = "obs"
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, getFun = getFun)

    def get_sbckPrior_timelist(self, start, end):
        getFun = self.get_sbckPrior
        timeList_type = "obs"
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, getFun = getFun)

    def get_sbckProc_timelist(self, start, end):
        getFun = self.get_sbckProc
        timeList_type = "obs"
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, getFun = getFun)

    #def get_truth_timelist(self, start, end)
    #    getFun = self.get_truth_emiss

    """
    def get_arr_points(self, func = None, **kwargs):
        arr, LON, LAT = func(**kwargs)
        pointsValues = []
        for point in self.obsPoints:
            point.anchor(LON, LAT)
            pointsValues.append(point.belinear(arr))
        pointsValues = np.array(pointsValues)
        return pointsValues
    """

    def get_arr_points(self, arr, LON, LAT, points_dic, interpMethod = "belinear"):
        pointsValues = []
        for pointName in points_dic:
            point = points_dic[pointName]
            point.anchor(LON, LAT)
            pointsValues.append(point.belinear(arr) if interpMethod.lower() == "belinear" else point.nearest(arr))
        pointsValues = np.array(pointsValues)
        return pointsValues

    def get_priorEmiss_points(self, time, sectorName = None, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        emiss, LON, LAT = self.get_prior_emiss(time, sectorName = sectorName, timeScale = timeScale)
        if points_dic is None:
            points_dic = self.obsPoints
        return self.get_arr_points(emiss, LON, LAT, points_dic, interpMethod = interpMethod)

    def get_priorSigma_points(self, time, sectorName = None, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        emiss, LON, LAT = self.get_prior_sigma(time, sectorName = sectorName, timeScale = timeScale)
        if points_dic is None:
            points_dic = self.obsPoints
        return self.get_arr_points(emiss, LON, LAT, points_dic, interpMethod = interpMethod)

    def get_posteriorEmiss_points(self, time, sectorName = None, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        emiss, LON, LAT = self.get_posterior_emiss(time, sectorName = sectorName, timeScale = timeScale)
        if points_dic is None:
            points_dic = self.obsPoints
        return self.get_arr_points(emiss, LON, LAT, points_dic, interpMethod = interpMethod)

    def get_posteriorSigma_points(self, time, sectorName = None, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        emiss, LON, LAT = self.get_posterior_sigma(time, sectorName = sectorName, timeScale = timeScale)
        if points_dic is None:
            points_dic = self.obsPoints
        return self.get_arr_points(emiss, LON, LAT, points_dic, interpMethod = interpMethod)

    def get_truthEmiss_points(self, time, sectorName = None, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        emiss, LON, LAT = self.get_truth_emiss(time, sectorName = sectorName, timeScale = timeScale)
        if points_dic is None:
            points_dic = self.obsPoints
        return self.get_arr_points(emiss, LON, LAT, points_dic, interpMethod = interpMethod)

    def get_priorEmissPoints_timelist(self, start, end, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        if points_dic is None:
            site_nameList = None
        else:
            site_nameList = list(points_dic.keys())
        getFun = self.get_priorEmiss_points
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = site_nameList, getFun = getFun, getFun_kwargs = {"timeScale": timeScale, "interpMethod": interpMethod, "points_dic": points_dic})

    def get_priorSigmaPoints_timelist(self, start, end, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        if points_dic is None:
            site_nameList = None
        else:
            site_nameList = list(points_dic.keys())
        getFun = self.get_priorSigma_points
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = site_nameList, getFun = getFun, getFun_kwargs = {"timeScale": timeScale, "interpMethod": interpMethod, "points_dic": points_dic})

    def get_posteriorEmissPoints_timelist(self, start, end, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        if points_dic is None:
            site_nameList = None
        else:
            site_nameList = list(points_dic.keys())
        getFun = self.get_posteriorEmiss_points
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = site_nameList, getFun = getFun, getFun_kwargs = {"timeScale": timeScale, "interpMethod": interpMethod, "points_dic": points_dic})

    def get_posteriorSigmaPoints_timelist(self, start, end, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        if points_dic is None:
            site_nameList = None
        else:
            site_nameList = list(points_dic.keys())
        getFun = self.get_posteriorSigma_points
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = site_nameList, getFun = getFun, getFun_kwargs = {"timeScale": timeScale, "interpMethod": interpMethod, "points_dic": points_dic})

    def get_truthEmissPoints_timelist(self, start, end, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        if points_dic is None:
            site_nameList = None
        else:
            site_nameList = list(points_dic.keys())
        getFun = self.get_truthEmiss_points
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = site_nameList, getFun = getFun, getFun_kwargs = {"timeScale": timeScale, "interpMethod": interpMethod, "points_dic": points_dic})

    def time_dimension_compute(self, func = None, start = None, end = None, dt = None, output = "total", **kwargs):

        if start is None:
            start = self.caseStart
        if end is None:
            end = self.caseEnd
        if dt is None:
            dt = self.caseDt

        assert output.lower() in ["total", "mean", "total2", "mean2"]
        
        result = 0
        nCount = 0

        current = start
        while(current <= end):
            temp = func(current, **kwargs)
            #temp = temp[0] if isinstance(temp, tuple) else temp
            if isinstance(temp, tuple):
                temp, LON, LAT = temp
            else:
                temp = temp
                LON = None
                LAT = None

            if output.lower() == "total2" or output.lower() == "mean2":
                result = result + temp**2
            else:
                result = result + temp


            nCount += 1
            current += dt

        if output.lower() == "total":
            result = result

        if output.lower() == "mean":
            result = result / nCount

        if output.lower() == "total2":
            result = np.sqrt(result)

        if output.lower() == "mean2":
            result = np.sqrt(result) / nCount

        if LON is None:
            return result
        else:
            return result, LON, LAT

    def time_dimension_total(self, **kwargs):
        return self.time_dimension_compute(output = "total", **kwargs)

    def time_dimension_mean(self, **kwargs):
        return self.time_dimension_compute(output = "mean", **kwargs)

    def time_dimension_total2(self, **kwargs):
        return self.time_dimension_compute(output = "total2", **kwargs)

    def time_dimension_mean2(self, **kwargs):
        return self.time_dimension_compute(output = "mean2", **kwargs)

    def spatial_dimension_compute(self, func = None, maskout = None, areaWeighted = True, areaMatrix = None, sector_areaS = None, output = "total", **kwargs):

        assert output.lower() in ["total", "mean", "total2", "mean2"]

        if sector_areaS is None:
            sector_areaS = self.sectors[0]

        result = func(**kwargs)
        result = result[0] if isinstance(result, tuple) else result

        if maskout is not None:
            assert result.shape == maskout.shape
            result = result * maskout
    
        if areaWeighted:
            if areaMatrix is None:
                areaMatrix = self.areaS_dic[sector_areaS]

            result  = result * areaMatrix

        if output.lower() == "total":
            result = np.nansum(result)

        if output.lower() == "mean":
            result = np.nanmean(result)

        if output.lower() == "total2":
            result = result **2
            result = np.nansum(result)
            result = np.sqrt(result)

        if output.lower() == "mean2":
            result = result **2
            result = np.nanmean(result)
            result = np.sqrt(result)
        
        if (output.lower() == "mean" or output.lower() == "mean2") and areaWeighted:
            result = result / np.nansum(areaMatrix)

        return result

    def spatial_dimension_total(self, **kwargs):
        return self.spatial_dimension_compute(output = "total", **kwargs)

    def spatial_dimension_total2(self, **kwargs):
        return self.spatial_dimension_compute(output = "total2", **kwargs)

    def spatial_dimension_mean(self, **kwargs):
        return self.spatial_dimension_compute(output = "mean", **kwargs)

    def spatial_dimension_mean2(self, **kwargs):
        return self.spatial_dimension_compute(output = "mean2", **kwargs)

    def get_priorEmiss_region(self, time, regions, sectorName = None, timeScale = "orig",  areaWeighted = True, sumOrMean = "sum"):
        #emiss, LON, LAT = self.get_prior_emiss(time, sectorName = sectorName, timeScale = timeScale)
        sumOrMean = sumOrMean.lower()
        assert sumOrMean in ["sum", "mean"]
        isSum = sumOrMean == "sum"
        isMean = sumOrMean == "mean"
        regionsValue = []
        for region in regions:
            objRegion = self.regionsMaskout[region]
            if isSum:
                emissValue = self.spatial_dimension_total(func = self.get_prior_emiss, maskout = objRegion.mask_array, time = time, sectorName = sectorName, timeScale = timeScale,  areaWeighted = areaWeighted)
            if isMean:
                emissValue = self.spatial_dimension_mean(func = self.get_prior_emiss, maskout = objRegion.mask_array, time = time, sectorName = sectorName, timeScale = timeScale,  areaWeighted = areaWeighted)
            regionsValue.append(emissValue)
        return regionsValue

    def get_priorSigma_region(self, time, regions, sectorName = None, timeScale = "orig", areaWeighted = True, sumOrMean = "sum"):
        #emiss, LON, LAT = self.get_prior_emiss(time, sectorName = sectorName, timeScale = timeScale)

        sumOrMean = sumOrMean.lower()
        assert sumOrMean in ["sum", "mean"]
        isSum = sumOrMean == "sum"
        isMean = sumOrMean == "mean"

        regionsValue = []
        for region in regions:
            objRegion = self.regionsMaskout[region]
            if isSum:
                emissValue = self.spatial_dimension_total2(func = self.get_prior_sigma, maskout = objRegion.mask_array, time = time, sectorName = sectorName, timeScale = timeScale, areaWeighted = areaWeighted)
            if isMean:
                emissValue = self.spatial_dimension_mean2(func = self.get_prior_sigma, maskout = objRegion.mask_array, time = time, sectorName = sectorName, timeScale = timeScale, areaWeighted = areaWeighted)
            regionsValue.append(emissValue)
        return regionsValue

    def get_posteriorEmiss_region(self, time, regions, sectorName = None, timeScale = "orig", areaWeighted = True, sumOrMean = "sum"):
        #emiss, LON, LAT = self.get_prior_emiss(time, sectorName = sectorName, timeScale = timeScale)

        sumOrMean = sumOrMean.lower()
        assert sumOrMean in ["sum", "mean"]
        isSum = sumOrMean == "sum"
        isMean = sumOrMean == "mean"

        regionsValue = []
        for region in regions:
            objRegion = self.regionsMaskout[region]
            if isSum:
                emissValue = self.spatial_dimension_total(func = self.get_posterior_emiss, maskout = objRegion.mask_array, time = time, sectorName = sectorName, timeScale = timeScale, areaWeighted = areaWeighted)
            if isMean:
                emissValue = self.spatial_dimension_mean(func = self.get_posterior_emiss, maskout = objRegion.mask_array, time = time, sectorName = sectorName, timeScale = timeScale, areaWeighted = areaWeighted)
            regionsValue.append(emissValue)
        return regionsValue

    def get_posteriorSigma_region(self, time, regions, sectorName = None, timeScale = "orig", areaWeighted = True, sumOrMean = "sum"):
        #emiss, LON, LAT = self.get_prior_emiss(time, sectorName = sectorName, timeScale = timeScale)

        sumOrMean = sumOrMean.lower()
        assert sumOrMean in ["sum", "mean"]
        isSum = sumOrMean == "sum"
        isMean = sumOrMean == "mean"

        regionsValue = []
        for region in regions:
            objRegion = self.regionsMaskout[region]
            if isSum:
                emissValue = self.spatial_dimension_total2(func = self.get_posterior_sigma, maskout = objRegion.mask_array, time = time, sectorName = sectorName, timeScale = timeScale, areaWeighted = areaWeighted)
            if isMean:
                emissValue = self.spatial_dimension_mean2(func = self.get_posterior_sigma, maskout = objRegion.mask_array, time = time, sectorName = sectorName, timeScale = timeScale, areaWeighted = areaWeighted)
            regionsValue.append(emissValue)
        return regionsValue

    def get_truthEmiss_region(self, time, regions, sectorName = None, timeScale = "orig", areaWeighted = True, sumOrMean = "sum"):

        sumOrMean = sumOrMean.lower()
        assert sumOrMean in ["sum", "mean"]
        isSum = sumOrMean == "sum"
        isMean = sumOrMean == "mean"

        regionsValue = []
        for region in regions:
            objRegion = self.regionsMaskout[region]
            if isSum:
                emissValue = self.spatial_dimension_total(func = self.get_truth_emiss, maskout = objRegion.mask_array, time = time, sectorName = sectorName, timeScale = timeScale, areaWeighted = areaWeighted)
            if isMean:
                emissValue = self.spatial_dimension_mean(func = self.get_truth_emiss, maskout = objRegion.mask_array, time = time, sectorName = sectorName, timeScale = timeScale, areaWeighted = areaWeighted)
            regionsValue.append(emissValue)
        return regionsValue


    def get_priorEmissRegion_timelist(self, start, end, regions, timeScale = "orig", areaWeighted = True, sumOrMean = "sum"):
        getFun = self.get_priorEmiss_region
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = regions, getFun = getFun, getFun_kwargs = {"regions": regions, "timeScale": timeScale, "areaWeighted": areaWeighted, "sumOrMean": sumOrMean})
        
    def get_priorSigmaRegion_timelist(self, start, end, regions, timeScale = "orig", areaWeighted = True, sumOrMean = "sum"):
        getFun = self.get_priorSigma_region
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = regions, getFun = getFun, getFun_kwargs = {"regions": regions, "timeScale": timeScale, "areaWeighted": areaWeighted, "sumOrMean": sumOrMean})

    def get_posteriorEmissRegion_timelist(self, start, end, regions, timeScale = "orig", areaWeighted = True, sumOrMean = "sum"):
        getFun = self.get_posteriorEmiss_region
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = regions, getFun = getFun, getFun_kwargs = {"regions": regions, "timeScale": timeScale, "areaWeighted": areaWeighted, "sumOrMean": sumOrMean})

    def get_posteriorSigmaRegion_timelist(self, start, end, regions, timeScale = "orig", areaWeighted = True, sumOrMean = "sum"):
        getFun = self.get_posteriorSigma_region
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = regions, getFun = getFun, getFun_kwargs = {"regions": regions, "timeScale": timeScale, "areaWeighted": areaWeighted, "sumOrMean": sumOrMean})

    def get_truthEmissRegion_timelist(self, start, end, regions, timeScale = "orig", areaWeighted = True, sumOrMean = "sum"):
        getFun = self.get_truthEmiss_region
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = regions, getFun = getFun, getFun_kwargs = {"regions": regions, "timeScale": timeScale, "areaWeighted": areaWeighted, "sumOrMean": sumOrMean})
