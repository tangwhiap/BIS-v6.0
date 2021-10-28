#!/usr/bin/env python

# A class combines the functions of BIS output analysis. #
# The base class Bayesian Inversion System output analysis. #

# Authors:
#   Wenhan TANG - 05/2021
#   ...

import numpy as np
import netCDF4 as nc
import datetime as dtm
import pandas as pd
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from ..utils.MaskTools import Maskout
import glob
import os
from pdb import set_trace

Prefix_processing = "BISOUT"
Prefix_sigma = "Sigma"
ShpDir = "/home/tangwh/datasets/china_shp/cnhimap.shp"
CasesDir = "/home/tangwh/modeldata/BIS_cases"

class BIS_Case_Analyst(object):

    def __init__(self, CaseDir = None, Start = None, End = None, CaseName = "BIS case", FolderName = None, dt = None, window = 6, n_tback = 24):

        if CaseDir is None:
            CaseDir = CasesDir + "/" + (CaseName if FolderName is None else FolderName)

        self.CaseDir = CaseDir
        self.CaseName = CaseName

        if Start is None:
            fileList = glob.glob(self.CaseDir + "/output/final/" + Prefix_processing + "_*")
            print(self.CaseDir + "/output/final/" + Prefix_processing + "_*")
            fileList.sort()
            assert len(fileList) >= 1
            Start = fileList[0].strip().split("/")[-1].split(".")[0].split("_")[1:]
            Start = Start[0] + "_" + Start[1]

        if End is None:
            fileList = glob.glob(self.CaseDir + "/output/final/" + Prefix_processing + "_*")
            fileList.sort()
            assert len(fileList) >= 1
            End = fileList[-1].strip().split("/")[-1].split(".")[0].split("_")[1:]
            End = End[0] + "_" + End[1]

        #print(Start, End)

        if dt is None:
            fileList = glob.glob(self.CaseDir + "/output/final/" + Prefix_processing + "_*")
            fileList.sort()
            assert len(fileList) >= 1
            file1 = fileList[0]; file2 = fileList[1]
            time1 = file1.strip().split("/")[-1].split(".")[0].split("_")[1:]
            time1 = time1[0] + "_" + time1[1]
            time2 = file2.strip().split("/")[-1].split(".")[0].split("_")[1:]
            time2 = time2[0] + "_" + time2[1]
            time1 = dtm.datetime.strptime(time1, "%Y-%m-%d_%H:%M:%S")
            time2 = dtm.datetime.strptime(time2, "%Y-%m-%d_%H:%M:%S")
            dt = (time2 - time1).seconds // 3600

        if type(Start) == type("string") and type(End) == type("string"):
            self.Start = dtm.datetime.strptime(Start, "%Y-%m-%d_%H:%M:%S")
            self.End = dtm.datetime.strptime(End, "%Y-%m-%d_%H:%M:%S")

        else:
            assert type(Start) == type(End) and type(End) == type(dtm.datetime(1,1,1))
            self.Start = Start
            self.End = End
 
        self.dt = dt
        self.dtout = dtm.timedelta(hours = dt) 
        self.window = window
        self.Twin = dtm.timedelta(hours = self.window)
        self.n_tback = n_tback
        self.tback = dtm.timedelta(hours = n_tback)
        self.prior_mean_2D = None
        self.posterior_mean_2D = None
        self.total_emiss_1D = None
        self.info_check()
        self.get_Domain_info()
        self.get_OBS_info()
        #self.mask_BJ = Maskout(lonlist = self.lon, latlist = self.lat, Region_list = ["Beijing"])
        self.Region_dic = {}

        print("Successfully create BIS-case analyst object for \"" + self.CaseName + "\" from " + str(Start) + " to " + str(End) + " .")


    def get_HQ_dir(self, RecepTime, back_time):
        return self.CaseDir + "/history_save/" + RecepTime.strftime("%Y%m%d%H") + "/HQ_" + back_time.strftime("%Y%m%d%H") + ".nc"

    def get_RQ_dir(self, RecepTime, back_time):
        return self.CaseDir + "/history_save/" + RecepTime.strftime("%Y%m%d%H") + "/RQ_" + back_time.strftime("%Y%m%d%H") + ".nc"

    def get_OBS_dir(self, RecepTime):
        return self.CaseDir + "/obs/OBS_" + RecepTime.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_INV_dir(self, RecepTime):
        return self.CaseDir + "/history_save/" + RecepTime.strftime("%Y%m%d%H") + "/INV.nc"

    def get_Hx_dir(self, RecepTime, HxType = "orig"):
        HxType = HxType.lower()
        assert HxType in ["orig", "pos", "final"]
        if HxType == "orig":
            Suffix = ""
        if HxType == "pos":
            Suffix = "_pos"
        if HxType == "final":
            Suffix = "_final"
        return self.CaseDir + "/history_save/" + RecepTime.strftime("%Y%m%d%H") + "/Hx" + Suffix + ".nc"

    def get_footprint_dir(self, RecepName, RecepTime):
        foot_dir = glob.glob(self.CaseDir + "/footprint/" + RecepName + "/by-id/" + RecepTime.strftime("%Y%m%d%H%M") + "_*_*_*/" + RecepTime.strftime("%Y%m%d%H%M") + "_*_*_*_foot.nc")
        assert len(foot_dir) == 1
        foot_dir = foot_dir[0]
        return foot_dir

    def get_posterior_dir(self, time):
        #Prefix_processing = "BISOUT"
        #return self.CaseDir + "/output/final/"  + Prefix_processing + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"
        # modify later
        return self.CaseDir + "/output/processing/"  + Prefix_processing + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"
    
    def get_priorSigma_dir(self, time):
        return self.CaseDir + "/sigma_prior/" + Prefix_sigma + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_Domain_info(self):

        def areaS(lon, lat):
            R = 6371000
            dlon = (lon[1:] - lon[:-1]).mean()
            dlat = (lat[1:] - lat[:-1]).mean()
            LON, LAT = np.meshgrid(lon, lat)
            S = (R**2) * np.deg2rad(dlon) * (np.sin(np.deg2rad(LAT + dlat/2 )) - np.sin(np.deg2rad(LAT - dlat/2 )))
            assert S.shape[0] == len(lat)
            assert S.shape[1] == len(lon)
            return S
        
        SampleFile = self.get_posterior_dir(self.Start)
        ncf = nc.Dataset(SampleFile, "r")
        self.nlon = ncf.dimensions["lon"].size
        self.nlat = ncf.dimensions["lat"].size
        self.Ng = self.nlon * self.nlat
        self.lon = ncf.variables["lon"][:].filled(np.nan)
        self.lat = ncf.variables["lat"][:].filled(np.nan)
        self.area = areaS(self.lon, self.lat)
        ncf.close()

    def get_OBS_info(self):

        #print(self.receptor_times)

        SampleFile = self.get_OBS_dir(self.receptor_times[0])
        ncf = nc.Dataset(SampleFile, "r")
        self.Ns = ncf.dimensions["sites"].size
        self.ObsName_list = ncf.variables["sites"][:].tolist()
        ncf.close()
        self.site_location_dic = {}
        for site in self.ObsName_list:
            glob_list = glob.glob(self.CaseDir + "/footprint/" + site + "/by-id/" + self.receptor_times[0].strftime("%Y%m%d%H%M") + "_*")
            glob_demo = glob_list[0].strip().split("/")[-1].split("_")
            site_lon = float(glob_demo[1])
            site_lat = float(glob_demo[2])
            site_height = float(glob_demo[3])
            self.site_location_dic[site] = (site_lon, site_lat, site_height)
            
            
    
    """    
    def get_prior_dir(self):
        return self.CaseDir + "/emiss/Prior.nc"
    """

    def get_prior_dir(self, time):
        return self.CaseDir + "/output/prior/Prior_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"
        
    """
    def get_prior_emiss(self, time):
        ncf = nc.Dataset(self.get_prior_dir(), "r")
        ftime = ncf.variables["time"][:].filled(np.nan)
        utime = ncf.variables["time"].units
        numtime = nc.date2num(time, units = utime)
        emiss = ncf.variables["emiss"][ftime == numtime, :, :].filled(np.nan)[0, :, :]
        assert emiss.shape[0] == self.nlat
        assert emiss.shape[1] == self.nlon
        return emiss
    """

    def get_prior_emiss(self, time):
        ncf = nc.Dataset(self.get_prior_dir(time), "r")
        emiss = ncf.variables["emiss"][:, :].filled(np.nan)
        assert emiss.shape[0] == self.nlat
        assert emiss.shape[1] == self.nlon
        return emiss

    def get_posterior_emiss(self, time):
        #print(self.get_posterior_dir(time))
        ncf = nc.Dataset(self.get_posterior_dir(time), "r")
        emiss = ncf.variables["emiss"][:, :].filled(np.nan)
        assert emiss.shape[0] == self.nlat
        assert emiss.shape[1] == self.nlon
        return emiss

    def get_RQ_dirList(self, time = None, RecepTime = None, assert_founded = False):
        #assert time is not None or RecepTime is not None
        FileDir = self.CaseDir + "/history_save/"
        FileDir += "*/" if RecepTime is None else RecepTime.strftime("%Y%m%d%H") + "/"
        FileDir += "RQ_" + ("*" if time is None else time.strftime("%Y%m%d%H")) + ".nc"
        dirList = glob.glob(FileDir)
        if assert_founded:
            assert len(dirList) > 0, FileDir + " couldn't be found!"
        #if len(dirList) == 1:
        #    dirList = dirList[0]
        return dirList
    
    def get_dEmiss_dirList(self, time = None, RecepTime = None, assert_founded = False):
        #assert time is not None or RecepTime is not None
        FileDir = self.CaseDir + "/history_save/"
        FileDir += "*/" if RecepTime is None else RecepTime.strftime("%Y%m%d%H") + "/"
        FileDir += "dEmiss_" + ("*" if time is None else time.strftime("%Y%m%d%H")) + ".nc"
        dirList = glob.glob(FileDir)
        if assert_founded:
            assert len(dirList) > 0, FileDir + " couldn't be found!"
        #if len(dirList) == 1:
        #    dirList = dirList[0]
        return dirList

    def get_priorSigma(self, time):
        ncf = nc.Dataset(self.get_priorSigma_dir(time), "r")
        emiss = ncf.variables["emiss"][:, :].filled(np.nan)
        assert emiss.shape[0] == self.nlat
        assert emiss.shape[1] == self.nlon
        return emiss
    
    def get_posteriorSigma(self, time):
        prior_sigma = self.get_priorSigma(time)
        try:
            reduced_unc = self.get_reduced_total_unc(time)[0]
        except:
            reduced_unc = prior_sigma * 0
        #if np.sum(reduced_unc<0) > 0:
        #    set_trace()
        posterior_variance = prior_sigma ** 2 - reduced_unc
        N_lower0 = np.sum(posterior_variance < 0)
        if N_lower0 > 0:
            #print("Warning! There are " + str(N_lower0) + " points where the posterior variance is lower than 0.")
            posterior_variance = np.where(posterior_variance > 0, posterior_variance, 0)
        posterior_sigma = posterior_variance ** 0.5
        return posterior_sigma

    def get_reduced_unc(self, time = None, RecepTime = None):
        dirList = self.get_RQ_dirList(time = time, RecepTime = RecepTime)
        #set_trace()
        RQ_list = []
        for ifile in dirList:

            time = dtm.datetime.strptime(ifile.strip().split("/")[-1], "RQ_%Y%m%d%H.nc")
            RecepTime = dtm.datetime.strptime(ifile.strip().split("/")[-2], "%Y%m%d%H")
            ncf = nc.Dataset(ifile, "r")
            RQ_value = ncf.variables["reduced"][:, :].filled(np.nan)
            RQ_list.append({"RQ": RQ_value, "time": time, "RecepTime": RecepTime})

        return RQ_list

    def get_dEmiss(self, time = None, RecepTime = None):
        dirList = self.get_dEmiss_dirList(time = time, RecepTime = RecepTime)
        #set_trace()
        dEmiss_list = []
        for ifile in dirList:
            time = dtm.datetime.strptime(ifile.strip().split("/")[-1], "dEmiss_%Y%m%d%H.nc")
            RecepTime = dtm.datetime.strptime(ifile.strip().split("/")[-2], "%Y%m%d%H")
            ncf = nc.Dataset(ifile, "r")
            dEmiss_value = ncf.variables["emiss"][:, :].filled(np.nan)
            dEmiss_list.append({"dEmiss": dEmiss_value, "time": time, "RecepTime": RecepTime})
        return dEmiss_list

    def get_reduced_total_unc(self, time):
        #if time is None:
        #    time = self.pos_emiss_times
        RQ_list = self.get_reduced_unc(time = time)
        Reduced_times = len(RQ_list)
        #RQ_total = np.zeros_like(RQ_list[0]["RQ"])
        RQ_total = np.zeros((self.nlat, self.nlon))
        for iRQ in RQ_list:
            RQ_total += iRQ["RQ"]
        return RQ_total, Reduced_times

    def get_reduced_total_unc_mean(self, Start = None, End = None, dt = None):
        if Start is None:
            Start = self.Start
        if End is None:
            End = self.End
        if dt is None:
            dt = self.dt
        tlist = pd.date_range(Start, End, freq = str(dt) + "H")
        #RQ_mean = np.full( (self.nlat, self.nlon), 0.0 )
        RQ_mean = []
        for itime in tlist:
            RQ, _ = self.get_reduced_total_unc(time = itime)
            RQ_mean.append(RQ)
        RQ_mean = np.array(RQ_mean)
        RQ_mean = RQ_mean.mean(axis = 0)
        return RQ_mean
        
        
    def total_RQ_timelist(self, Start = None, End = None, dt = None):
        if Start is None:
            Start = self.Start
        if End is None:
            End = self.End
        if dt is None:
            dt = self.dt
        tlist = pd.date_range(Start, End, freq = str(dt) + "H")
        RQtot_list = []
        for itime in tlist:
            assert itime in self.pos_emiss_times
            RQ_total, _ = self.get_reduced_total_unc(time = itime)
            #emiss_total = (emiss * self.area).sum() * 3600 / 1000000 * 44 / 1000000
            RQtot_list.append((RQ_total * self.area).sum() * 3600 / 1000000 * 44 / 1000000)
        RQtot_list = np.array(RQtot_list)

        return tlist, RQtot_list
            
    
    def get_footprint(self, Receptors_list = None, time_list = None, RecepTime_list = None):
        #assert time is not None or RecepTime is not None and Receptor is not None
        spc_Receptors = Receptors_list is not None
        spc_time = time_list is not None
        spc_RecepTime = RecepTime_list is not None
        if not spc_Receptors:
            Receptors_list = self.ObsName_list
        if not spc_time:
            time_list = self.pos_emiss_times
        if not RecepTime_list:
            RecepTime_list = self.receptor_times
        foot_receptors = []
        for receptor in Receptors_list:
            foot_receptor_mean = []
            for RecepTime in RecepTime_list:
                foot_dir = self.get_footprint_dir(RecepName = receptor, RecepTime = RecepTime)
                #print("Open ", foot_dir)
                ncf_foot = nc.Dataset(foot_dir, "r")
                Valid_foot_time_range = self.recep_backtime_dic[RecepTime.strftime("%Y%m%d%H")]
                Valid_time_range = []
                for time in time_list:
                    #set_trace()
                    if time in Valid_foot_time_range:
                        Valid_time_range.append(time)
                if len(Valid_time_range) == 0:
                    #print("Not contain, go on")
                    continue
                #Valid_time_range = time_list[ [time in Valid_foot_time_range for time in time_list] ]
                foot = np.full( (self.nlat, self.nlon), 0.0 )
                for time in Valid_time_range:
                    tunits = ncf_foot.variables["time"].units
                    tcalendar = ncf_foot.variables["time"].calendar
                    tv = nc.date2num(time, units = tunits, calendar = tcalendar)
                    ftv = ncf_foot.variables["time"][:].filled(np.nan)
                    tid = np.where(ftv == tv)[0]
                    if len(tid) == 0:
                        #print("pass")
                        continue
                    foot += ncf_foot.variables["foot"][np.where(ftv == tv)[0], :, :].filled(np.nan)[0, :, :]
                    #print("get")
                ncf_foot.close()
                foot_receptor_mean.append(foot)
            foot_receptor_mean = np.array(foot_receptor_mean)
            #set_trace()
            foot_receptor_mean = foot_receptor_mean.mean(axis = 0)  
            foot_receptors.append(foot_receptor_mean)
        foot_receptors = np.array(foot_receptors)
        foot_receptors = foot_receptors.sum(axis = 0)
        return foot_receptors
        
        
    def get_OBS(self, timeList = None, receptors_list = None):
        #print("int get_OBS", time)
        if timeList is None:
            timeList = self.receptor_times
        if isinstance(timeList, dtm.datetime) or isinstance(timeList, pd.Timestamp):
            timeList = [timeList]
    
        if receptors_list is None:
            receptors_list = self.ObsName_list

        if type(receptors_list) == type("string"):
            receptors_list = [receptors_list]
        #sid = list(map(lambda x: x in receptors_list, sites_list))
        obs_dic = {}
        recep_ind_list = []
        for receptor in receptors_list:
            recep_ind = self.ObsName_list.index(receptor)
            recep_ind_list.append(recep_ind)
            obs_dic[receptor] = []
        recep_ind_list = np.array(recep_ind_list)

        for time in timeList:
        
            ncf_obs = nc.Dataset(self.CaseDir + "/obs/OBS_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc", "r")
            #sites_list = ncf_obs.variables["sites"][:]
            obs = ncf_obs.variables["obs"][recep_ind_list].filled(np.nan)
            ncf_obs.close()
            for irecep, receptor in enumerate(receptors_list):
                obs_dic[receptor].append(obs[irecep])
        for receptor in receptors_list:
            obs_dic[receptor] = np.array(obs_dic[receptor])
            #if receptors_list is not None:
            #    obs = ncf_obs.variables["obs"][sid].filled(np.nan)
            #else:
            #    obs = ncf_obs.variables["obs"][:].filled(np.nan)
            #    receptors_list = ncf_obs.variables["sites"][:]
        return timeList, obs_dic

    def get_Hx(self, receptors_list = None, HxType = "orig"):
        HxType = HxType.lower()
        assert HxType in ["orig", "pos", "final"]
        if isinstance(receptors_list, str):
            receptors_list = [receptors_list]
        recep_ind_list = []
        Hx_dic = {}
        for receptor in receptors_list:
            recep_ind = self.ObsName_list.index(receptor)
            recep_ind_list.append(recep_ind)
            Hx_dic[receptor] = []
        recep_ind_list = np.array(recep_ind_list)
        for time in self.receptor_times:
            HxDir = self.get_Hx_dir(time, HxType = HxType)
            ncf = nc.Dataset(HxDir, "r")
            Hx = ncf.variables["Hx"][recep_ind_list]
            ncf.close()
            for irecep, receptor in enumerate(receptors_list):
                Hx_dic[receptor].append(Hx[irecep])
        for receptor in receptors_list:
            Hx_dic[receptor] = np.array(Hx_dic[receptor])

        return np.array(self.receptor_times), Hx_dic



    def draw_OBS_scatter(self, time, fig, receptors_list = None, lonlist = None, latlist = None, cmap = "jet", subplot = (1,1,1), draw_gridline = True, draw_lonlat_label = True, draw_cbar = True, size = 100, overlap_size_ratio = 0.3, cbar_orientation = "horizontal", cbar_label = "ppm", vmin = None, vmax = None, **kwarg):

        #def alpha_vary_cmap(colormap):
        #    cmap = plt.get_cmap(colormap)
        #    my_cmap = cmap(np.arange(cmap.N))
        #    my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
        #    my_cmap = ListedColormap(my_cmap)
        #    return my_cmapobs_dict = self.get_OBS(time, receptors_list)

        sites_lon = []
        sites_lat = []
        sites_obs = []
        is_overlap = []
        #print("in draw_OBS_scatter", time)
        _, obs_dict = self.get_OBS(time, receptors_list)
        for site in obs_dict:
            if np.isnan(obs_dict[site]):
                continue
            slon = self.site_location_dic[site][0]
            slat = self.site_location_dic[site][1]
            is_overlap.append( ~(np.sum(np.abs(slon - np.array(sites_lon) < 0.02) & (np.abs(slat - np.array(sites_lat)) < 0.02)) == 0) )
            sites_lon.append(slon)
            sites_lat.append(slat)
            sites_obs.append(obs_dict[site])

        sites_lon = np.array(sites_lon)
        sites_lat = np.array(sites_lat)
        sites_obs = np.array(sites_obs)
        is_overlap = np.array(is_overlap)
        #set_trace()
        #print(is_overlap)

        #set_trace()
        if lonlist is None:
            lonlist = self.lon
        if latlist is None:
            latlist = self.lat
        
        #assert len(Arr.shape) == 2 and Arr.shape[0] == len(latlist) and Arr.shape[1] == len(lonlist)
        proj = ccrs.PlateCarree()
        geo_ax = fig.add_subplot(subplot[0], subplot[1], subplot[2], projection = proj)
        lon_s = lonlist[0]; lon_e = lonlist[-1]
        lat_s = latlist[0]; lat_e = latlist[-1]
        geo_ax.set_extent([lon_s, lon_e, lat_s, lat_e], crs = proj)
        if draw_gridline:
            gl = geo_ax.gridlines(crs = proj, linestyle = "--", alpha = 0.5, draw_labels = True)
            gl.top_labels = False
            gl.right_labels = False
            if not draw_lonlat_label:
                gl.left_lables = False
                gl.bottom_lables = False
        geo_ax.add_geometries(Reader(ShpDir).geometries(), proj, facecolor = "none", edgecolor = "k", linewidth = 1)
        if vmin is None:
            vmin = sites_obs.min()
        if vmax is None:
            vmax = sites_obs.max()
        cs1 = geo_ax.scatter(sites_lon[~is_overlap], sites_lat[~is_overlap], s = size, c = sites_obs[~is_overlap], cmap = cmap, vmin = vmin, vmax = vmax, edgecolor = "white", **kwarg)
        cs2 = geo_ax.scatter(sites_lon[is_overlap], sites_lat[is_overlap], s = size * overlap_size_ratio, c = sites_obs[is_overlap], cmap = cmap, vmin = vmin, vmax = vmax, edgecolor = "white", **kwarg)
        if draw_cbar:
            cbar = fig.colorbar(cs1, ax = geo_ax, orientation = cbar_orientation)
            cbar.set_label(cbar_label)
        #if cmap_va:
        #    cmap = alpha_vary_cmap(cmap)
        return geo_ax, cs1, cs2

    def info_check(self):
        assert os.path.exists(self.CaseDir), "Case directory: " + self.CaseDir + " is not exist!" 
        assert os.path.exists(self.CaseDir + "/output/final"), self.CaseDir + "/output/final is not exist!"
        assert os.path.exists(self.CaseDir + "/history_save"), self.CaseDir + "/history_save is not exist!"
        self.receptor_times = pd.date_range(self.Start + self.tback, self.End, freq = str(self.window) + "H")
        self.pos_emiss_times = pd.date_range(self.Start, self.End, freq = str(self.dt) + "H")

        self.recep_backtime_dic = {}
        for itime in self.receptor_times:
            assert os.path.exists(self.CaseDir + "/history_save/" + itime.strftime("%Y%m%d%H")), self.CaseDir + "/history_save/" + itime.strftime("%Y%m%d%H") + " is not exist!"
            
            back_time_list = pd.date_range(itime.to_pydatetime() - self.tback, itime.to_pydatetime(), freq = str(self.dt) + "H")[:-1]
            #print(len(back_time_list))
            self.recep_backtime_dic[itime.strftime("%Y%m%d%H")] = back_time_list
            for back_time in back_time_list:
                assert os.path.exists(self.get_HQ_dir(itime, back_time)), self.get_HQ_dir(itime, back_time) + " is not exist!"
                #assert os.path.exists(self.get_RQ_dir(itime, back_time))
            assert os.path.exists(self.get_INV_dir(itime))
            assert os.path.exists(self.get_OBS_dir(itime))

        for itime in self.pos_emiss_times:
            assert os.path.exists(self.get_posterior_dir(itime)), self.get_posterior_dir(itime) + " is not exist!"

        self.opti_recepTime_4eachTime = {}
        for itime in self.pos_emiss_times:
            self.opti_recepTime_4eachTime[itime.strftime("%Y%m%d%H")] = []
        for itime in self.receptor_times:
            backtime = self.recep_backtime_dic[itime.strftime("%Y%m%d%H")]
            for btime in backtime:
                self.opti_recepTime_4eachTime[btime.strftime("%Y%m%d%H")].append(itime)
                
    def prior_mean(self, start = None, end = None, dt = None):

        if self.prior_mean_2D is not None:
            return self.prior_mean_2D

        if start == None:
            start = self.Start

        if end == None:
            end = self.End

        if dt == None:
            dt = self.dt

        emiss_total = np.full( (self.nlat, self.nlon), 0.0 ) 
        Ntime = 0

        dt = dtm.timedelta(hours = dt)
        Current = start
        while( Current <= end ):
            emiss_total += self.get_prior_emiss(Current)
            Ntime += 1
            Current += dt 

        emiss_mean = emiss_total / Ntime
        self.prior_mean_2D = emiss_mean

        return emiss_mean

    def posterior_mean(self, start = None, end = None, dt = None):
        if self.posterior_mean_2D is not None:
            return self.posterior_mean_2D

        if start == None:
            start = self.Start
        if end == None:
            end = self.End
        if dt == None:
            dt = self.dt
        emiss_total = np.full( (self.nlat, self.nlon), 0.0 ) 
        Ntime = 0
        dt = dtm.timedelta(hours = dt)
        Current = start
        while( Current <= end ):
            emiss_total += self.get_posterior_emiss(Current)
            Ntime += 1
            Current += dt 
        emiss_mean = emiss_total / Ntime
        self.posterior_mean_2D = emiss_mean
        return emiss_mean

    def emiss_point_value(self, Time, plon, plat, EmissType = "posterior", get_variance = False):
        EmissType = EmissType.lower()
        assert EmissType in ["prior", "posterior"]
        dif_lon = np.abs(plon - self.lon)
        ix = int(np.where(dif_lon == dif_lon.min())[0][0])
        dif_lat = np.abs(plat - self.lat)
        iy = int(np.where(dif_lat == dif_lat.min())[0][0])
        getFun = self.get_posterior_emiss if EmissType == "posterior" else self.get_prior_emiss
        emiss = getFun(Time)[iy, ix]
        if get_variance:
            sigma = self.get_priorSigma(Time)[iy, ix]
            variance = sigma ** 2
            if EmissType == "posterior":
                try:
                   unc = self.get_reduced_total_unc(Time)[0][iy, ix]
                except:
                    print("Warning! Reduced uncertainty on ", Time, " couln'd be found, using 0 alternatively.")
                    unc = 0
                variance = variance - unc
                if variance <0:
                    print("Warning: variance of point: (", plon, ", ", plat, ") is ", variance, ", which is lower than 0, change to 0.")
                    variance = 0
                
        return emiss, variance if get_variance else None

    def get_SiteLoc_emiss(self, start = None, end = None, dt = None, EmissType = "posterior", get_variance = False):
        EmissType = EmissType.lower()
        assert EmissType in ["prior", "posterior"]
        if start == None:
            start = self.Start
        if end == None:
            end = self.End
        if dt == None:
            dt = self.dt

        timeList = []
        emiss_dic = {}
        if get_variance:
            variance_dic = {}
        for station in self.site_location_dic:
            emiss_dic[station] = []
            if get_variance:
                variance_dic[station] = []

        Current = start
        while(Current <= end):
            timeList.append(Current)
            for station in self.site_location_dic:
                plon = self.site_location_dic[station][0]
                plat = self.site_location_dic[station][1]
                emiss_value, variance_value = self.emiss_point_value(Time = Current, plon = plon, plat = plat, EmissType = EmissType, get_variance = get_variance)
                emiss_dic[station].append(emiss_value)
                if get_variance:
                    variance_dic[station].append(variance_value)
            Current += dtm.timedelta(hours = dt)

        for station in self.site_location_dic:
            emiss_dic[station] = np.array(emiss_dic[station])
            if get_variance:
                variance_dic[station] = np.array(variance_dic[station])
        timeList = np.array(timeList)
        return timeList, emiss_dic, variance_dic if get_variance else None

        
    def total_emiss_timelist(self, emiss_type = None, get_emiss_fun = None, start = None, end = None, dt = None, Region = None, sigma_out = False):
        if self.total_emiss_1D is not None:
            return self.total_emiss_1D

        if start == None:
            start = self.Start
        if end == None:
            end = self.End
        if dt == None:
            dt = self.dt
        if isinstance(start, str):
            start = dtm.datetime.strptime(start, "%Y-%m-%d_%H:%M:%S")
        if isinstance(end, str):
            end = dtm.datetime.strptime(end, "%Y-%m-%d_%H:%M:%S")
        if type(emiss_type) == type("string"):
            emiss_type = emiss_type.lower()
        isfun = not(emiss_type in ["posterior", "prior"])
        if isfun:
            assert get_emiss_fun is not None
            assert sigma_out == False # Sigma_out option is prepared only for posterior and prior emission.
        dt = dtm.timedelta(hours = dt)
        if Region is None:
            mask_array = np.full( (self.nlat, self.nlon), 1 )
        else:
            if Region not in self.Region_dic:
                self.Region_dic[Region] = Maskout(lonlist = self.lon, latlist = self.lat, Region_list = [Region])
            mask_array = self.Region_dic[Region].mask_array_out()

        xtime = []
        yemiss = []
        if sigma_out:
            ysigma = []
        Current = start 
        while( Current <= end ):
            if isfun:
                loc = locals()
                exec("emiss = self." + get_emiss_fun + "(Current)")
                emiss = loc["emiss"]
            else:
                emiss = self.get_posterior_emiss(Current) if emiss_type == "posterior" else self.get_prior_emiss(Current)
            if sigma_out:
                sigma = self.get_posteriorSigma(Current) if emiss_type == "posterior" else self.get_priorSigma(Current)
            emiss = emiss * mask_array
            if sigma_out:
                sigma = sigma * mask_array
            # umol/m^2/s --> tCO2
            emiss_total = (emiss * self.area).sum() * 3600 / 1000000 * 44 / 1000000
            sigma_total = (sigma * self.area).sum() * 3600 / 1000000 * 44 / 1000000
            xtime.append(np.datetime64(Current))
            yemiss.append(emiss_total)
            ysigma.append(sigma_total)
            Current += dt
        xtime = np.array(xtime)
        yemiss = np.array(yemiss)
        ysigma = np.array(ysigma)
        if sigma_out:
            return xtime, yemiss, ysigma
        else:
            return xtime, yemiss
    
    def plot2D(self, Arr, fig, lonlist = None, latlist = None, subplot = (1, 1, 1), cmap = "jet", cmap_va = False, draw_gridline = True, draw_lonlat_label = True, **plot_kwargs):

        if lonlist is None:
            lonlist = self.lon
        if latlist is None:
            latlist = self.lat

        def alpha_vary_cmap(colormap):
            cmap = plt.get_cmap(colormap)
            my_cmap = cmap(np.arange(cmap.N))
            my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
            my_cmap = ListedColormap(my_cmap)
            return my_cmap
        assert len(Arr.shape) == 2 and Arr.shape[0] == len(latlist) and Arr.shape[1] == len(lonlist)
        proj = ccrs.PlateCarree()
        geo_ax = fig.add_subplot(subplot[0], subplot[1], subplot[2], projection = proj)
        lon_s = lonlist[0]; lon_e = lonlist[-1]
        lat_s = latlist[0]; lat_e = latlist[-1]
        geo_ax.set_extent([lon_s, lon_e, lat_s, lat_e], crs = proj)
        if draw_gridline:
            gl = geo_ax.gridlines(crs = proj, linestyle = "--", alpha = 0.5, draw_labels = True)
            gl.top_labels = False
            gl.right_labels = False
            if not draw_lonlat_label:
                gl.left_lables = False
                gl.bottom_lables = False
        geo_ax.add_geometries(Reader(ShpDir).geometries(), proj, facecolor = "none", edgecolor = "k", linewidth = 1)
        if cmap_va:
            cmap = alpha_vary_cmap(cmap)
        cs = geo_ax.pcolormesh(lonlist, latlist, Arr, cmap = cmap, **plot_kwargs)
        return geo_ax, cs
        
    def plot_tlist(self, timelist, Arr, fig = None, ax = None, subplot = (1, 1, 1), **plot_kwargs):
        assert len(timelist) == len(Arr)
        if ax is None:
            assert fig is not None, "At least one of the keyword \"fig\" or \"ax\" should be given"
            ax = fig.add_subplot( subplot[0], subplot[1], subplot[2] )
        line = ax.plot(timelist, Arr, **plot_kwargs)
        return ax, line

    def percentile_range(self, arr, p = 100):
        perc1 = (100 - p) / 2
        perc2 = 100 - perc1
        vr1 = np.percentile(arr, perc1)
        vr2 = np.percentile(arr, perc2)
        return vr1, vr2

    def center0_adjust(self, vrange, force = False):
        assert len(vrange) == 2 and vrange[0] <= vrange[1]
        if force or (vrange[0] * vrange[1] < 0):
            v1 = np.abs(vrange[0])
            v2 = np.abs(vrange[1])
            vm = np.max([v1, v2])
            vrange = (-vm, vm)
        else:
            print("Warning!, vrange:", vrange, " isn't need to be adjust, you can set force = True")
        return vrange

    def get_vrange_OBS(self, timeList = None, receptors_list = None, p = 100):
        """
        if timeList is None:
            timeList = self.receptor_times
        if receptors_list is None:
            receptors_list = self.ObsName_list
        """
        _, obs_dic = self.get_OBS(timeList = timeList, receptors_list = receptors_list)
        total_list = []
        for obs in obs_dic:
            obs = obs_dic[obs].flatten()
            valid_ind = np.where(~np.isnan(obs))
            obs = obs[valid_ind]
            total_list += obs.tolist()
        total_list = np.array(total_list)
        return self.percentile_range(total_list, p)

    def get_vrange_foot(self, start = None, end = None, dt = None, p = 100, islog = True):    
        if start is None:
            Start = self.Start
        if end is None:
            End = self.End
        if dt is None:
            dt = self.dt
        if isinstance(Start, str):
            Start = dtm.datetime.strptime(Start, "%Y-%m-%d_%H:%M:%S")
        if isinstance(End, str):
            End = dtm.datetime.strptime(End, "%Y-%m-%d_%H:%M:%S")
        if isinstance(dt, int):
            dt = dtm.timedelta(hours = dt)
        total_list = []
        Current = Start
        while(Current <= End):
            opti_time_list = self.opti_recepTime_4eachTime[Current.strftime("%Y%m%d%H")]
            N_reduced = len(opti_time_list)
            for i_reduced in range(N_reduced):
                recepTime = opti_time_list[i_reduced]
                footprint = self.get_footprint(time_list = [Current], RecepTime_list = [recepTime])
                if islog:
                    footprint = np.log10(footprint)
                    valid_ind = np.where((~np.isnan(footprint)) & (~np.isinf(footprint)))
                    footprint  = footprint[valid_ind]
                total_list += footprint.flatten().tolist()
            Current += dt
        total_list = np.array(total_list)
        #set_trace()
        return self.percentile_range(total_list, p)
    
    def get_vrange_dEmiss(self, time = None, RecepTime = None, p = 100, center0 = False):
        dEmiss_list = self.get_dEmiss(time, RecepTime)
        total_list = []
        for dEmiss_dic in dEmiss_list:
            total_list.append(dEmiss_dic["dEmiss"].flatten())
        total_list = np.array(total_list)
        vrange = self.percentile_range(total_list, p)
        if center0:
            vrange = self.center0_adjust(vrange)
        return vrange

    def get_vrange_reduced_unc(self, time = None, RecepTime = None, p = 100, islog = True):

        RQ_list = self.get_reduced_unc(time, RecepTime)
        total_list = []
        for RQ in RQ_list:
            RQ_temp = RQ["RQ"].flatten()
            #set_trace()
            if islog:
                RQ_temp = np.log10(RQ_temp)
                valid_ind = np.where((~np.isnan(RQ_temp)) & (~np.isinf(RQ_temp)))
                RQ_temp = RQ_temp[valid_ind]
            RQ_temp = RQ_temp.tolist()
            total_list += RQ_temp
        total_list = np.array(total_list)
        #vrange = self.percentile_range(total_list, p)
        return self.percentile_range(total_list, p)

    def get_vrange_emiss(self, EmissType = "Posterior", start = None, end = None, dt = None, p = 100, center0 = False):

        EmissType = EmissType.lower()
        assert EmissType in ["prior", "posterior"]
        getFun = self.get_prior_emiss if EmissType == "prior" else self.get_posterior_emiss
        if start is None:
            Start = self.Start
        if end is None:
            End = self.End
        if dt is None:
            dt = self.dt
        if isinstance(Start, str):
            Start = dtm.datetime.strptime(Start, "%Y-%m-%d_%H:%M:%S")
        if isinstance(End, str):
            End = dtm.datetime.strptime(End, "%Y-%m-%d_%H:%M:%S")
        if isinstance(dt, int):
            dt = dtm.timedelta(hours = dt)
        total_list = []
        Current = Start
        while(Current <= End):
            emiss = getFun(Current)    
            total_list.append(emiss.flatten())
            Current += dt
        total_list = np.array(total_list)
        vrange = self.percentile_range(total_list, p)
        if center0:
            vrange = self.center0_adjust(vrange, force = True)
        return vrange
