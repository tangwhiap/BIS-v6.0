#!/usr/bin/env python
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import datetime as dtm
import pandas as pd
from .BIS_Analyst import BIS_Case_Analyst


class REAL_Analyst(BIS_Case_Analyst):

    def __init__(self, *argv, **kwargs):
        BIS_Case_Analyst.__init__(self, *argv, **kwargs)
        self.CaseType = "REAL"

    def get_BCK(self, timeList = None, receptors_list = None, BCKType = "orig"):
        #print("intget_OBS", time)
        BCKType = BCKType.lower()
        assert BCKType in ["orig", "optimized"]
        if timeList is None:
            timeList = self.receptor_times
        if isinstance(timeList, dtm.datetime) or isinstance(timeList, pd.Timestamp):
            timeList = [timeList]
    
        if receptors_list is None:
            receptors_list = self.ObsName_list

        if type(receptors_list) == type("string"):
            receptors_list = [receptors_list]
        #sid = list(map(lambda x: x in receptors_list, sites_list))
        bck_dic = {}
        recep_ind_list = []
        for receptor in receptors_list:
            recep_ind = self.ObsName_list.index(receptor)
            recep_ind_list.append(recep_ind)
            bck_dic[receptor] = []
        recep_ind_list = np.array(recep_ind_list)

        for time in timeList:
        
            ncf_bck = nc.Dataset(self.CaseDir + "/history_save/" + time.strftime("%Y%m%d%H") + "/Qbck.nc", "r")
            #sites_list = ncf_bck.variables["sites"][:]
            sid = 0 if BCKType == "orig" else -1
            bck = ncf_bck.variables["BCK"][sid, recep_ind_list].filled(np.nan)
            ncf_bck.close()
            for irecep, receptor in enumerate(receptors_list):
                bck_dic[receptor].append(bck[irecep])
        for receptor in receptors_list:
            bck_dic[receptor] = np.array(bck_dic[receptor])
            #if receptors_list is not None:
            #    bck = ncf_bck.variables["bck"][sid].filled(np.nan)
            #else:
            #    bck = ncf_bck.variables["bck"][:].filled(np.nan)
            #    receptors_list = ncf_bck.variables["sites"][:]
        return timeList, bck_dic

class REAL_Analyst_Independent(REAL_Analyst):

    def __init__(self, *argv, **kwargs):
        REAL_Analyst.__init__(self, *argv, **kwargs)

