#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from ..REAL_base.main_class import ExpREAL
from ....utils.module2dic import module2dic

#from . import benchplots_configure as benchConfig

import numpy as np
import pandas as pd
from pdb import set_trace

class ExpREAL_independent(ExpREAL):

    def __init__(self, *args, **kwargs):
        self.load_indicators()
        ExpREAL.__init__(self, *args, **kwargs)
        config = kwargs["configure"]
        #assert config["EXPNAME"] == "independent"
        assert len(self.sectors) == 1, "This experiment cannot be compatible with multiple sectors"
        self.myConfig = dict(**self.myConfig, **config["ExpConfig"]["independent"])

        #-- Create observation and background object --#
        self.build_observation_object()
        self.objOBS.get_data_offline(start = self.Start, end = self.End, dt = self.dt)

        self.build_background_object()
    
        #-- Load benchplots configures --#
        #self.benchConfig = module2dic(benchConfig)

        LON, LAT = self.get_LONLAT()
        self.newLON = LON
        self.newLAT = LAT
        self.make_emissSetting()
        self.make_HSetting()
        self.make_ESetting()
        self.make_DSetting()
        self.make_EbckSetting()
        self.make_DbckSetting()
        self.make_RSetting()
        print("Successfully create ExpREAL_independent objcet.")

    def load_indicators(self):
        from .indicators import Sigma_Indicator, X_Indicator, dX_Indicator, dSigma_Indicator
        self.indicators = {"Sigma": Sigma_Indicator, "X": X_Indicator, "dX": dX_Indicator, "dSigma": dSigma_Indicator}
        
    def get_LONLAT(self):

        Htype = list(self.HtypeToFootprint.keys())[0]
        LON, LAT = self.objFootDic[Htype].get_dom_config()

        return LON, LAT

    def make_emissSetting(self):
        self.sectors_PriorSetting = {}
        self.sectors_SigmaSetting = {}
        lon = self.newLON[0, :]
        lat = self.newLAT[:, 0]
        for sector in self.sectors:
            getFun_kwargs_Emiss = {"lon_s": lon[0], "lon_e": lon[-1], "lat_s": lat[0], "lat_e": lat[-1], "LON": self.newLON, "LAT": self.newLAT}
            getFun_kwargs_Prior = {**getFun_kwargs_Emiss, **self.myConfig["PriorKwargs"]}
            getFun_kwargs_Sigma = {
                "funName": self.myConfig["SigmaMethod"], "sectorName": sector, "getFun_prior": self.get_prior_emiss, "LON": self.newLON, "LAT": self.newLAT,
                **self.myConfig["SigmaKwargs"],
            }
            
            self.sectors_PriorSetting[sector] = {
                "emiss": self.myConfig["PriorName"],
                "getFun_kwargs": getFun_kwargs_Prior,
                "LON": self.newLON,
                "LAT": self.newLAT,
            }
            self.sectors_SigmaSetting[sector] = {
                "getSigmaFun": self.compute_sigma,
                "getFun_kwargs": getFun_kwargs_Sigma,
                "isRegrid": False,
            }

    def make_HSetting(self):
        self.types_HSetting = {}
        for type_ in self.typeToSector_H:
            self.types_HSetting[type_] = {
                "backtime_j2i": self.objCore.backtime_j2i,
                "recepTimeList": self.objCore.objIter.obsTimeList,
            }

    def make_ESetting(self):
        self.types_ESetting = {}
        for type_ in self.typeToSector_E:
            self.types_ESetting[type_] = {
                "method": self.myConfig["EMethod"],
                "method_kwargs": {"Ls": self.myConfig["EKwargs"]["Ls"], "LON": self.newLON, "LAT": self.newLAT},
            }

    def make_DSetting(self):

        DTList = pd.date_range(self.Start, self.End, freq = str(self.dtHrs) + "H")
        dim1TList = dim2TList = [time.to_pydatetime() for time in DTList]
        self.types_DSetting = {}
        for type_ in self.typeToSector_D:
            self.types_DSetting[type_] = {
                "method": self.myConfig["DMethod"],
                "method_kwargs": {"Lt": self.myConfig["DKwargs"]["Lt"], "dim1TList": dim1TList, "dim2TList": dim2TList},
            }

    def make_EbckSetting(self):
        receptors = self.objFootDic[list(self.typeToSector_H.keys())[0]].receptors
        lonList = []
        latList = []
        for receptor in receptors:
            lonList.append(float(receptors[receptor]["lon"]))
            latList.append(float(receptors[receptor]["lat"]))
        lonList = np.array(lonList)
        latList = np.array(latList)
        self.EbckSetting = {
            "method": self.myConfig["EbckMethod"],
            "method_kwargs": {"LsBck": self.myConfig["EbckKwargs"]["LsBck"], "lonList": lonList, "latList": latList},
        }

    def make_DbckSetting(self):
        dim1TList = dim2TList = self.objCore.objIter.obsTimeList
        self.DbckSetting = {
            "method": self.myConfig["DbckMethod"],
            "method_kwargs": {"LtBck": self.myConfig["DbckKwargs"]["LtBck"], "dim1TList": dim1TList, "dim2TList": dim2TList},
        }

    def make_RSetting(self):
        #RMethod = self.myConfig["RMethod"]
        #error_inReceptors = self.myConfig["RKwargs"]["error_inReceptors"]
        self.RSetting = {}
        self.RSetting = {
            "method": self.myConfig["RMethod"], 
            "method_kwargs": self.myConfig["RKwargs"],
        }
