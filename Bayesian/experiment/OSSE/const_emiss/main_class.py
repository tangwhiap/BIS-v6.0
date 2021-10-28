#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from ..OSSE_base.main_class import ExpOSSE
from ....utils.netcdf_io import get_ncvar
from ....utils.distance import areaS

import numpy as np
import pandas as pd
import multiprocessing as mtp
from pdb import set_trace

myName = "const_emiss"

class ExpOSSE_const_emiss(ExpOSSE):

    """
        *** Only for SINGLE sector ***
        *** Only for normal lat-lon grid ***
    """

    def __init__(self, *args, **kwargs):
        
        self.load_indicators()
        ExpOSSE.__init__(self, *args, **kwargs)
        config = kwargs["configure"]
        assert config["EXPNAME"] == myName or config["BASEEXP"] == myName
        assert len(self.sectors) == 1, "This experiment cannot be compatible with multiple sectors"
        self.myConfig = dict(**self.myConfig, **config["ExpConfig"][myName])
        LON, LAT = self.get_LONLAT()
        self.newLON = LON
        self.newLAT = LAT
        self.make_emissSetting()
        self.make_HSetting()
        self.make_ESetting()
        self.make_DSetting()
        self.make_RSetting()
        print("Successfully create ExpOSSE_const_emiss objcet.")


    def load_indicators(self):
        from .indicators import Sigma_Indicator, X_Indicator, dX_Indicator, dSigma_Indicator
        self.indicators = {"Sigma": Sigma_Indicator, "X": X_Indicator, "dX": dX_Indicator, "dSigma": dSigma_Indicator}

    def get_LONLAT(self):

        Htype = list(self.HtypeToFootprint.keys())[0]
        LON, LAT = self.objFootDic[Htype].get_dom_config()

        return LON, LAT

    def make_emissSetting(self):
        self.sectors_PriorSetting = {}
        self.sectors_TruthSetting = {}
        self.sectors_SigmaSetting = {}
        lon = self.newLON[0, :]
        lat = self.newLAT[:, 0]
        for sector in self.sectors:
            getFun_kwargs_Emiss = {"lon_s": lon[0], "lon_e": lon[-1], "lat_s": lat[0], "lat_e": lat[-1], "LON": self.newLON, "LAT": self.newLAT}
            #getFun_kwargs_Prior = {**getFun_kwargs_Emiss, **self.myConfig["PriorKwargs"]}
            getFun_kwargs_Truth = {**getFun_kwargs_Emiss, **self.myConfig["TruthKwargs"]}

            #getFun_kwargs_Prior = {**getFun_kwargs_Emiss, **self.myConfig["PriorKwargs"]}
            getFun_kwargs_Sigma = {
                "funName": self.myConfig["SigmaMethod"], "sectorName": sector, "getFun_prior": self.get_prior_emiss, "getFun_truth": self.get_truth_emiss, "LON": self.newLON, "LAT": self.newLAT,
                "Start": self.Start, "End": self.End, "dt": self.dt, **self.myConfig["SigmaKwargs"],
            }
            
            #self.sectors_PriorSetting[sector] = {
            #    "emiss": self.myConfig["PriorName"],
            #    "getFun_kwargs": getFun_kwargs_Prior,
            #    "LON": self.newLON,
            #    "LAT": self.newLAT,
            #}

            self.sectors_TruthSetting[sector] = {
                "emiss": self.myConfig["TruthName"],
                "getFun_kwargs": getFun_kwargs_Truth,
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
    def make_RSetting(self):
        #RMethod = self.myConfig["RMethod"]
        #error_inReceptors = self.myConfig["RKwargs"]["error_inReceptors"]
        self.RSetting = {}
        self.RSetting = {
            "method": self.myConfig["RMethod"], 
            "method_kwargs": self.myConfig["RKwargs"],
        }

    def make_prior_emiss(self, *args, **kwargs):
        make_prior_emiss(self, *args, **kwargs)

def make_prior_emiss(self, sectors_PriorSetting, start = None, end = None, dt = None, regridFun = None):
    print("make_prior!!!")

    #def areaS(LON, LAT):
    #    R = 6371000
    #    dlon = (LON[0, 1:] - LON[0, :-1]).mean()
    #    dlat = (LAT[1:, 0] - LAT[:-1, 0]).mean()
    #    S = (R**2) * np.deg2rad(dlon) * (np.sin(np.deg2rad(LAT + dlat/2 )) - np.sin(np.deg2rad(LAT - dlat/2 )))
    #    assert S.shape == LON.shape
    #    assert S.shape == LAT.shape
    #    return S

    if start is None:
        start = self.Start
    if end is None:
        end = self.End
    if dt is None:
        dt = self.dt

    area = areaS(self.newLON, self.newLAT)
    for sectorName in self.sectors:
        emiss_total = 0
        current = start
        icount = 0
        while(current <= end):
            emiss = get_ncvar(self.get_truthFile_name(current, sectorName), "emiss")
            assert emiss.shape == self.newLON.shape
            emiss_total = emiss_total + (emiss * area).sum()
            icount += 1
            current += dt
        emiss_ave = emiss_total / icount / (area.sum())

        from ....emiss.constant import interface

        ParallelArgs = []
        current = start
        while(current <= end):
            outFile = self.get_priorFile_name(current, sectorName)
            processEmissKwargs = {
                "Current": current, "getDataFun": interface, "getDataFun_kwargs": {"LON": self.newLON, "LAT": self.newLAT, "const": emiss_ave}, 
                "OutFile": outFile, "writeFun": self.X_Sigma_Indicator_write, "isRegrid": False
            }
            self.process_emiss(processEmissKwargs)
            ParallelArgs.append(processEmissKwargs)
            current += dt


        #pool = mtp.Pool(self.nProc)
        #pool.map(self.process_emiss, ParallelArgs)
        #pool.close()
        #pool.join()
        


