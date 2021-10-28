#!/usr/bin/env python

# Authors:
#    Wenhan TANG - 08/2021
#    ...

from ..independent.main_class import ExpREAL_independent
from ....utils.netcdf_io import get_ncvar
from ....utils.distance import areaS

import numpy as np
import multiprocessing as mtp


class ExpREAL_independent_priorFnet(ExpREAL_independent):

    def __init__(self, *args, **kwargs):
        ExpREAL_independent.__init__(self, *args, **kwargs)
        config = kwargs["configure"]
        self.myConfig = dict(**self.myConfig, **config["ExpConfig"]["independent_priorFnet"])
        self.directoriesToMake += ["ffeDir", "cftaDir", "ffeSigmaDir", "cftaSigmaDir"]
        print("Successfully create ExpREAL_independent_priorFnet objcet.")

    def get_ffePriorFile_name(self, time, sectorName):
        return self.myConfig["ffeDir"] + "/" + self.myConfig["ffe_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_cftaPriorFile_name(self, time, sectorName):
        return self.myConfig["cftaDir"] + "/" + self.myConfig["cfta_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_ffeSigmaPriorFile_name(self, time, sectorName):
        return self.myConfig["ffeSigmaDir"] + "/" + self.myConfig["ffeSigma_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_cftaSigmaPriorFile_name(self, time, sectorName):
        return self.myConfig["cftaSigmaDir"] + "/" + self.myConfig["cftaSigma_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"
    
    def get_ffeHourlyFile_name(self, time, sectorName):
        return self.myConfig["outHourlyDir"] + "/" + self.myConfig["ffeHourly_Prefix"] + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_sffeHourlyFile_name(self, time, sectorName):
        return self.myConfig["outHourlyDir"] + "/" + self.myConfig["sffeHourly_Prefix"] + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_cftaHourlyFile_name(self, time, sectorName):
        return self.myConfig["outHourlyDir"] + "/" + self.myConfig["cftaHourly_Prefix"] + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_scftaHourlyFile_name(self, time, sectorName):
        return self.myConfig["outHourlyDir"] + "/" + self.myConfig["scftaHourly_Prefix"] + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_ffeDailyFile_name(self, time, sectorName):
        return self.myConfig["outDailyDir"] + "/" + self.myConfig["ffeDaily_Prefix"] + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_sffeDailyFile_name(self, time, sectorName):
        return self.myConfig["outDailyDir"] + "/" + self.myConfig["sffeDaily_Prefix"] + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_cftaDailyFile_name(self, time, sectorName):
        return self.myConfig["outDailyDir"] + "/" + self.myConfig["cftaDaily_Prefix"] + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_scftaDailyFile_name(self, time, sectorName):
        return self.myConfig["outDailyDir"] + "/" + self.myConfig["scftaDaily_Prefix"] + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_ffeWeeklyFile_name(self, time, sectorName):
        return self.myConfig["outWeeklyDir"] + "/" + self.myConfig["ffeWeekly_Prefix"] + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_sffeWeeklyFile_name(self, time, sectorName):
        return self.myConfig["outWeeklyDir"] + "/" + self.myConfig["sffeWeekly_Prefix"] + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_cftaWeeklyFile_name(self, time, sectorName):
        return self.myConfig["outWeeklyDir"] + "/" + self.myConfig["cftaWeekly_Prefix"] + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_scftaWeeklyFile_name(self, time, sectorName):
        return self.myConfig["outWeeklyDir"] + "/" + self.myConfig["scftaWeekly_Prefix"] + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_ffeMonthlyFile_name(self, time, sectorName):
        return self.myConfig["outMonthlyDir"] + "/" + self.myConfig["ffeMonthly_Prefix"] + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_sffeMonthlyFile_name(self, time, sectorName):
        return self.myConfig["outMonthlyDir"] + "/" + self.myConfig["sffeMonthly_Prefix"] + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_cftaMonthlyFile_name(self, time, sectorName):
        return self.myConfig["outMonthlyDir"] + "/" + self.myConfig["cftaMonthly_Prefix"] + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_scftaMonthlyFile_name(self, time, sectorName):
        return self.myConfig["outMonthlyDir"] + "/" + self.myConfig["scftaMonthly_Prefix"] + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_ffeAllFile_name(self, sectorName):
        return self.myConfig["outAllDir"] + "/" + self.myConfig["ffeAll_Prefix"] + ".nc"

    def get_sffeAllFile_name(self, sectorName):
        return self.myConfig["outAllDir"] + "/" + self.myConfig["sffeAll_Prefix"] + ".nc"

    def get_cftaAllFile_name(self, sectorName):
        return self.myConfig["outAllDir"] + "/" + self.myConfig["cftaAll_Prefix"] + ".nc"

    def get_scftaAllFile_name(self, sectorName):
        return self.myConfig["outAllDir"] + "/" + self.myConfig["scftaAll_Prefix"] + ".nc"


    def emiss_file_add(self, args):
        time, sectorName = args
        ffePriorFile = self.get_ffePriorFile_name(time, sectorName)
        cftaPriorFile = self.get_cftaPriorFile_name(time, sectorName)
        [ffe_data, ffe_emiss, ffe_LON, ffe_LAT] = get_ncvar(ffePriorFile, ["data", "emiss", "LON", "LAT"])
        [cfta_data, cfta_emiss, cfta_LON, cfta_LAT] = get_ncvar(cftaPriorFile, ["data", "emiss", "LON", "LAT"])
        assert np.sum(np.abs(ffe_LON - cfta_LON)) <= 1e-6
        assert np.sum(np.abs(ffe_LAT - cfta_LAT)) <= 1e-6
        emiss = ffe_emiss + cfta_emiss
        LON = ffe_LON; LAT = ffe_LAT
        priorFile = self.get_priorFile_name(time, sectorName)
        self.X_Sigma_Indicator_write(priorFile, emiss, LON, LAT)

    def sigma_file_add(self, args):
        time, sectorName = args
        ffeSigmaPriorFile = self.get_ffeSigmaPriorFile_name(time, sectorName)
        cftaSigmaPriorFile = self.get_cftaSigmaPriorFile_name(time, sectorName)
        [ffe_data, ffe_emiss, ffe_LON, ffe_LAT] = get_ncvar(ffeSigmaPriorFile, ["data", "emiss", "LON", "LAT"])
        [cfta_data, cfta_emiss, cfta_LON, cfta_LAT] = get_ncvar(cftaSigmaPriorFile, ["data", "emiss", "LON", "LAT"])
        assert np.sum(np.abs(ffe_LON - cfta_LON)) <= 1e-6
        assert np.sum(np.abs(ffe_LAT - cfta_LAT)) <= 1e-6
        emiss = np.sqrt(ffe_emiss **2 + cfta_emiss **2)
        LON = ffe_LON; LAT = ffe_LAT
        priorFile = self.get_sigmaPriorFile_name(time, sectorName)
        self.X_Sigma_Indicator_write(priorFile, emiss, LON, LAT)

    def move_ffe_cfta(self, ffe_or_cfta, emiss_or_sigma):

        assert ffe_or_cfta.lower() in ["ffe", "cfta"]
        assert emiss_or_sigma.lower() in ["emiss", "sigma"]

        if emiss_or_sigma.lower() == "emiss":
            get_File_from = self.get_priorFile_name
        
            if ffe_or_cfta.lower() == "ffe":
                get_File_to = self.get_ffePriorFile_name

            if ffe_or_cfta.lower() == "cfta":
                get_File_to = self.get_cftaPriorFile_name
            
        elif emiss_or_sigma.lower() == "sigma":
            get_File_from = self.get_sigmaPriorFile_name

            if ffe_or_cfta.lower() == "ffe":
                get_File_to = self.get_ffeSigmaPriorFile_name

            if ffe_or_cfta.lower() == "cfta":
                get_File_to = self.get_cftaSigmaPriorFile_name

        current = self.Start
        while(current <= self.End):
            for sector in self.sectors:
                kwargs = {"time": current, "sectorName": sector}
                self.move_files(get_File_from, get_File_to, fun_from_kwargs = kwargs, fun_to_kwargs = kwargs)

            current += self.dt   

    def get_ffe_emiss(self, time, sectorName, outType = "orig"):
        assert outType in ["vector", "orig"]
        if outType == "orig":
            emiss, LON, LAT = get_ncvar(self.get_ffePriorFile_name(time, sectorName), ["emiss", "LON", "LAT"])
            return emiss, LON, LAT
        if outType == "vector":
            data = get_ncvar(self.get_ffePriorFile_name(time, sectorName), ["data"])
            return data
    
    def get_ffe_sigma(self, time, sectorName, outType = "orig"):
        assert outType in ["vector", "orig"]
        if outType == "orig":
            emiss, LON, LAT = get_ncvar(self.get_ffeSigmaPriorFile_name(time, sectorName), ["emiss", "LON", "LAT"])
            return emiss, LON, LAT
        if outType == "vector":
            data = get_ncvar(self.get_ffeSigmaPriorFile_name(time, sectorName), ["data"])
            return data

    def get_cfta_emiss(self, time, sectorName, outType = "orig"):
        assert outType in ["vector", "orig"]
        if outType == "orig":
            emiss, LON, LAT = get_ncvar(self.get_cftaPriorFile_name(time, sectorName), ["emiss", "LON", "LAT"])
            return emiss, LON, LAT
        if outType == "vector":
            data = get_ncvar(self.get_cftaPriorFile_name(time, sectorName), ["data"])
            return data

    def get_cfta_sigma(self, time, sectorName, outType = "orig"):
        assert outType in ["vector", "orig"]
        if outType == "orig":
            emiss, LON, LAT = get_ncvar(self.get_cftaSigmaPriorFile_name(time, sectorName), ["emiss", "LON", "LAT"])
            return emiss, LON, LAT
        if outType == "vector":
            data = get_ncvar(self.get_cftaSigmaPriorFile_name(time, sectorName), ["data"])
            return data


    def make_prior_emiss(self, sectors_PriorSetting):
        lon = self.newLON[0, :]
        lat = self.newLAT[:, 0]
        getFun_kwargs_Emiss = {"lon_s": lon[0], "lon_e": lon[-1], "lat_s": lat[0], "lat_e": lat[-1], "LON": self.newLON, "LAT": self.newLAT}
        for sector in self.sectors:
            sectors_PriorSetting[sector]["emiss"] = self.myConfig["ffeName"]
            sectors_PriorSetting[sector]["getFun_kwargs"] = {**getFun_kwargs_Emiss, **self.myConfig["ffeKwargs"]}
        ExpREAL_independent.make_prior_emiss(self, sectors_PriorSetting)    
        self.move_ffe_cfta(ffe_or_cfta = "ffe", emiss_or_sigma = "emiss")

        for sector in self.sectors:
            sectors_PriorSetting[sector]["emiss"] = self.myConfig["cftaName"]
            sectors_PriorSetting[sector]["getFun_kwargs"] = {**getFun_kwargs_Emiss, **self.myConfig["cftaKwargs"]}
        ExpREAL_independent.make_prior_emiss(self, sectors_PriorSetting)
        self.move_ffe_cfta(ffe_or_cfta = "cfta", emiss_or_sigma = "emiss")

        parallelArgs = []
        current = self.Start
        while(current <= self.End):
            for sector in self.sectors:
                #self.emiss_file_add(current, sectorName = sector)
                parallelArgs.append((current, sector))
            current += self.dt

        pool = mtp.Pool(self.nProc)
        pool.map(self.emiss_file_add, parallelArgs)
        pool.close()
        pool.join()

    def make_prior_sigma(self, sectors_SigmaSetting):

        for sector in self.sectors:
            getFun_kwargs_Sigma = {
                    "funName": self.myConfig["ffeSigmaMethod"], "sectorName": sector, "getFun_prior": self.get_ffe_emiss, "LON": self.newLON, "LAT": self.newLAT,
                    **self.myConfig["ffeSigmaKwargs"],
                }
            sectors_SigmaSetting[sector]["getFun_kwargs"] = getFun_kwargs_Sigma
        ExpREAL_independent.make_prior_sigma(self, sectors_SigmaSetting)
        self.move_ffe_cfta(ffe_or_cfta = "ffe", emiss_or_sigma = "sigma")

        for sector in self.sectors:
            getFun_kwargs_Sigma = {
                    "funName": self.myConfig["cftaSigmaMethod"], "sectorName": sector, "getFun_prior": self.get_cfta_emiss, "LON": self.newLON, "LAT": self.newLAT,
                    **self.myConfig["cftaSigmaKwargs"],
                }
            sectors_SigmaSetting[sector]["getFun_kwargs"] = getFun_kwargs_Sigma
        ExpREAL_independent.make_prior_sigma(self, sectors_SigmaSetting)
        self.move_ffe_cfta(ffe_or_cfta = "cfta", emiss_or_sigma = "sigma")

        parallelArgs = []
        current = self.Start
        while(current <= self.End):
            for sector in self.sectors:
                #self.emiss_file_add(current, sectorName = sector)
                parallelArgs.append((current, sector))
            current += self.dt

        pool = mtp.Pool(self.nProc)
        pool.map(self.sigma_file_add, parallelArgs)
        pool.close()
        pool.join()

    
    def build_output_dic(self):

        ExpREAL_independent.build_output_dic(self)

        self.make_output_dic["ffeEmiss"] = {
            "inputFun": self.get_ffe_emiss,
            "outNameFun_hourly": self.get_ffeHourlyFile_name,
            "outNameFun_daily": self.get_ffeDailyFile_name,
            "outNameFun_weekly": self.get_ffeWeeklyFile_name,
            "outNameFun_monthly": self.get_ffeMonthlyFile_name,
            "outNameFun_all": self.get_ffeAllFile_name,
            #--- umol/m^2/s => tCO2/cell/hour --#
            "transFun": trans_emiss_unit,
            "outFun": straight,
        }

        self.make_output_dic["cftaEmiss"] = {
            "inputFun": self.get_cfta_emiss,
            "outNameFun_hourly": self.get_cftaHourlyFile_name,
            "outNameFun_daily": self.get_cftaDailyFile_name,
            "outNameFun_weekly": self.get_cftaWeeklyFile_name,
            "outNameFun_monthly": self.get_cftaMonthlyFile_name,
            "outNameFun_all": self.get_cftaAllFile_name,
            #--- umol/m^2/s => tCO2/cell/hour --#
            "transFun": trans_emiss_unit,
            "outFun": straight,
        }

        self.make_output_dic["ffeSigma"] = {
            "inputFun": self.get_ffe_sigma,
            "outNameFun_hourly": self.get_sffeHourlyFile_name,
            "outNameFun_daily": self.get_sffeDailyFile_name,
            "outNameFun_weekly": self.get_sffeWeeklyFile_name,
            "outNameFun_monthly": self.get_sffeMonthlyFile_name,
            "outNameFun_all": self.get_sffeAllFile_name,
            #--- convert sigma to square ---#
            #--- umol/m^2/s => (tCO2/cell/hour)^2 --#
            "transFun": trans_sigma_unit,
            #--- convert square to sigma ---#
            "outFun": np.sqrt
        }

        self.make_output_dic["cftaSigma"] = {
            "inputFun": self.get_cfta_sigma,
            "outNameFun_hourly": self.get_scftaHourlyFile_name,
            "outNameFun_daily": self.get_scftaDailyFile_name,
            "outNameFun_weekly": self.get_scftaWeeklyFile_name,
            "outNameFun_monthly": self.get_scftaMonthlyFile_name,
            "outNameFun_all": self.get_scftaAllFile_name,
            #--- convert sigma to square ---#
            #--- umol/m^2/s => (tCO2/cell/hour)^2 --#
            "transFun": trans_sigma_unit,
            #--- convert square to sigma ---#
            "outFun": np.sqrt
        }



def trans_emiss_unit(arr, LON, LAT):
    return arr * areaS(LON, LAT) * 3600 / 1.0e6 * 44 / 1.0e6

def trans_sigma_unit(arr, LON, LAT):
    return (arr * areaS(LON, LAT) * 3600 / 1.0e6 * 44 / 1.0e6) ** 2

def straight(arr, *args, **kwargs):
    return arr
