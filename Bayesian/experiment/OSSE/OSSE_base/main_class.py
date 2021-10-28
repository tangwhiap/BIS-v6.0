#!/usr/bin/env python

#-> Base class for OSSE experiment object.<-#

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from ...exp_base.main_class import ExpBIS, trans_emiss_unit, trans_sigma_unit, straight
from ....utils.netcdf_io import obs_write, get_ncvar
from ....utils.module2dic import module2dic

#from . import benchplots_configure as benchConfig

import numpy as np
import multiprocessing as mtp
from pdb import set_trace

class ExpOSSE(ExpBIS):

    """
        The base class of OSSE experiment, derivated from ExpBIS class.

        Add methods:

            * get_truthFile_name: Get truth emission files name.

            * make_truth_emiss: Processing and writing the truth emission.

            * get_truth_emiss: Get truth emission data from their files.

            * make_obs -> make_simobs:  Computing observation data, calling the make_simobs method, which can simulate the observation by truth emission and footprint (H matrix)

            * exp_initialize: Initialization of this OSSE experiment.

            * exp_run: The main function without any arguments.

    """

    def __init__(self, *args, **kwargs):

        self.isOSSE = True
        self.isREAL = False

        self.hasBCK = False

        ExpBIS.__init__(self, *args, **kwargs)
        self.myConfig = dict(**self.myConfig, **kwargs["configure"]["ExpConfig"]["OSSE"])
        self.directoriesToMake += ["emissTruthDir"]

        #self.benchConfig = module2dic(benchConfig)

    def get_truthFile_name(self, time, sectorName):
        return self.myConfig["emissTruthDir"] + "/" + self.myConfig["emissTruth_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_truthHourlyFile_name(self, time, sectorName):
        return self.myConfig["outHourlyDir"] + "/" + self.myConfig["truthHourly_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_truthDailyFile_name(self, time, sectorName):
        return self.myConfig["outDailyDir"] + "/" + self.myConfig["truthDaily_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_truthWeeklyFile_name(self, time, sectorName):
        return self.myConfig["outWeeklyDir"] + "/" + self.myConfig["truthWeekly_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_truthMonthlyFile_name(self, time, sectorName):
        return self.myConfig["outMonthlyDir"] + "/" + self.myConfig["truthMonthly_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_truthAllFile_name(self, sectorName):
        return self.myConfig["outAllDir"] + "/" + self.myConfig["truthAll_Prefix"] + "_" + sectorName + ".nc"

    def make_truth_emiss(self, sectors_EmissSetting, start = None, end = None, dt = None, regridFun = None):

        if start is None:
            start = self.Start
        if end is None:
            end = self.End
        if dt is None:
            dt = self.dt
        
        for sectorName in sectors_EmissSetting:
            EmissSetting = sectors_EmissSetting[sectorName]
            emissName = EmissSetting["emiss"]

            if "LON" in EmissSetting:
                newLON = EmissSetting["LON"]
            else:
                newLON = None

            if "LAT" in EmissSetting:
                newLAT = EmissSetting["LAT"]
            else:
                newLAT = None

            if "getFun_kwargs" in EmissSetting:
                getFun_kwargs = EmissSetting["getFun_kwargs"]
            else:
                getFun_kwargs = {}

            getFun_kwargs["sector"] = sectorName

            loc = locals()
            exec("from ....emiss." + emissName + " import interface, hasInterpolated, interpolate_method")
            getDataFun = loc["interface"]
            hasInterpolated = loc["hasInterpolated"]
            interpolate_method = loc["interpolate_method"]
            isRegrid = not(hasInterpolated)
            #getDataFun = lambda time: interface(time = time, **getFun_kwargs)
            #def getDataFun(time):
            #    return interface(time = time, **getFun_kwargs)
            parallelArgs = []
            current = start
            while(current <= end):
                #print(current)
                outFile = self.myConfig["emissTruthDir"] + "/" + self.myConfig["emissTruth_Prefix"] + "_" + sectorName + "_" + current.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"
                # Arguments of function process_emiss:
                # Current, getDataFun, OutFile, writeFun, isRegrid, regridMethod, regridFun, newLON, newLAT 
                processEmissKwargs = {
                    "Current": current, "getDataFun": getDataFun, "getDataFun_kwargs": getFun_kwargs, "OutFile": outFile, "writeFun": self.X_Sigma_Indicator_write,
                    "isRegrid": isRegrid, "regridMethod": interpolate_method, "regridFun": regridFun, "newLON": newLON, "newLAT": newLAT,
                }
                parallelArgs.append(processEmissKwargs)
                #self.process_emiss(processEmissKwargs)
                current += dt

        pool = mtp.Pool(self.nProc)
        pool.map(self.process_emiss, parallelArgs)
        pool.close()
        pool.join()

    def get_truth_emiss(self, time, sectorName, outType = "orig"):
        assert outType in ["vector", "orig"]
        if outType == "orig":
            emiss, LON, LAT = get_ncvar(self.get_truthFile_name(time, sectorName), ["emiss", "LON", "LAT"])
            return emiss, LON, LAT
        if outType == "vector":
            data = get_ncvar(self.get_truthFile_name(time, sectorName), ["data"])
            return data

    def make_obs(self):
        obsTimeList = self.objCore.objIter.obsTimeList
        #ParallelArgs = []
        #for obsTime in obsTimeList:
        #    print(obsTime)
        #    self.make_simobs(recepTime = obsTime)
        pool = mtp.Pool(self.nProc)
        pool.map(self.make_simobs, obsTimeList)
        pool.close()
        pool.join()

    def make_simobs(self, recepTime):
        #print("make obs:", recepTime)
        backtime_j2i = self.objCore.backtime_j2i
        sumHx = self.objCore.Hx(timeList = backtime_j2i(recepTime), Xtype = "Truth", recepTime = recepTime, class_XIndicator = self.indicators["X"], objExp = self)
        obs = self.objCore.Hx2obs(sumHx = sumHx, method = "sum")
        recepError = self.myConfig["RKwargs"]["error_inReceptors"]
        errorList = [recepError[receptor] for receptor in recepError]
        randomError = [np.random.normal(0, error) for error in errorList]
        randomError = np.array(randomError).reshape(len(randomError), 1)
        obs = obs + randomError
        obs = np.where(obs >= 0, obs, 0)
        obs_write(outFile = self.get_obsFile_name(recepTime = recepTime), obs = obs)

    def exp_initialize(self, *args, **kwargs):

        """
            The following things should be done in initialization:
            1. Create directories of this case.
            2. Make truth emission files.
            3. Make prior emission files.
            4. Make prior emission uncertainty files.
            5. Copy the prior emission files to their "proc" files.
            6. Same for prior emission uncertainty files.
            7. Make H files.
            8. Make observation data files.
            9. Make E file.
            10. Make D file.
            11. Make R file.
        """

        self.case_mkdir()
        self.make_truth_emiss(self.sectors_TruthSetting)
        self.make_prior_emiss(self.sectors_PriorSetting)
        self.make_prior_sigma(self.sectors_SigmaSetting)
        self.copy_emiss_sigma("emiss")
        self.copy_emiss_sigma("sigma")
        self.make_H_offline(self.types_HSetting, self.objFootDic)
        self.make_obs()
        self.make_E(self.types_ESetting)
        self.make_D(self.types_DSetting)
        if self.hasBCK:
            self.make_Ebck(self.EbckSetting)
            self.make_Dbck(self.DbckSetting)
        self.make_R_offline(self.RSetting)
        self.compute_init_concentration()

    def build_output_dic(self):
        outputType = self.myConfig["outputType"].lower()
        assert outputType in ["sum", "mean"] 
        isSum = outputType == "sum"
        isMean = outputType == "mean"
        ExpBIS.build_output_dic(self)
        self.make_output_dic["truthEmiss"] = {
            "inputFun": self.get_truth_emiss,
            "outNameFun_hourly": self.get_truthHourlyFile_name,
            "outNameFun_daily": self.get_truthDailyFile_name,
            "outNameFun_weekly": self.get_truthWeeklyFile_name,
            "outNameFun_monthly": self.get_truthMonthlyFile_name,
            "outNameFun_all": self.get_truthAllFile_name,
            #--- umol/m^2/s => tCO2/cell/hour --#
            "transFun": trans_emiss_unit if isSum else straight,
            "outFun": straight,
        }
    """
    def exp_run(self):

        #-- Initialization ---#
        self.exp_initialize()

        #-- Main loop of optimization --#
        for recepTime in self.objCore.objIter:

            #-- Calling optimize method of optimization core --#
            print("Optimize: ", recepTime)
            self.objCore.optimize(recepTime = recepTime, objExp = self)

    """
    
        
        
