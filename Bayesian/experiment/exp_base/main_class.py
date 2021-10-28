#!/usr/bin/env python

#-> Bayesian Inversion System EXPERIMENT class <-#

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from .compute_sigma import compute_sigma
from .compute_E import compute_E
from .compute_D import compute_D
from .compute_R import compute_R
from .compute_Ebck import compute_Ebck
from .compute_Dbck import compute_Dbck
from .concentration import compute_init_concentration, compute_final_concentration
from .make_output import make_output

from .indicators import Sigma_Indicator, X_Indicator
from ...utils.netcdf_io import X_Sigma_Indicator_write, get_ncvar
from ...utils.sparse_matrix import array_to_sparse_to_nc
from ...utils.distance import areaS


import numpy as np
import netCDF4 as nc
from scipy import interpolate
import datetime as dtm
import os
from pdb import set_trace
import multiprocessing as mtp

class ExpBIS(object):

    """

        The object of this class represents an experiment of Bayesian Inversion system.

        The main function of this object is dealing the complex and diverse inversion case into a uniform form of
        "H", "X", "E", "D", "Sigma" ... vectors or matrices, combined with multiple sectors or types,
        which can be used by optimization core to perform inversion with the common bayesian inversion method.

        The "H", "X", "E", "D", "Sigma", "dX", "dSigma" here is called "indicator".
        The derivative class of "X", "dX", "Sigma", "dSigma" have been defined in module "indicators" (indicators.py)


        * The attribute of this class:

            config: module of "configure" (configure.py)

            myConfig: configure of ExpConfig "BIS".

            footConfig: Footprint configure list for each type of H.

            Start: Start time of experiment.

            End: End time of experiment.

            dtHrs: Time interval of emission data. (int, unit: hour)

            dt: Time interval of emission data. (datetime.timedelta)

            nProc: CPU core number to use.

            coreName: Name of optimization core.

            sectors: Namelist of emission sectors.

            typeToSector_H: Types of H and their correspondent sector.

            typeToSector_E: Types of E and their correspondent sector.

            typeToSector_D: Types of D and their correspondent sector.

            directoriesToMake: Directories in case folder.

            objCore: Object of optimization core.

            objFootDic: Dictionary of objects of footprint for each type of H.
            
        * The method of this class:
            __init__: Constructor function.

            case_mkdir: Create the case directories defined in self.directoriesToMake.

            exp_initialize: Initialization of experiment. (not defined)

            exp_postproc: Post processing of experiment.

            get_<type>File_name: Function of getting the directory and name of specific data file.
                <type> can be "prior", "proc", "sigmaPrior", "sigmaProc", "H", "obs", "D", "E", "R"

            make funcions: (make_prior_emiss, make_prior_sigma, make_H_offline, make_obs, make_E, make_D, make_R_offline), used to create files of prior emission and uncertainty, footprint, observation, spatial and temporal correlation of emission, and R matrix.

            get_prior_emiss: Function of getting prior emission data.

            get_prior_sigma: Function of getting prior uncertainty data.
            
            copy_files / copy_emiss_sigma / copy_emiss_final: Used to copy files (prior emiss -> proc emiss, prior sigma -> proc sigma, proc emiss -> final emiss)

            X_Sigma_Indicator_write: Alias to the same function in Bayesian.utils.netcdf_io, used to write emission and sigma to nc files.

            process_emiss: Emission processing function, used for regridding and writing emission and sigma data.

            comptute_sigma: Uncertainty computing function, details in module "compute_sigma".

        * The method must (may) be defined in derivative class:

            #exp_run: (necessary) The main function of running the optimization algorithm of this experiment.
            (This function has been defined in ExpBIS class)

            exp_initialize: (optional) The initialization of this experiment. (calling the functions like make_prior_emiss, make_prior_sigma, make_H_offline etc.)

            exp_postproc: (optional) The post processing of this experiment. (like drawing benchplots...)

            make_obs: (necessary) The function of processing or computing observation data.
                    (P.S. The main difference between OSSE and REAL experiment is the definition of this function)

    """

    def __init__(self, *args, **kwargs):
        expbis_constructor(self, *args, **kwargs)


    def case_mkdir(self, *args, **kwargs):

        for dirName in self.directoriesToMake:
            dir_ = self.myConfig[dirName]
            if not os.path.exists(dir_):
                os.makedirs(dir_)
            else:
                print("Directory: " + dir_ + " already exists.")

    def exp_initialize(self, *args, **kwargs):
        """
           * Things should be done before running the optimization core.
           * Something like:
                Creating the subdirectories for the case folder.
                Creating and processing the prior emission.
                Computing the uncertainty of prior emission.
                Creating and processing the truth emission. (for OSSE experiment)
                Creating the observation operator (i.e. H matrix). (for offline mode)
                Preparing the observation data. (for offline mode)
        """
        pass

    def exp_postproc(self, *args, **kwargs):

        """
            * Things should be done after running the optimization core.
            * Something like:
                Copying or linking processing emission files to posterior emission files.
                BenchPlots.
                ...
        """
        self.copy_emiss_final(copy_or_link = self.myConfig["postEmissCoL"])
        self.compute_final_concentration()
        self.make_output()

    def get_priorFile_name(self, time, sectorName):
        return self.myConfig["emissPriorDir"] + "/" + self.myConfig["emissPrior_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_procFile_name(self, time, sectorName):
        return self.myConfig["emissProcDir"] + "/" + self.myConfig["emissProc_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_postFile_name(self, time, sectorName):
        return self.myConfig["emissPostDir"] + "/" + self.myConfig["emissPost_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_sigmaPriorFile_name(self, time, sectorName):
        return self.myConfig["sigmaPriorDir"] + "/" + self.myConfig["sigmaPrior_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_sigmaProcFile_name(self, time, sectorName):
        return self.myConfig["sigmaProcDir"] + "/" + self.myConfig["sigmaProc_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_HFile_name(self, time, recepTime, typeName):
        return self.myConfig["HDir"] + "/" + self.myConfig["H_Prefix"] + "_" + typeName + "_t" + time.strftime("%Y-%m-%d_%H:%M:%S") + "_r" + recepTime.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_obsFile_name(self, recepTime):
        return self.myConfig["obsDir"] + "/" + self.myConfig["obs_Prefix"] + "_" + recepTime.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_DFile_name(self, typeName, **kwargs):
        return self.myConfig["DDir"] + "/" + self.myConfig["D_Prefix"] + "_" + typeName + ".nc"

    def get_EFile_name(self, typeName, **kwargs):
        return self.myConfig["EDir"] + "/" + self.myConfig["E_Prefix"] + "_" + typeName + ".nc"

    def get_RFile_offline_name(self, recepTime, **kwargs):
        return self.myConfig["RDir"] + "/" + self.myConfig["R_Prefix"] + "_" + recepTime.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_DbckFile_name(self):
        return self.myConfig["DbckDir"] + "/" + self.myConfig["Dbck_Prefix"] + ".nc"

    def get_EbckFile_name(self):
        return self.myConfig["EbckDir"] + "/" + self.myConfig["Ebck_Prefix"] + ".nc"

    def get_bckPriorFile_name(self, recepTime):
        return self.myConfig["bckPriorDir"] + "/" + self.myConfig["bckPrior_Prefix"] + "_" + recepTime.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_bckProcFile_name(self, recepTime):
        return self.myConfig["bckProcDir"] + "/" + self.myConfig["bckProc_Prefix"] + "_" + recepTime.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_sbckPriorFile_name(self, recepTime):
        return self.myConfig["sbckPriorDir"] + "/" + self.myConfig["sbckPrior_Prefix"] + "_" + recepTime.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_sbckProcFile_name(self, recepTime):
        return self.myConfig["sbckProcDir"] + "/" + self.myConfig["sbckProc_Prefix"] + "_" + recepTime.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"


    def get_initHxFile_name(self, recepTime):
        return self.myConfig["initConcDir"] + "/" + self.myConfig["initHx_Prefix"] + "_" + recepTime.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_initHxBckFile_name(self, recepTime):
        return self.myConfig["initConcDir"] + "/" + self.myConfig["initHxBck_Prefix"] + "_" + recepTime.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_finalHxFile_name(self, recepTime):
        return self.myConfig["finalConcDir"] + "/" + self.myConfig["finalHx_Prefix"] + "_" + recepTime.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_finalHxBckFile_name(self, recepTime):
        return self.myConfig["finalConcDir"] + "/" + self.myConfig["finalHxBck_Prefix"] + "_" + recepTime.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_postHourlyFile_name(self, time, sectorName):
        return self.myConfig["outHourlyDir"] + "/" + self.myConfig["postHourly_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_spostHourlyFile_name(self, time, sectorName):
        return self.myConfig["outHourlyDir"] + "/" + self.myConfig["spostHourly_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_postDailyFile_name(self, time, sectorName):
        return self.myConfig["outDailyDir"] + "/" + self.myConfig["postDaily_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_spostDailyFile_name(self, time, sectorName):
        return self.myConfig["outDailyDir"] + "/" + self.myConfig["spostDaily_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_postWeeklyFile_name(self, time, sectorName):
        return self.myConfig["outWeeklyDir"] + "/" + self.myConfig["postWeekly_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_spostWeeklyFile_name(self, time, sectorName):
        return self.myConfig["outWeeklyDir"] + "/" + self.myConfig["spostWeekly_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"
    
    def get_postMonthlyFile_name(self, time, sectorName):
        return self.myConfig["outMonthlyDir"] + "/" + self.myConfig["postMonthly_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_spostMonthlyFile_name(self, time, sectorName):
        return self.myConfig["outMonthlyDir"] + "/" + self.myConfig["spostMonthly_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_postAllFile_name(self, sectorName):
        return self.myConfig["outAllDir"] + "/" + self.myConfig["postAll_Prefix"] + "_" + sectorName + ".nc"

    def get_spostAllFile_name(self, sectorName):
        return self.myConfig["outAllDir"] + "/" + self.myConfig["spostAll_Prefix"] + "_" + sectorName + ".nc"

    def get_priorHourlyFile_name(self, time, sectorName):
        return self.myConfig["outHourlyDir"] + "/" + self.myConfig["priorHourly_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_spriorHourlyFile_name(self, time, sectorName):
        return self.myConfig["outHourlyDir"] + "/" + self.myConfig["spriorHourly_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_priorDailyFile_name(self, time, sectorName):
        return self.myConfig["outDailyDir"] + "/" + self.myConfig["priorDaily_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_spriorDailyFile_name(self, time, sectorName):
        return self.myConfig["outDailyDir"] + "/" + self.myConfig["spriorDaily_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_priorWeeklyFile_name(self, time, sectorName):
        return self.myConfig["outWeeklyDir"] + "/" + self.myConfig["priorWeekly_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_spriorWeeklyFile_name(self, time, sectorName):
        return self.myConfig["outWeeklyDir"] + "/" + self.myConfig["spriorWeekly_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"
    
    def get_priorMonthlyFile_name(self, time, sectorName):
        return self.myConfig["outMonthlyDir"] + "/" + self.myConfig["priorMonthly_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_spriorMonthlyFile_name(self, time, sectorName):
        return self.myConfig["outMonthlyDir"] + "/" + self.myConfig["spriorMonthly_Prefix"] + "_" + sectorName + "_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_priorAllFile_name(self, sectorName):
        return self.myConfig["outAllDir"] + "/" + self.myConfig["priorAll_Prefix"] + "_" + sectorName + ".nc"

    def get_spriorAllFile_name(self, sectorName):
        return self.myConfig["outAllDir"] + "/" + self.myConfig["spriorAll_Prefix"] + "_" + sectorName + ".nc"

    def make_prior_emiss(self, sectors_PriorSetting, start = None, end = None, dt = None, regridFun = None):

        """
            sectors_PriorSetting:
              sectors_PriorSetting is a keyword dictionaries. For each sector whose dictionary contains:
                * emiss: Emiss name like "MEIC", "NDRC", "FFDAS" ...
                * getFun_kwargs: keyword arguments dictionary for the interface of emiss data API.

                If the "hasInterpolated" in emiss data API is False, the following is necessary:
                    * LON: New grid longitude (2D)
                    * LAT: New grid latitude (2D)

            start: Start time of emission data processing.

            end: End time of emission data processing.

            dt: Time interval of emission data processing.

            regridFun:
               User defined regrid function handle, being necessary if "interpolate_method" in emission data API is "user_defined".

        """

        if start is None:
            start = self.Start
        if end is None:
            end = self.End
        if dt is None:
            dt = self.dt
       
        parallelArgs = []
        for sectorName in sectors_PriorSetting:
            PriorSetting = sectors_PriorSetting[sectorName]
            emissName = PriorSetting["emiss"]

            if "LON" in PriorSetting:
                newLON = PriorSetting["LON"]
            else:
                newLON = None

            if "LAT" in PriorSetting:
                newLAT = PriorSetting["LAT"]
            else:
                newLAT = None

            if "getFun_kwargs" in PriorSetting:
                getFun_kwargs = PriorSetting["getFun_kwargs"]
            else:
                getFun_kwargs = {}

            getFun_kwargs["sector"] = sectorName

            loc = locals()
            exec("from ...emiss." + emissName + " import interface, hasInterpolated, interpolate_method")
            getDataFun = loc["interface"]
            hasInterpolated = loc["hasInterpolated"]
            interpolate_method = loc["interpolate_method"]
            isRegrid = not(hasInterpolated)

            current = start
            while(current <= end):

                #print(current)
                
                #getDataFun = lambda time: interface(time, **getFun_kwargs)
                outFile = self.get_priorFile_name(current, sectorName)
                # Arguments of function process_emiss:
                # Current, getDataFun, OutFile, writeFun, isRegrid, regridMethod, regridFun, newLON, newLAT 
                processEmissKwargs = {
                    "Current": current, "getDataFun": getDataFun, "getDataFun_kwargs": getFun_kwargs, "OutFile": outFile, "writeFun": self.X_Sigma_Indicator_write,
                    "isRegrid": isRegrid, "regridMethod": interpolate_method, "regridFun": regridFun, "newLON": newLON, "newLAT": newLAT,
                }
                parallelArgs.append(processEmissKwargs)
                current += dt

        pool = mtp.Pool(self.nProc)
        pool.map(self.process_emiss, parallelArgs)
        pool.close()
        pool.join()

    

    def make_prior_sigma(self, sectors_SigmaSetting, start = None, end = None, dt = None):

        """
                sectors_SigmaSetting:
                  sectors_SigmaSetting is a keyword dictionaries. For each sector whose dictionary contains:

                    * getSigmaFun: Function handle using for computing error standard deviation (sigma) of prior emission.
                    * getFun_kwargs: keyword arguments dictionary for function "getSigmaFun".

                    * isRegrid: True if sigma value computed by getSigmaFun should be regrided.

                    If isRegrid is True, the following is necessary:
                        * LON: New grid longitude (2D)
                        * LAT: New grid latitude (2D)
                        * regridMethod: regrid method ("linear", "nearest, "cubic", "user_defined")

                        If regridMethod is "user_defined", the following is necessary:
                            * regridFun: Handle of the user-defined regridding function.

                start: Start time of sigma values computing.

                end: End time of sigma values computing.

                dt: Time interval of sigma values computing.

        """

        if start is None:
            start = self.Start
        if end is None:
            end = self.End
        if dt is None:
            dt = self.dt

        parallelArgs = []

        for sectorName in sectors_SigmaSetting:

            SigmaSetting = sectors_SigmaSetting[sectorName]
            getSigmaFun = SigmaSetting["getSigmaFun"]

            if "getFun_kwargs" in SigmaSetting:
                getFun_kwargs = SigmaSetting["getFun_kwargs"]
            else:
                getFun_kwargs = {}

            if "isRegrid" in SigmaSetting:
                isRegrid = SigmaSetting["isRegrid"]
            else:
                isRegrid = False

            if "LON" in SigmaSetting:
                newLON = SigmaSetting["LON"]
            else:
                newLON = None

            if "LAT" in SigmaSetting:
                newLAT = SigmaSetting["LAT"]
            else:
                newLAT = None

            if "regridMethod" in SigmaSetting:
                regridMethod = SigmaSetting["regridMethod"]
            else:
                regridMethod = "linear"

            if "regridFun" in SigmaSetting:
                regridFun = SigmaSetting["regridFun"]
            else:
                regridFun = None

            current = start
            while(current <= end):
                #print(current)
                #print(getFun_kwargs)
                #set_trace()
                #getDataFun = lambda time: getSigmaFun(time = time, **getFun_kwargs)
                outFile = self.myConfig["sigmaPriorDir"] + "/" + self.myConfig["sigmaPrior_Prefix"] + "_" + sectorName + "_" + current.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"
                # Arguments of function process_emiss:
                # Current, getDataFun, OutFile, writeFun, isRegrid, regridMethod, regridFun, newLON, newLAT 
                processEmissKwargs = {
                    "Current": current, "getDataFun": getSigmaFun, "getDataFun_kwargs": getFun_kwargs, "OutFile": outFile, "writeFun": self.X_Sigma_Indicator_write,
                    "isRegrid": isRegrid, "regridMethod": regridMethod, "regridFun": regridFun, "newLON": newLON, "newLAT": newLAT,
                }
                #print(processEmissKwargs["getDataFun_kwargs"])
                #self.process_emiss(processEmissKwargs)
                parallelArgs.append(processEmissKwargs)
                current += dt

        pool = mtp.Pool(self.nProc)
        pool.map(self.process_emiss, parallelArgs)
        pool.close()
        pool.join()

    def copy_files(self, get_File_from, get_File_to, fun_from_kwargs, fun_to_kwargs):
        from_dir = get_File_from(**fun_from_kwargs)
        to_dir = get_File_to(**fun_to_kwargs)
        os.system("cp -p " + from_dir + " " + to_dir)

    def move_files(self, get_File_from, get_File_to, fun_from_kwargs, fun_to_kwargs):
        from_dir = get_File_from(**fun_from_kwargs)
        to_dir = get_File_to(**fun_to_kwargs)
        os.system("mv " + from_dir + " " + to_dir)

    def link_files(self, get_File_from, get_File_to, fun_from_kwargs, fun_to_kwargs):
        from_dir = get_File_from(**fun_from_kwargs)
        to_dir = get_File_to(**fun_to_kwargs)
        os.system("ln -sf " + from_dir + " " + to_dir)

    def copy_emiss_sigma(self, emiss_or_sigma):
        assert emiss_or_sigma.lower() in ["emiss", "sigma"]

        if emiss_or_sigma.lower() == "emiss":
            get_File_from = self.get_priorFile_name
            get_File_to = self.get_procFile_name

        if emiss_or_sigma.lower() == "sigma":
            get_File_from = self.get_sigmaPriorFile_name
            get_File_to = self.get_sigmaProcFile_name

        current = self.Start
        while(current <= self.End):
            for sector in self.sectors:
                kwargs = {"time": current, "sectorName": sector}
                self.copy_files(get_File_from, get_File_to, fun_from_kwargs = kwargs, fun_to_kwargs = kwargs)
            current += self.dt

    def copy_emiss_final(self, copy_or_link):

        assert copy_or_link.lower() in ["copy", "link"]

        get_File_from = self.get_procFile_name
        get_File_to = self.get_postFile_name

        current = self.Start
        while(current <= self.End):
            for sector in self.sectors:
                kwargs = {"time": current, "sectorName": sector}
                if copy_or_link.lower() == "copy":
                    self.copy_files(get_File_from, get_File_to, fun_from_kwargs = kwargs, fun_to_kwargs = kwargs)
                if copy_or_link.lower() == "link":
                    self.link_files(get_File_from, get_File_to, fun_from_kwargs = kwargs, fun_to_kwargs = kwargs)

            current += self.dt
        

    def make_H_offline(self, types_HSetting, objFootDic):

        """
            Arguments:
                * types_HSetting: Dictionaries for each type of H.
                    for each type (i.e. types_HSetting[type_]):
                        "backtime_j2i": Function handle provided in optimization core.
                        "recepTimeList": timeList of receptors.

                * objFootDic: Footprint objects dictionary for each type of H.
        """

        parallelArgs = []
        for type_ in types_HSetting:

            HSetting = types_HSetting[type_]
            assert "backtime_j2i" in HSetting
            #assert "objFootDic" in HSetting

            backtime_j2i = HSetting["backtime_j2i"]
            recepTimeList = HSetting["recepTimeList"]
            objFoot = objFootDic[type_]

            if recepTimeList is None:
                recepTimeList = objFoot.recepTimeList

            for recepTime in recepTimeList:
                #print("Receptor time", recepTime)

                timeList = backtime_j2i(recepTime)
                for time in timeList:

                    process_H_kwargs = {
                        "time": time, "recepTime": recepTime, "objFoot": objFoot,
                        "Htype": type_, "getHFileName": self.get_HFile_name,
                    }
                    parallelArgs.append(process_H_kwargs)
                    #process_H(process_H_kwargs)

        pool = mtp.Pool(self.nProc)
        pool.map(self.process_H, parallelArgs)
        pool.close()
        pool.join()
            
            

    def get_prior_emiss(self, time, sectorName, outType = "orig"):
        assert outType in ["vector", "orig"]
        if outType == "orig":
            emiss, LON, LAT = get_ncvar(self.get_priorFile_name(time, sectorName), ["emiss", "LON", "LAT"])
            return emiss, LON, LAT
        if outType == "vector":
            data = get_ncvar(self.get_priorFile_name(time, sectorName), ["data"])
            return data
    
    def get_prior_sigma(self, time, sectorName, outType = "orig"):
        assert outType in ["vector", "orig"]
        if outType == "orig":
            emiss, LON, LAT = get_ncvar(self.get_sigmaPriorFile_name(time, sectorName), ["emiss", "LON", "LAT"])
            return emiss, LON, LAT
        if outType == "vector":
            data = get_ncvar(self.get_sigmaPriorFile_name(time, sectorName), ["data"])
            return data

    def get_posterior_emiss(self, time, sectorName, outType = "orig"):
        assert outType in ["vector", "orig"]
        if outType == "orig":
            emiss, LON, LAT = get_ncvar(self.get_postFile_name(time, sectorName), ["emiss", "LON", "LAT"])
            return emiss, LON, LAT
        if outType == "vector":
            data = get_ncvar(self.get_postFile_name(time, sectorName), ["data"])
            return data
    
    def get_posterior_sigma(self, time, sectorName, outType = "orig"):
        assert outType in ["vector", "orig"]
        if outType == "orig":
            emiss, LON, LAT = get_ncvar(self.get_sigmaProcFile_name(time, sectorName), ["emiss", "LON", "LAT"])
            return emiss, LON, LAT
        if outType == "vector":
            data = get_ncvar(self.get_sigmaProcFile_name(time, sectorName), ["data"])
            return data

    def make_obs(self, *args, **kwargs):
        pass

    def make_E(self, types_ESetting):
        """
            type_ESetting: dictionaries for each type of E matrix:
                ESetting: (type_ESetting["site"] for example)
                    * method: Method of E matrix computing algorithm.
                    * method_kwargs: Keyword dictionary of compute_E function.
        """
        for type_ in types_ESetting:
            ESetting = types_ESetting[type_]
            assert "method" in ESetting
            assert "method_kwargs" in ESetting
            method = ESetting["method"]
            method_kwargs = ESetting["method_kwargs"]
            EDir = self.get_EFile_name(typeName = type_)
            if os.path.exists(EDir):
                print("E matrix file: " + EDir + " already exists, skip...")
            else:
                compute_E(objExp = self, method = method, typeName = type_, **method_kwargs)

    def make_D(self, types_DSetting):
        """
            type_DSetting: dictionaries for each type of D matrix:
                DSetting: (type_DSetting["site"] for example)
                    * method: Method of D matrix computing algorithm.
                    * method_kwargs: Keyword dictionary of compute_D function.
        """
        for type_ in types_DSetting:
            DSetting = types_DSetting[type_]
            assert "method" in DSetting
            assert "method_kwargs" in DSetting
            method = DSetting["method"]
            method_kwargs = DSetting["method_kwargs"]
            DDir = self.get_DFile_name(typeName = type_)
            compute_D(objExp = self, method = method, typeName = type_, **method_kwargs)

    def make_Ebck(self, EbckSetting):
        assert "method" in EbckSetting
        assert "method_kwargs" in EbckSetting
        method = EbckSetting["method"]
        method_kwargs = EbckSetting["method_kwargs"]
        compute_Ebck(objExp = self, method = method, **method_kwargs)

    def make_Dbck(self, DbckSetting):
        assert "method" in DbckSetting
        assert "method_kwargs" in DbckSetting
        method = DbckSetting["method"]
        method_kwargs = DbckSetting["method_kwargs"]
        compute_Dbck(objExp = self, method = method, **method_kwargs)

    def make_R_offline(self, RSetting):
        assert "method" in RSetting
        assert "method_kwargs" in RSetting
        method = RSetting["method"]
        method_kwargs = RSetting["method_kwargs"]
        for recepTime in self.objCore.objIter.iterTimeList:
            compute_R(recepTime = recepTime, funName = method, get_RFile_name = self.get_RFile_offline_name, **method_kwargs)
    
    def X_Sigma_Indicator_write(self, *args, **kwargs):
        X_Sigma_Indicator_write(*args, **kwargs)
    
    def process_emiss(self, kwargs):
        process_emiss(self, kwargs)

    def process_H(self, kwargs):
        process_H(self, kwargs)

    def compute_sigma(self, time, **kwargs):
        return compute_sigma(time = time, **kwargs)

    def compute_init_concentration(self, *args, **kwargs):
        compute_init_concentration(self, *args, **kwargs)

    def compute_final_concentration(self, *args, **kwargs):
        compute_final_concentration(self, *args, **kwargs)

    def make_output(self, *args, **kwargs):
        self.timeScale_timeList = make_output(self, *args, **kwargs)

    def build_output_dic(self, *args, **kwargs):
        build_output_dic(self, *args, **kwargs)

    def exp_run(self):
        #"""
        #--- Initialization ---#
        self.exp_initialize()

        #--- Main loop of optimization ---#
        for recepTime in self.objCore.objIter:

            #-- Calling optimize method of optimization core --#
            print("Optimize: ", recepTime)
            self.objCore.optimize(recepTime = recepTime, objExp = self)
        #"""
        #--- Post processing ---#
        self.exp_postproc()


def expbis_constructor(self, *args, **kwargs):

    """
    Keyword:
      configure: A module loaded from Bayesian.main.configure. (necessary)

    Main step:
    1. Load configurations from given \"configure\" module.
    2. Create an object of optimization core.
    3. Create footprint objects for each type of H matrix.
    4. ...

    """

    # Loading some relating configuration.
    assert "configure" in kwargs, "Error! Keyword \"configure\" couldn't be found."
    config = kwargs["configure"]
    self.myConfig = config["ExpConfig"]["BIS"]
    self.coreName = self.myConfig["core"]
    self.sectors = self.myConfig["sectors"]
    self.typeToSector_H = self.myConfig["typeToSector_H"]
    self.typeToSector_E = self.myConfig["typeToSector_E"]
    self.typeToSector_D = self.myConfig["typeToSector_D"]
    self.HtypeToFootprint = self.myConfig["HtypeToFootprint"]
    self.nProc = self.myConfig["nProc"]
    self.Start = dtm.datetime.strptime(self.myConfig["start"], "%Y-%m-%d_%H:%M:%S")
    self.End = dtm.datetime.strptime(self.myConfig["end"], "%Y-%m-%d_%H:%M:%S")
    self.dtHrs = self.myConfig["dtHrs"]
    self.dt = dtm.timedelta(hours = self.dtHrs)
    #self.optBck = self.myConfig["optBck"]

    #if "ExpIndicators" not in kwargs or kwargs["ExpIndicators"] is None:
    #    self.indicators = {"Sigma": Sigma_Indicator, "X": X_Indicator}
   
    # Define the directories to make.
    self.directoriesToMake = [
            "emissDir", "emissPriorDir", "emissProcDir", "emissPostDir",
            "sigmaDir", "sigmaPriorDir", "sigmaProcDir", "HDir", "EDir",
            "DDir", "RDir", "obsDir", "initConcDir", "finalConcDir",
            "outHourlyDir", "outDailyDir", "outWeeklyDir", "outMonthlyDir",
            "outAllDir",
        ]


    # Build the object of optimization core.
    loc = locals()
    # from ...core.core_<coreName>.main_class import Core<coreName> as core
    exec("from ...core.core_" + self.coreName + ".main_class import Core" + self.coreName + " as core")
    core = loc["core"]
    self.objCore = core(start = self.Start, end = self.End, dtHrs = self.dtHrs, coreConfig = config["CoreConfig"], ExpIndicators = self.indicators, nProc = self.nProc)

    self.objFootDic = {}
    self.footConfig = {}
    for type_ in self.HtypeToFootprint:
        footName = self.HtypeToFootprint[type_]
        self.footConfig[type_] = config["FootConfig"][footName]
        loc = locals()
        exec("from ...footprint.footprint_" + footName + " import " + footName + " as footprint")
        footprint = loc["footprint"]
        self.objFootDic[type_] = footprint(**self.footConfig[type_])

    self.build_output_dic()
    

def process_emiss(objExp, kwargs):
    #Current, getDataFun, OutFile, writeFun, isRegrid, regridMethod, regridFun, newLON, newLAT 
    assert "Current" in kwargs
    assert "getDataFun" in kwargs
    assert "OutFile" in kwargs
    assert "isRegrid" in kwargs
    assert "writeFun" in kwargs

    Current = kwargs["Current"]
    getDataFun = kwargs["getDataFun"]
    outFile = kwargs["OutFile"]
    isRegrid = kwargs["isRegrid"]
    writeFun = kwargs["writeFun"]

    if "getDataFun_kwargs" in kwargs:
        getDataFun_kwargs = kwargs["getDataFun_kwargs"]
    else:
        getDataFun_kwargs = {}

    origEmiss, origLON, origLAT = getDataFun(Current, **getDataFun_kwargs)
    assert origEmiss.shape == origLON.shape and origLON.shape == origLAT.shape

    # Interpolate the original emission
    #origLON, origLAT = np.meshgrid(origLon, origLat)
    if isRegrid:
        assert "newLON" in kwargs
        assert "newLAT" in kwargs
        assert "regridMethod" in kwargs
        newLON = kwargs["newLON"]
        newLAT = kwargs["newLAT"]
        #newLON, newLAT = np.meshgrid(newLON, newLAT)

        regridMethod = kwargs["regridMethod"]
        if regridMethod.lower() == "user_defined":
            assert "regridFun" in kwargs
            regridFun = kwargs["regridFun"]
            emissOutput = regridFun(origLON = origLON, origLAT = origLAT, origEmiss = origEmiss, newLON = newLON, newLAT = newLAT)
        else:
            emissOutput = interpolate.griddata((origLON.flatten(), origLAT.flatten()), origEmiss.flatten(), (newLON, newLAT), method = regridMethod)
    else:
        #assert np.sum(np.abs(origLon - lon)) <= 1E-8
        #assert np.sum(np.abs(origLat - lat)) <= 1E-8
        emissOutput = origEmiss
        newLON = origLON
        newLAT = origLAT
        #newLON, newLAT = np.meshgrid(newLON, newLAT)
    #print(emissOutput.shape)

    emissOutput = np.where(np.isnan(emissOutput), 0, emissOutput)
    

    writeFun(outFile, emissOutput, newLON, newLAT)

    #if config.isShowProg:
    #    show_progress(iprog.value, Nprog.value)
    #    iprog.value += 1

def process_H(objExp, kwargs):
    # time, recepTime, Htype
    assert "time" in kwargs
    assert "recepTime" in kwargs
    assert "Htype" in kwargs
    assert "objFoot" in kwargs
    assert "getHFileName" in kwargs

    time = kwargs["time"]
    recepTime = kwargs["recepTime"]
    Htype = kwargs["Htype"]
    objFoot = kwargs["objFoot"]
    get_HFile_name = kwargs["getHFileName"]
    #print(time, recepTime, Htype)

    arrayH = []
    for receptor in objFoot.receptors:
        foot, _, _ = objFoot.get_footprint(time = time, recepTime = recepTime, receptorName = receptor)
        arrayH.append(foot.flatten())
    arrayH = np.array(arrayH)
    array_to_sparse_to_nc(arrayH, get_HFile_name(time = time, recepTime = recepTime, typeName = Htype))


def trans_emiss_unit(arr, LON, LAT):
    #--- umol/m^2/s => tCO2/cell/hour --#
    return arr * areaS(LON, LAT) * 3600 / 1.0e6 * 44 / 1.0e6


def trans_sigma_unit(arr, LON, LAT):
    return (arr * areaS(LON, LAT) * 3600 / 1.0e6 * 44 / 1.0e6) ** 2

def straight(arr, *args, **kwargs):
    return arr

def build_output_dic(objExp):
    outputType = objExp.myConfig["outputType"].lower()
    assert outputType in ["sum", "mean"] 
    isSum = outputType == "sum"
    isMean = outputType == "mean"
    objExp.make_output_dic = {
        "priorEmiss": 
        {
            "inputFun": objExp.get_prior_emiss,
            "outNameFun_hourly": objExp.get_priorHourlyFile_name,
            "outNameFun_daily": objExp.get_priorDailyFile_name,
            "outNameFun_weekly": objExp.get_priorWeeklyFile_name,
            "outNameFun_monthly": objExp.get_priorMonthlyFile_name,
            "outNameFun_all": objExp.get_priorAllFile_name,
            #--- umol/m^2/s => tCO2/cell/hour --#
            "transFun": trans_emiss_unit if isSum else straight,
            "outFun": straight,
        },
        "postEmiss":
        {
            "inputFun": objExp.get_posterior_emiss,
            "outNameFun_hourly": objExp.get_postHourlyFile_name,
            "outNameFun_daily": objExp.get_postDailyFile_name,
            "outNameFun_weekly": objExp.get_postWeeklyFile_name,
            "outNameFun_monthly": objExp.get_postMonthlyFile_name,
            "outNameFun_all": objExp.get_postAllFile_name,
            #--- umol/m^2/s => tCO2/cell/hour --#
            "transFun": trans_emiss_unit if isSum else straight,
            "outFun": straight,
        },
        "priorSigma":
        {
            "inputFun": objExp.get_prior_sigma,
            "outNameFun_hourly": objExp.get_spriorHourlyFile_name,
            "outNameFun_daily": objExp.get_spriorDailyFile_name,
            "outNameFun_weekly": objExp.get_spriorWeeklyFile_name,
            "outNameFun_monthly": objExp.get_spriorMonthlyFile_name,
            "outNameFun_all": objExp.get_spriorAllFile_name,
            #--- convert sigma to square ---#
            #--- umol/m^2/s => (tCO2/cell/hour)^2 --#
            "transFun": trans_sigma_unit if isSum else straight,
            #--- convert square to sigma ---#
            "outFun": np.sqrt,
        },
        "postSigma":
        {
            "inputFun": objExp.get_posterior_sigma,
            "outNameFun_hourly": objExp.get_spostHourlyFile_name,
            "outNameFun_daily": objExp.get_spostDailyFile_name,
            "outNameFun_weekly": objExp.get_spostWeeklyFile_name,
            "outNameFun_monthly": objExp.get_spostMonthlyFile_name,
            "outNameFun_all": objExp.get_spostAllFile_name,
            #--- convert sigma to square ---#
            #--- umol/m^2/s => (tCO2/cell/hour)^2 --#
            "transFun": trans_sigma_unit if isSum else straight,
            #--- convert square to sigma ---#
            "outFun": np.sqrt,
        },
    }

    
        
