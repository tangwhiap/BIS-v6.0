#!/usr/bin/env python

#-> Base class for REAL experiment object. <-#

# Authors:
#   Wenhan TANG - 07/2021
#   ...


from ...exp_base.main_class import ExpBIS
from ....utils.netcdf_io import obs_write, get_ncvar, nc_write
from ....utils.module2dic import module2dic

#from . import benchplots_configure as benchConfig

import numpy as np
import multiprocessing as mtp
from pdb import set_trace

class ExpREAL(ExpBIS):

    def __init__(self, *args, **kwargs):

        self.isOSSE = False
        self.isREAL = True

        self.hasBCK = True

        ExpBIS.__init__(self, *args, **kwargs)
        self.directoriesToMake += ["bckPriorDir", "bckProcDir", "sbckPriorDir", "sbckProcDir"]

        #-- Load configuration of REAL experiment --#
        self.myConfig = dict(**self.myConfig, **kwargs["configure"]["ExpConfig"]["REAL"])

        #-- Load observation configures --#
        self.obsConfig = kwargs["configure"]["obsConfig"]

        #-- Load backgroudn configures --#
        self.bckConfig = kwargs["configure"]["bckConfig"]

        ##-- Create observation and background object --#
        #self.build_observation_object()
        #self.build_background_object()

        #-- Load benchplots configures --#
        #self.benchConfig = module2dic(benchConfig)

    def build_observation_object(self):

        #-- Create object of observation --#
        obsName = self.myConfig["OBS"]
        loc = locals()
        exec("from ....observation." + obsName + ".main_class import OBS_" + obsName + " as obsClass")
        self.objOBS = loc["obsClass"](**self.obsConfig[obsName])

    def build_background_object(self):

        #-- Create object of background --#
        bckName = self.myConfig["BCK"]
        loc = locals()
        exec("from ....background." + bckName + ".main_class import BCK_" + bckName + " as bckClass")
        self.objBCK = loc["bckClass"](**self.bckConfig[bckName])

    def make_obs(self):
        obsTimeList = self.objCore.objIter.obsTimeList
        for time in obsTimeList:
            #print("make obs:", time)
            obs = self.objOBS.get_obs_proc(time)#, sites = self.objFoot.receptors)
            obs = obs.reshape(len(obs), 1)
            obs_write(outFile = self.get_obsFile_name(recepTime = time), obs = obs)
    

    #-- Create a background file for a single time --#
    def create_background_file(self, time):
        bck = self.objBCK.compute_BCK(time)
        if len(bck.shape) == 1:
            bck = bck.reshape(len(bck), 1)
        assert len(bck.shape) == 2
        nc_write(self.get_bckPriorFile_name(time), bck)


    #-- Create a background sigma value file for a single time --#
    def create_background_sigma_file(self, time):
        sbck = self.objBCK.compute_sigma_BCK(time)
        if len(sbck.shape) == 1:
            sbck = sbck.reshape(len(sbck), 1)
        assert len(sbck.shape) == 2
        nc_write(self.get_sbckPriorFile_name(time), sbck)

    def make_background(self):

        # Parallel running for makeing background files. --#
        parallelArgs = [time for time in self.objCore.objIter.obsTimeList]
        pool = mtp.Pool(self.nProc)
        pool.map(self.create_background_file, parallelArgs)
        pool.close()
        pool.join()

        # Parallel running for makeing background sigma files. --#
        pool = mtp.Pool(self.nProc)
        pool.map(self.create_background_sigma_file, parallelArgs)
        pool.close()
        pool.join()

    def copy_background_sigma(self, bck_or_sbck):
        assert bck_or_sbck.lower() in ["bck", "sbck"]
        if bck_or_sbck.lower() == "bck":
            get_File_from = self.get_bckPriorFile_name
            get_File_to = self.get_bckProcFile_name

        if bck_or_sbck.lower() == "sbck":
            get_File_from = self.get_sbckPriorFile_name
            get_File_to = self.get_sbckProcFile_name

        for time in self.objCore.objIter.obsTimeList:
            kwargs = {"recepTime": time}
            self.copy_files(get_File_from, get_File_to, fun_from_kwargs = kwargs, fun_to_kwargs = kwargs)


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
            11. Make E and D for background. (if hasBCK is True)
            12. Make R file.
        """
        
        
        #--- Create directories of this case. ---#
        self.case_mkdir()

        #--- Make prior emission files. ---#
        self.make_prior_emiss(self.sectors_PriorSetting)

        #--- Make prior emission uncertainty files. ---#
        self.make_prior_sigma(self.sectors_SigmaSetting)

        #--- Copy the prior emission files to their "proc" files. ---#
        self.copy_emiss_sigma("emiss")

        #-- Same for prior emission uncertainty files. ---#
        self.copy_emiss_sigma("sigma")

        #--- Make footprint (H) files. ---#
        self.make_H_offline(self.types_HSetting, self.objFootDic)

        #--- Make obervation data files. ---#
        self.make_obs()

        #--- Make Background data files. ---#
        self.make_background()

        #--- Copy Background data files to their "proc" files. ---#
        self.copy_background_sigma("bck")

        #--- Copy Background uncertainty files to their "proc" files. ---#
        self.copy_background_sigma("sbck")

        #--- Make spatial correlation (E matrix) file. ---#
        self.make_E(self.types_ESetting)

        #--- Make temporal correlation (D matrix) file. ---#
        self.make_D(self.types_DSetting)

        #--- If the case has background part. ---#
        if self.hasBCK:
            #-- Make spatial correlation file for background. --#
            self.make_Ebck(self.EbckSetting)
            #-- Make temporal correlation file for background. --#
            self.make_Dbck(self.DbckSetting)

        #--- Make model-data mismatch (R matrix) file. ---#
        self.make_R_offline(self.RSetting)

        #--- Compute site concentration alias of prior CO2 flux. ---#
        self.compute_init_concentration()

#def process_background(time, objBCK):

        

    

