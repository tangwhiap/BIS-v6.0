#!/usr/bin/env python

# Auhtors:
#   Wenhan TANG - 07/2021
#   ...


from .compute_Hx import Hx, Hx2obs, compute_Hx
from .compute_HQ import compute_HQ, compute_HQHt
from .compute_HdQ import compute_HdQ
from .compute_HQHt_R import compute_HQHt_R
#from .compute_INV import compute_INV
from .compute_d import compute_d
from .inversion import inversion
from .background import background_optimize

import os

class CoreOptimize(object):
    
    def __init__(self, start, end, coreConfig, ExpIndicators, dtHrs = 1, nProc = 1):

        self.start = start
        self.end = end
        self.dtHrs = 1
        self.nProc = nProc
        #self.recepTime = recepTime
        self.myConfig = coreConfig["optimize"]
        self.optBck = self.myConfig["optBck"]
        self.ExpIndicators = ExpIndicators
        self.directoriesToMake = [
            "workDir"#, "HQ_Dir",
        ]
        self.subWorkDirs = [
            "HQ", 
            "HdQ",
            "reduced",
            "INV_HQ",
            "Hi_HQti",
        ]
        self.make_dirs()

    def make_dirs(self):
        for dirName in self.directoriesToMake:
            dir_ = self.myConfig[dirName]
            if not os.path.exists(dir_):
                os.makedirs(dir_)

    def get_workdir(self, recepTime):
        return self.myConfig["workDir"] + "/" + recepTime.strftime("%Y%m%d%H")

    def make_workdir(self, recepTime):
        workdir = self.get_workdir(recepTime = recepTime)
        if not os.path.exists(workdir):
            os.makedirs(workdir)
        for dirName in self.subWorkDirs:
            if not os.path.exists(workdir + "/" + dirName):
                os.makedirs(workdir + "/" + dirName)

        return workdir

    def get_HQFile_name(self, jtime, recepTime, sectorName):
        return self.get_workdir(recepTime) + "/HQ/" + self.myConfig["HQ_Prefix"] + "_" + sectorName + "_" + jtime.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_HdQFile_name(self, jtime, recepTime, sectorName):
        return self.get_workdir(recepTime) + "/HdQ/" + self.myConfig["HdQ_Prefix"] + "_" + sectorName + "_" + jtime.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_HQHtFile_name(self, recepTime, sectorName):
        return self.get_workdir(recepTime) + "/" +  self.myConfig["HQHt_Prefix"] + "_" + sectorName + ".nc"

    def get_HQHtRFile_name(self, recepTime):
        return self.get_workdir(recepTime) + "/" +  self.myConfig["HQHtR_Prefix"] +  ".nc"

    #def get_QbckFile_name(self, recepTime):
    #    return self.get_workdir(recepTime) + "/" + self.myConfig["Qbck_Prefix"] + ".nc"

    def get_INVFile_name(self, recepTime):
        return self.get_workdir(recepTime) + "/" + self.myConfig["INV_Prefix"] + ".nc"

    def get_INVHQFile_name(self, jtime, recepTime, sectorName):
        return self.get_workdir(recepTime) + "/INV_HQ/" + self.myConfig["INV_HQ_Prefix"] + "_" + sectorName + "_" + jtime.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_HiHQtiFile_name(self, recepTime, hisTime, sectorName):
        return self.get_workdir(recepTime) + "/Hi_HQti/" + self.myConfig["HiHQti_Prefix"] + "_" + sectorName + "_" + hisTime.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def get_priorHxFile_name(self, recepTime):
        return self.get_workdir(recepTime) + "/" + self.myConfig["HxPrior_Prefix"] + ".nc"

    def get_procHxFile_name(self, recepTime):
        return self.get_workdir(recepTime) + "/" + self.myConfig["HxProc_Prefix"] + ".nc"

    def get_priorHxBckFile_name(self, recepTime):
        return self.get_workdir(recepTime) + "/" + self.myConfig["HxBckPrior_Prefix"] + ".nc"

    def get_procHxBckFile_name(self, recepTime):
        return self.get_workdir(recepTime) + "/" + self.myConfig["HxBckProc_Prefix"] + ".nc"

    def get_dFile_name(self, recepTime):
        return self.get_workdir(recepTime) + "/" + self.myConfig["d_Prefix"] + ".nc"

    def get_dsigma2File_name(self, time, recepTime, sectorName):
        return self.get_workdir(recepTime) + "/reduced/" + self.myConfig["dsigma2_Prefix"] + "_" + sectorName + "_"+ time.strftime("%Y-%m-%d_%H:%M:%S") + ".nc"

    def Hx(self, *args, **kwargs):
        return Hx(*args, **kwargs)

    def Hx2obs(self, *args, **kwargs):
        return Hx2obs(*args, **kwargs)

    def compute_Hx(self, *args, **kwargs):
        return compute_Hx(*args, **kwargs)

    def compute_HQ(self, *args, **kwargs):
        compute_HQ(self, *args, **kwargs)

    def compute_HdQ(self, *args, **kwargs):
        compute_HdQ(self, *args, **kwargs)

    def compute_HQHt(self, *args, **kwargs):
        compute_HQHt(self, *args, **kwargs)

    def compute_HQHt_R(self, *args, **kwargs):
        compute_HQHt_R(self, *args, **kwargs)

    def compute_INV(self, *args, **kwargs):
        compute_INV(self, *args, **kwargs)

    def compute_d(self, *args, **kwargs):
        compute_d(self, *args, **kwargs)

    def background_optimize(self, *args, **kwargs):
        background_optimize(self, *args, **kwargs)

    def inversion(self, *args, **kwargs):
        inversion(self, *args, **kwargs)
        
    #def optimize(self, recepTime, objExp):
    #    self.make_workdir(recepTime)

