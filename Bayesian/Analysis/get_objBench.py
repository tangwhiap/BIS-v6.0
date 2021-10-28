#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...



from .get_benchClass import get_benchClass
from ..main.main_class import BisCase

import os

procDir = "/".join(__file__.split("/")[:-1])  + "/proc"

def get_objBench(benchDir, caseDir, caseName):

    configDir = caseDir + "/configure.py"
    configPackDir = caseDir + "/stations"
    assert os.path.exists(configDir)
    assert os.path.exists(configPackDir)

    configNewDir = procDir + "/" + caseName

    if not os.path.exists(configNewDir):
        os.makedirs(configNewDir)

    os.system("cp -p " + configDir + " " + configNewDir)
    os.system("cp -rp " + configPackDir + " " + configNewDir)

    loc = locals()
    exec("from .proc." + caseName + " import configure as config")
    caseConfig = loc["config"]

    objCase = BisCase(caseConfig, directoryReserve = True)

    assert caseName == caseConfig.CASENAME, "Conflict case name between: " + caseName + " and " + caseConfig.CASENAME
    expName = caseConfig.EXPNAME
    expType = caseConfig.EXPTYPE
    isFree = caseConfig.ISFREE
    
    benchClass = get_benchClass(expName = expName, expType = expType, isFree = isFree)

    objBench = benchClass(benchDir = benchDir, objCase = objCase)
    
    return objBench

    

