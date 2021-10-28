#!/usr/bin/env python

#-> Application module (app) <-#
# The Beginning of the whole system

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from .main import configure as conf
from .main.main_class import BisCase

# Create an object of BIS case.
case = BisCase(conf)

# Showing the completation.
def show_complete():

    print("")
    print("===================================")
    print("!     Successful Complete BIS     !")
    print("===================================")
    print("If you have any question, please contact Wenhan Tang at:")
    print("-> tangwh@mail.iap.ac.cn")
    print("or -> tangwenhanisusu@gmail.com")
    print("")

# Drawing BenchPlots.
def run_benchplots():

    global case

    # Weather it should be run.
    isDrawBench = conf.isDrawBench
    if not isDrawBench:
        return

    # Import Bayeisan Inversion Analysis and Display System (BIADS).
    from .Analysis.get_benchClass import get_benchClass
    import os

    # Load experiment name.
    expName = conf.EXPNAME
    # Load experiment type.
    expType = conf.EXPTYPE
    # Is this experiment a free-type-experiment 
    isFree = conf.ISFREE
    # Load case name.
    caseName = conf.CASENAME
    # Load BenchPlots modules name list.
    drawList = conf.DRAWLIST
    # Load spatial plots time scale list.
    spatial_timeScaleList = conf.SPATIAL_TIMESCALELIST
    # Load points/regions plots time scale list.
    prTs_timeScaleList = conf.PRTS_TIMESCALELIST

    # Define the main directory of benchplots.
    benchDir = conf.BENCHDIR + "/" + caseName

    # Create benchplots directory.
    if not os.path.exists(benchDir):
        os.makedirs(benchDir)
    # Copy (and cover) the configure file for the experiment.
    os.system("cp -p " + conf.CASEDIR + "/configure.py " + benchDir + "/configure.txt")
    os.system("cp -rp " + conf.CASEDIR + "/stations " + benchDir + "/stations")

    # Load BenchPlots class based on the type pf experiment.
    BenchClass = get_benchClass(expName = expName, expType = expType, isFree = isFree)

    # Create an object of BenchPlots class.
    objBench = BenchClass(benchDir = benchDir, objCase = case)

    # Begin to draw benchplots
    objBench.run(drawList = drawList, spatial_timeScaleList = spatial_timeScaleList, prTs_timeScaleList = prTs_timeScaleList)

    # Creating tar file for benchplots (*.pdf, *.mp4).
    print("Make tar file ...")
    os.chdir(benchDir) 
    cmd = "tar -zcvf BenchPlots_" + caseName + ".tar.gz *.pdf *.mp4 configure.txt 2> /dev/null"
    os.system(cmd)

# The running function, which equals to the build-in function "app_run" of BisCase object.
def app_run(*args, **kwargs):

    global case

    # Begin at here.
    case.app_run(*args, **kwargs)

    # Drawing benchplots.
    run_benchplots()

    # Showing the completation.
    show_complete()

