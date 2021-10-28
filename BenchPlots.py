#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...


#--- Import BenchPlots package. ---#
from Bayesian.Analysis.get_objBench import get_objBench

import os

#--- CASE name ---#
#caseName = "afternoon_mean_test"
caseName = "PIS_Ls3_4DVar_trad"

#--- CASE directory ---#
caseDir = "/home/tangwh/modeldata/BIS_cases/" + caseName


#--- A list to tell me which part of benchplots you want to have. ---#
drawList = [
    #"movie",  # Spatial movie.
    #"spatial_plots", # Spatial plots.
    "site_timeSeries", # Time series of site.
    #"pointRegion_timeSeries", # Time series of some points/regions. (The points/regions is defined in benchplots_configure.py)
]

#--- Different time frequency for spatial plots ---#
# ["daily", "weekly", "monthly", "all"]
spatial_timeScaleList = [
    #"daily",
    "weekly",
]

#--- Different time frequency for points/regions time series plots ---#
# ["hourly", "daily", "weekly", "monthly"]
prTs_timeScaleList = [
    "hourly",
    #"daily", 
    #"weekly",
] 

isMakeTar = True

#--- BenchPlots root directory ---#
benchplotsDir = "/home/tangwh/public_html/BIS_show"

#--- Directory of this experiment's benchplots ---#
benchDir = benchplotsDir + "/" + caseName

#--- Create benchplots directory ---#
if not os.path.exists(benchDir):
    os.makedirs(benchDir)

#--- Copy (and cover) the configure file for the experiment. ---#
print("cp -p " + caseDir + "/configure.py " + benchDir + "/configure.txt")
os.system("cp -p " + caseDir + "/configure.py " + benchDir + "/configure.txt")


#--- Create an object of BenchPlots class. ---#
objBench = get_objBench(benchDir = benchDir, caseDir = caseDir, caseName = caseName)

#--- Begin to draw benchplots. ---#
objBench.run(drawList = drawList, spatial_timeScaleList = spatial_timeScaleList, prTs_timeScaleList = prTs_timeScaleList)

#--- Creating tar file for benchplots (*.pdf, *.mp4). ---#
if isMakeTar:
    print("Make tar file ...")
    os.chdir(benchDir) 
    cmd = "tar -zcvf BenchPlots_" + caseName + ".tar.gz *.pdf *.mp4 configure.txt 2> /dev/null"
    os.system(cmd)

