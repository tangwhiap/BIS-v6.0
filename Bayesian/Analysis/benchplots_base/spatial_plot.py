#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...

from .draw_object_class import AnalysisField, AnalysisStations

from ..utils.Bench import BenchPlots_png as BenchPlots

import numpy as np
import matplotlib.pyplot as plt
import datetime as dtm
import multiprocessing as mtp
import os

from pdb import set_trace

def progress_bar(time, start, end, fig, axRec = [0.1, 0.05, 0.85, 0.2]):
    ax = fig.add_axes(axRec)#, axes_class = axisartist.Axes)
    #ax.set_zorder(0)
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_position(('data', 0))
    ax.set_yticks([])
    #ax.axis["y=0"].label.set_text("y = 0")
    ax.plot(time, 0, marker = "o", color = "red")
    ax.set_xlim(start, end)
    fig.autofmt_xdate()
    return ax

def stations_local_plot(objBench, staLon, staLat, gs, setting, fig = None, suptitle = ""):
    if fig is None:
        fig = plt.figure()

    if "stationsPlot_kwargs" not in setting:
        setting["stationsPlot_kwargs"] = {}

    if "scatter_kwargs" not in setting:
        setting["scatter_kwargs"] = {}


    shpDir = objBench.myConfig["shpDir"]
    #gses = fig.add_gridspec(ncols = ncols, nrows = nrows, **gs_kwargs)

    stationsPlot_kwargs = setting["stationsPlot_kwargs"]
    scatter_kwargs = setting["scatter_kwargs"]
    lon_s = objBench.objAna.sampleLON.min()
    lon_e = objBench.objAna.sampleLON.max()
    lat_s = objBench.objAna.sampleLAT.min()
    lat_e = objBench.objAna.sampleLAT.max()

    # Atloc setting
    #lon_s = 115.2
    #lon_e = 117.6
    #lat_s = 39.4
    #lat_e = 41.1
    ################

    spatialRange = {"lon_s": lon_s, "lon_e": lon_e, "lat_s": lat_s, "lat_e": lat_e}

    objStaPlot = AnalysisStations(lonList = staLon, latList = staLat, title = setting["title"])
    objStaPlot.stations_plot(fig = fig, gs = gs, shpDir = shpDir, spatialRange = spatialRange, stationsPlot_kwargs = stationsPlot_kwargs, scatter_kwargs = scatter_kwargs)

    return fig

def bis_spatial_plot(objBench, field_dic, nrows, ncols, setting, fig = None, suptitle = "", gs_kwargs = {}):
    if fig is None:
        fig = plt.figure()

    if "spatial_kwargs" not in setting:
        setting["spatial_kwargs"] = {}

    if "pcolormesh_kwargs" not in setting:
        setting["pcolormesh_kwargs"] = {}

    if "cbar_kwargs" not in setting:
        setting["cbar_kwargs"] = {}

    shpDir = objBench.myConfig["shpDir"]
    gses = fig.add_gridspec(ncols = ncols, nrows = nrows, **gs_kwargs)

    spatial_kwargs = setting["spatial_kwargs"]
    pcolormesh_kwargs = setting["pcolormesh_kwargs"]
    cbar_kwargs = setting["cbar_kwargs"]

    #setting = objBench.myConfig["movie_kwargs"]
    for gs, fieldName in zip(gses, field_dic):
        #print(fieldName)
        array, LON, LAT = field_dic[fieldName]

        fieldSetting = setting[fieldName + "_kwargs"]
        
        objField = AnalysisField(array = array, LON = LON, LAT = LAT, title = fieldSetting["title"], unit = fieldSetting["unit"], vmin = fieldSetting["vminmax"][0], vmax = fieldSetting["vminmax"][1])

        fig, geo_ax, cs, cbar = objField.spatial_plot(fig = fig, gs = gs, shpDir = shpDir, spatial_kwargs = spatial_kwargs, pcolormesh_kwargs = pcolormesh_kwargs, cbar_kwargs = cbar_kwargs)

    fig.suptitle(suptitle)

    return fig

def draw_fieldPlots_snap(objBench, fig, time, timeScale = "daily"):

    fieldPlots_kwargs = objBench.myConfig["fieldPlots_kwargs"]
    vminmax_dic = objBench.myConfig["fieldPlots_vminmax"][timeScale if objBench.isSum else "mean"]
    shpDir = objBench.myConfig["shpDir"]
    field_dic = objBench.build_field_dic(time, timeScale)
    timeUnit_dic = objBench.myConfig["timeUnit_dic"]

    for fieldName in field_dic:
        if objBench.isSum:
            if fieldName == "reducedUncertainty":
                fieldPlots_kwargs[fieldName + "_kwargs"]["unit"] = "(KtCO2/km^2/" + timeUnit_dic[timeScale] + ")^2"
            else:
                fieldPlots_kwargs[fieldName + "_kwargs"]["unit"] = "KtCO2/km^2/" + timeUnit_dic[timeScale]
        if objBench.isMean:
            fieldPlots_kwargs[fieldName + "_kwargs"]["unit"] = "umol/m^2/s"
        fieldPlots_kwargs[fieldName + "_kwargs"]["vminmax"] = vminmax_dic[fieldName]

    suptitle = time.strftime("%Y-%m-%d %H:%M:%S (UTC)")
    nrows = (len(field_dic) + 1) // 2
    objBench.bis_spatial_plot(field_dic = field_dic, nrows = nrows, ncols = 2, fig = fig, suptitle = suptitle, setting = fieldPlots_kwargs)

def draw_fieldPlots(objBench, drawStart, drawEnd, timeScale = "daily"):
    assert timeScale.lower() in ["daily", "monthly", "weekly", "all"]
    myConfig = objBench.myConfig["fieldPlots_kwargs"]
    bcp = BenchPlots(figsize = myConfig["figsize"], save_name = objBench.benchDir + "/" + myConfig["PlotsName_Prefix"] + "_" + timeScale)
    timeList = objBench.objAna.timeList_dic[timeScale]
    timeList = np.array(timeList)[(np.array(timeList) >= drawStart) & (np.array(timeList) <= drawEnd)]
    for time in timeList:
        print("drawing " + time.strftime("%Y-%m-%d %H:%M:%S (UTC)"))
        objBench.draw_fieldPlots_snap(fig = bcp.cfig, time = time, timeScale = timeScale)
        bcp.next_page()
    bcp.close()
        



def draw_movie_snap(objBench, args):

    time, movieStart, movieEnd, snapDir = args
    print(time)
    movie_kwargs = objBench.myConfig["movie_kwargs"]
    vminmax_dic = objBench.myConfig["movie_vminmax"]
    shpDir = objBench.myConfig["shpDir"]

    field_dic = objBench.build_movie_field_dic(time)

    fig = plt.figure(figsize = movie_kwargs["figsize"])

    progress_bar(time, movieStart, movieEnd, fig, axRec = [0.1, 0.01, 0.85, 0.15])

    suptitle = time.strftime("%Y-%m-%d %H:%M:%S (UTC)")

    nrows = (len(field_dic) + 1) // 2

    for fieldName in field_dic:
        movie_kwargs[fieldName + "_kwargs"]["vminmax"] = vminmax_dic[fieldName]

    objBench.bis_spatial_plot(field_dic = field_dic, nrows = nrows, ncols = 2, fig = fig, suptitle = suptitle, setting = movie_kwargs)
    
    fig.savefig(snapDir + "/SpatialMovie_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".png", dpi = movie_kwargs["dpi"])

    plt.close(fig)

def make_cartoon(objBench, movieStart, movieEnd, movieDtHrs):

    movie_kwargs = objBench.myConfig["movie_kwargs"]
    movieDir = objBench.benchDir
    movieName = movie_kwargs["movieName"]
    snapDir = objBench.benchDir + "/movie_snaps"
    dt = dtm.timedelta(hours = movieDtHrs)

    if not os.path.exists(snapDir):
        os.makedirs(snapDir)

    parallelArgs = []
    current = movieStart
    while(current <= movieEnd):
        parallelArgs.append((current, movieStart, movieEnd, snapDir))
        current += dt

    pool = mtp.Pool(objBench.nProc)
    pool.map(objBench.draw_movie_snap, parallelArgs)
    pool.close()
    pool.join()
    
    

    if os.path.exists(movieDir + "/" + movieName + ".mp4"):
        print("Deleting " + movieDir + "/" + movieName + ".mp4")
        os.system("rm -f " + movieDir + "/" + movieName + ".mp4")
    #cmd = "ffmpeg -framerate 20 -pattern_type glob -i \"" + snapDir + "/spatial_*.png\" -c:v libx264 -preset ultrafast -pix_fmt yuv420p -r 10 " + movieDir + "/" + movieName + ".mp4"
    cmd = "ffmpeg -framerate 20 -pattern_type glob -i \"" + snapDir + "/SpatialMovie_*.png\" -pix_fmt yuv420p -r 10 " + movieDir + "/" + movieName + ".mp4"
    print(cmd)
    os.system(cmd)
    print("rm -f " + snapDir + "/SpatialMovie_*.png")
    #os.system("rm -f " + snapDir + "/SpatialMovie_*.png")
    #os.system("rmdir " + snapDir)

def build_movie_field_dic(objBench, time):

    objAna = objBench.objAna
    field_dic = {}
    field_dic["priorEmiss"] = objAna.get_prior_emiss(time)
    field_dic["posteriorEmiss"] = objAna.get_posterior_emiss(time)
    field_dic["reducedUncertainty"] = objAna.get_reduced_uncertainty_time(time)
    field_dic["posteriorMprior"] = field_dic["posteriorEmiss"][0] - field_dic["priorEmiss"][0], field_dic["priorEmiss"][1], field_dic["priorEmiss"][2]
    return field_dic

def build_field_dic(objBench, time, timeScale):

    objAna = objBench.objAna
    field_dic = {}
    priorEmiss, LON, LAT = objAna.get_prior_emiss(time, timeScale = timeScale)
    posteriorEmiss, _, _ = objAna.get_posterior_emiss(time, timeScale = timeScale)
    priorSigma, _, _ = objAna.get_prior_sigma(time, timeScale = timeScale)
    posteriorSigma, _, _ = objAna.get_posterior_sigma(time, timeScale = timeScale)
    reducedUncertainty = priorSigma - posteriorSigma
    posteriorMprior = posteriorEmiss - priorEmiss
    
    areaS = objAna.sampleAreaS

    if objBench.isSum:
        priorEmiss *= (0.001 * areaS / 1.0e6)
        posteriorEmiss *= (0.001 * areaS / 1.0e6)
        reducedUncertainty *= ((0.001 * areaS / 1.0e6) ** 2)
        posteriorMprior *= (0.001 * areaS / 1.0e6)

    field_dic["priorEmiss"] = priorEmiss, LON, LAT
    field_dic["posteriorEmiss"] = posteriorEmiss, LON, LAT
    field_dic["reducedUncertainty"] =  reducedUncertainty, LON, LAT
    field_dic["posteriorMprior"] = posteriorMprior, LON, LAT

    return field_dic
