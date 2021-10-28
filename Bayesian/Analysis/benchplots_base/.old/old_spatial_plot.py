#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...

import numpy as np
import netCDF4 as nc
import datetime as dtm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
import glob
import os
import multiprocessing as mtp
from pdb import set_trace


#ShpDir = "/home/tangwh/datasets/china_shp/cnhimap.shp"

def spatial_plot(objBench, Arr, LON, LAT, fig, gs, cmap = "jet", cmap_va = False, draw_gridline = True, draw_lonlat_label = True, **pcolormesh_kwargs):

    def alpha_vary_cmap(colormap):
        cmap = plt.get_cmap(colormap)
        my_cmap = cmap(np.arange(cmap.N))
        my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
        my_cmap = ListedColormap(my_cmap)
        return my_cmap

    assert Arr.shape == LON.shape and LON.shape == LAT.shape

    proj = ccrs.PlateCarree()
    geo_ax = fig.add_subplot(gs, projection = proj)
    lon_s = LON.min(); lon_e = LON.max() 
    lat_s = LAT.min(); lat_e = LAT.max()

    geo_ax.set_extent([lon_s, lon_e, lat_s, lat_e], crs = proj)

    if draw_gridline:
        gl = geo_ax.gridlines(crs = proj, linestyle = "--", alpha = 0.5, draw_labels = True)
        gl.top_labels = False
        gl.right_labels = False
        if not draw_lonlat_label:
            gl.left_lables = False
            gl.bottom_lables = False

    geo_ax.add_geometries(Reader(objBench.myConfig["ShpDir"]).geometries(), proj, facecolor = "none", edgecolor = "k", linewidth = 1)

    if cmap_va:
        cmap = alpha_vary_cmap(cmap)
    cs = geo_ax.pcolormesh(LON, LAT, Arr, cmap = cmap, **pcolormesh_kwargs)
    return geo_ax, cs

def bis_spatial_plot(objBench, LON, LAT, priorEmiss, posteriorEmiss, reducedUncertainty, posteriorMprior, truthEmiss = None, priorMtruth = None, plotTruth = False, fig = None, suptitle = "", vminmax = None, cmap = "coolwarm", gs_kwargs = {}, pcolormesh_kwargs = {}, cbar_kwargs = {}, **kwargs):

    if fig is None:
        fig = plt.figure()
    
    gs = fig.add_gridspec(ncols = 2, nrows = 3 if plotTruth else 2, **gs_kwargs)
    vmin = vminmax["priorEmiss"][0]
    vmax = vminmax["priorEmiss"][1]
    ax_prior, cs_prior = objBench.spatial_plot(priorEmiss, LON, LAT, gs = gs[0, 0], fig = fig, cmap = cmap, vmin = vmin, vmax = vmax, **pcolormesh_kwargs)

    vmin = vminmax["posteriorEmiss"][0]
    vmax = vminmax["posteriorEmiss"][1]
    ax_posterior, cs_posterior = objBench.spatial_plot(posteriorEmiss, LON, LAT, gs = gs[0, 1], fig = fig, cmap = cmap, vmin = vmin, vmax = vmax, **pcolormesh_kwargs)

    vmin = vminmax["reducedUncertainty"][0]
    vmax = vminmax["reducedUncertainty"][1]
    ax_reduced, cs_reduced = objBench.spatial_plot(reducedUncertainty, LON, LAT, gs = gs[1, 0], fig = fig, cmap = cmap, vmin = vmin, vmax = vmax, **pcolormesh_kwargs)

    vmin = vminmax["posteriorMprior"][0]
    vmax = vminmax["posteriorMprior"][1]
    ax_pMp, cs_pMp = objBench.spatial_plot(posteriorMprior, LON, LAT, gs = gs[1, 1], fig = fig, cmap = cmap, vmin = vmin, vmax = vmax, **pcolormesh_kwargs)

    if plotTruth:
        vmin = vminmax["truthEmiss"][0]
        vmax = vminmax["truthEmiss"][1]
        ax_truth, cs_truth = objBench.spatial_plot(posteriorMprior, LON, LAT, gs = gs[2, 0], fig = fig, cmap = cmap, vmin = vmin, vmax = vmax, **pcolormesh_kwargs)

        vmin = vminmax["priorMtruth"][0]
        vmax = vminmax["priorMtruth"][1]
        ax_pMt, cs_pMt = objBench.spatial_plot(priorMtruth, LON, LAT, gs = gs[2, 1], fig = fig, cmap = cmap, vmin = vmin, vmax = vmax, **pcolormesh_kwargs)


    ax_prior.set_title("prior")
    ax_posterior.set_title("posterior")
    ax_reduced.set_title("reduced uncertainty")
    ax_pMp.set_title("posterior - prior")

    if plotTruth:
        ax_truth.set_title("truth")
        ax_pMt.set_title("prior - truth")

    cbar_prior = fig.colorbar(cs_prior, ax = ax_prior, **cbar_kwargs)
    cbar_prior.set_label("umol/m^2/s")
    cbar_posterior = fig.colorbar(cs_posterior, ax = ax_posterior, **cbar_kwargs)
    cbar_posterior.set_label("umol/m^2/s")
    cbar_reduced = fig.colorbar(cs_reduced, ax = ax_reduced, **cbar_kwargs)
    cbar_reduced.set_label("(umol/m^2/s)^2")
    cbar_pMp = fig.colorbar(cs_pMp, ax = ax_pMp, **cbar_kwargs)
    cbar_pMp.set_label("umol/m^2/s")
    if plotTruth:
        cbar_truth = fig.colorbar(cs_truth, ax = ax_truth, **cbar_kwargs)
        cbar_truth.set_label("umol/m^2/s")
        cbar_pMt = fig.colorbar(cs_pMt, ax = ax_pMt, **cbar_kwargs)
        cbar_pMt.set_label("umol/m^2/s")

    fig.suptitle(suptitle)

    return fig

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

def draw_snap(objBench, args):

    time, movieStart, movieEnd = args

    
    movie_kwargs = objBench.myConfig["movie_kwargs"]
    snapDir = movie_kwargs["snapDir"]

    print("Movie snap: " + time.strftime("%Y-%m-%d %H:%M:%S (UTC)"))
    objAna = objBench.objAna
    isOSSE = objAna.objCase.objExp.isOSSE
    priorEmiss, LON, LAT = objAna.get_prior_emiss(time)
    posteriorEmiss, _, _ = objAna.get_posterior_emiss(time)
    reducedUncertainty, _, _ = objAna.get_reduced_uncertainty_time(time)
    posteriorMprior = posteriorEmiss - priorEmiss
    if isOSSE:
        truthEmiss, _, _ = objAna.get_truth_emiss(time)
        priorMtruth = priorEmiss - truthEmiss
    else:
        truthEmiss = priorMtruth = None

    fig = plt.figure(figsize = movie_kwargs["figsize"])
    progress_bar(time, movieStart, movieEnd, fig, axRec = [0.1, 0.02, 0.85, 0.18])

    suptitle = time.strftime("%Y-%m-%d %H:%M:%S (UTC)")

    fig = objBench.bis_spatial_plot(LON, LAT, priorEmiss, posteriorEmiss, reducedUncertainty, posteriorMprior, truthEmiss, priorMtruth, plotTruth = isOSSE, fig = fig, suptitle = suptitle, gs_kwargs = {"bottom": 0.15}, **movie_kwargs)

    fig.savefig(snapDir + "/SpatialMovie_" + time.strftime("%Y-%m-%d_%H:%M:%S") + ".png", dpi = movie_kwargs["dpi"])

def make_cartoon(objBench, movieStart, movieEnd, movieDtHrs):

    movie_kwargs = objBench.myConfig["movie_kwargs"]
    movieDir = movie_kwargs["movieDir"]
    movieName = movie_kwargs["movieName"]
    snapDir = movie_kwargs["snapDir"]
    dt = dtm.timedelta(hours = movieDtHrs)

    if not os.path.exists(snapDir):
        os.makedirs(snapDir)

    parallelArgs = []
    current = movieStart
    while(current <= movieEnd):
        parallelArgs.append((current, movieStart, movieEnd))
        current += dt

    pool = mtp.Pool(objBench.nProc)
    pool.map(objBench.draw_snap, parallelArgs)
    pool.close()
    pool.join()
    
    

    if os.path.exists(movieDir + "/" + movieName + ".mp4"):
        print("Deleting " + movieDir + "/" + movieName + ".mp4")
        os.system("rm -f " + movieDir + "/" + movieName + ".mp4")
    #cmd = "ffmpeg -framerate 20 -pattern_type glob -i \"" + snapDir + "/spatial_*.png\" -c:v libx264 -preset ultrafast -pix_fmt yuv420p -r 10 " + movieDir + "/" + movieName + ".mp4"
    cmd = "ffmpeg -framerate 20 -pattern_type glob -i \"" + snapDir + "/SpatialMovie_*.png\" -pix_fmt yuv420p -r 10 " + movieDir + "/" + movieName + ".mp4"
    os.system(cmd)
    #os.system("rm -f " + snapDir + "/SpatialMovie_*_.png")
        
    


#elf.objAna.objCase.objExp.benchConfig
