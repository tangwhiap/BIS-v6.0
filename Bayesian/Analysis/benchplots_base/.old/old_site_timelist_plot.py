#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...

from .utils.Bench import BenchPlots_png as BenchPlots

import numpy as np
import matplotlib.pyplot as plt

def plot_concentrations(timeList, obsConc, initConc, priorConc, postConc, finalConc, ax, **kwargs):
    #myConfig = objAna.myConfig["site_timeSeries_kwargs"]

    ax.plot(timeList, obsConc, label = "OBS", **kwargs["obsConc_kwargs"])
    ax.plot(timeList, initConc, label = "prior", **kwargs["initConc_kwargs"])
    ax.plot(timeList, postConc, label = "post_primary", **kwargs["postConc_kwargs"])
    ax.plot(timeList, finalConc, label = "post_final", **kwargs["finalConc_kwargs"])

    ax.set_ylim(kwargs["vrange"])
    ax.set_ylabel("ppm")
    ax.grid()
    ax.legend()
    return ax

def plot_Hx(timeList, initHx, priorHx, postHx, finalHx, ax, **kwargs):

    ax.plot(timeList, initHx, label = "Hx_prior", **kwargs["priorHx_kwargs"])
    ax.plot(timeList, postHx, label = "Hx_post", **kwargs["postHx_kwargs"])
    ax.plot(timeList, finalHx, label = "Hx_final", **kwargs["finalHx_kwargs"])
    ax.plot(timeList, np.full((len(timeList)), 0), color = "black", linestyle = "--")

    ax.set_ylim(kwargs["vrange"])
    ax.set_ylabel("ppm")
    ax.grid()
    ax.legend()
    return ax
    
def plot_bck(timeList, bckPrior, sbckPrior, bckProc, sbckProc, ax, **kwargs):

    ax.plot(timeList, bckPrior, label = "BCK_orig", **kwargs["bckPrior_kwargs"])
    ax.fill_between(timeList, bckPrior - sbckPrior, bckPrior + sbckPrior, **kwargs["sbckPrior_kwargs"])

    ax.plot(timeList, bckProc, label = "BCK_optimized", **kwargs["bckProc_kwargs"])
    ax.fill_between(timeList, bckProc - sbckProc, bckProc + sbckProc, **kwargs["sbckProc_kwargs"])

    ax.set_ylim(kwargs["vrange"])
    ax.set_ylabel("ppm")
    ax.grid()
    ax.legend()
    return ax

def plot_point_emiss(timeList, priorEmiss, priorSigma, posteriorEmiss, posteriorSigma, ax, isOSSE, truthEmiss = None, **kwargs):

    ax.plot(timeList, priorEmiss, label = "Prior", **kwargs["priorEmiss_kwargs"])
    ax.fill_between(timeList, priorEmiss - priorSigma, priorEmiss + priorSigma, **kwargs["priorSigma_kwargs"])

    ax.plot(timeList, posteriorEmiss, label = "Post", **kwargs["posteriorEmiss_kwargs"])
    ax.fill_between(timeList, posteriorEmiss - posteriorSigma, posteriorEmiss + posteriorSigma, **kwargs["posteriorSigma_kwargs"])

    if isOSSE:
        ax.plot(timeList, truthEmiss, label = "Truth", **kwargs["truthEmiss_kwargs"]) 

    ax.set_ylim(kwargs["vrange"])
    ax.set_ylabel("umol/m^2/s")
    ax.grid()
    ax.legend()
    return ax

def point_timeSeries(objBench, drawStart, drawEnd):

    objAna = objBench.objAna

    myConfig = objBench.myConfig["point_timeSeries_kwargs"]

    ds_obs = objAna.get_obs_timelist(drawStart, drawEnd)

    ds_initHx = objAna.get_initHx_timelist(drawStart, drawEnd)
    ds_initHxBck = objAna.get_initHxBck_timelist(drawStart, drawEnd)

    ds_priorHx = objAna.get_priorHx_timelist(drawStart, drawEnd)
    ds_procHx = objAna.get_procHx_timelist(drawStart, drawEnd)

    ds_priorHxBck = objAna.get_priorHxBck_timelist(drawStart, drawEnd)
    ds_procHxBck = objAna.get_procHxBck_timelist(drawStart, drawEnd)

    ds_finalHx = objAna.get_finalHx_timelist(drawStart, drawEnd)
    ds_finalHxBck = objAna.get_finalHxBck_timelist(drawStart, drawEnd)

    ds_bckPrior = objAna.get_bckPrior_timelist(drawStart, drawEnd)
    ds_bckProc = objAna.get_bckProc_timelist(drawStart, drawEnd)

    ds_sbckPrior = objAna.get_sbckPrior_timelist(drawStart, drawEnd)
    ds_sbckProc = objAna.get_sbckProc_timelist(drawStart, drawEnd)

    ds_priorEmiss = objAna.get_priorEmissPoints_timelist(drawStart, drawEnd)
    ds_priorSigma = objAna.get_priorSigmaPoints_timelist(drawStart, drawEnd)

    ds_posteriorEmiss = objAna.get_posteriorEmissPoints_timelist(drawStart, drawEnd)
    ds_posteriorSigma = objAna.get_posteriorSigmaPoints_timelist(drawStart, drawEnd)
    
    if objAna.isOSSE:
        ds_truthEmiss = objAna.get_truthEmissPoints_timelist(drawStart, drawEnd)

    bcp = BenchPlots(figsize = myConfig["figsize"], save_name = myConfig["BenchName"])

    for siteName in objAna.obsLoc:
        #gs = bcp.cfig.add_gridspec(ncols = 1, nrows = 4)
        axs = bcp.cfig.subplots(4,1,sharex =  True)
        bcp.cfig.subplots_adjust(hspace = 0.1) 

        ax = axs[0]
        timeList = ds_initHx["time"].values
        obsConc = ds_obs.sel(site = siteName)["data"].values
        initConc = ds_initHxBck.sel(site = siteName)["data"].values
        priorConc = ds_priorHxBck.sel(site = siteName)["data"].values
        postConc = ds_procHxBck.sel(site = siteName)["data"].values
        finalConc = ds_finalHxBck.sel(site = siteName)["data"].values
        plot_concentrations(timeList, obsConc, initConc, priorConc, postConc, finalConc, ax, **myConfig["plot_concentration_kwargs"])

        ax = axs[1] 
        initHx = ds_initHx.sel(site = siteName)["data"].values
        priorHx = ds_priorHx.sel(site = siteName)["data"].values
        postHx = ds_procHx.sel(site = siteName)["data"].values
        finalHx = ds_finalHx.sel(site = siteName)["data"].values
        plot_Hx(timeList, initHx, priorHx, postHx, finalHx, ax, **myConfig["plot_Hx_kwargs"])

        ax = axs[2]
        bckPrior = ds_bckPrior.sel(site = siteName)["data"].values
        sbckPrior = ds_sbckPrior.sel(site = siteName)["data"].values
        bckProc = ds_bckProc.sel(site = siteName)["data"].values
        sbckProc = ds_sbckProc.sel(site = siteName)["data"].values
        plot_bck(timeList, bckPrior, sbckPrior, bckProc, sbckProc, ax, **myConfig["plot_bck_kwargs"])

        ax = axs[3]
        priorEmiss = ds_priorEmiss.sel(site = siteName)["data"].values
        priorSigma = ds_priorSigma.sel(site = siteName)["data"].values
        posteriorEmiss = ds_posteriorEmiss.sel(site = siteName)["data"].values
        posteriorSigma = ds_posteriorSigma.sel(site = siteName)["data"].values

        isOSSE = objAna.isOSSE
        if isOSSE:
            truthEmiss = ds_truthEmiss.sel(site = siteName)["data"].values
        else:
            truthEmiss = None

        plot_point_emiss(timeList, priorEmiss, priorSigma, posteriorEmiss, posteriorSigma, ax, isOSSE, truthEmiss, **myConfig["plot_point_emiss_kwargs"])

        bcp.cfig.autofmt_xdate()
        bcp.cfig.suptitle(siteName)

        bcp.next_page()

    bcp.close()