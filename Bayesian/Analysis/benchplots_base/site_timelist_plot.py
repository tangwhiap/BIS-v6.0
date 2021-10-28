#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...

from .draw_object_class import AnalysisTimeSeries, AnalysisLinePlots, ds2ATS
from .spatial_plot import stations_local_plot
from ..utils.Bench import BenchPlots_png as BenchPlots

from pdb import set_trace

def site_timeSeries(objBench, drawStart, drawEnd):

    
    """
    if objAna.isOSSE:
        ds_truthEmiss = objAna.get_truthEmissPoints_timelist(drawStart, drawEnd)
    """

    myConfig = objBench.myConfig["site_timeSeries_kwargs"]
    staPlotConfig = myConfig["station_local_kwargs"]
    obsConfig = myConfig["obs_timeSeries_kwargs"]

    bcp = BenchPlots(figsize = myConfig["figsize"], save_name = objBench.benchDir + "/" + myConfig["BenchName"])

    sites_observation_dic = objBench.build_sites_observation_dic(drawStart = drawStart, drawEnd = drawEnd)
    bcp.page_reverse()
    gs = (bcp.cfig.add_gridspec(ncols = 1, nrows = 2))[0]
    staLon = [objBench.objAna.obsLoc[staName]["lon"] for staName in objBench.objAna.obsLoc]
    staLat = [objBench.objAna.obsLoc[staName]["lat"] for staName in objBench.objAna.obsLoc]
    stations_local_plot(objBench, staLon, staLat, gs, setting = staPlotConfig, fig = bcp.cfig)



    ax = bcp.cfig.add_subplot(2, 1, 2)
    objALP = AnalysisLinePlots(sites_observation_dic, LinePlots_kwargs = obsConfig)
    ax = objALP.plot(ax = ax)

    bcp.cfig.autofmt_xdate()
    bcp.cfig.suptitle("observation")

    bcp.next_page()
    bcp.page_reverse()

    for siteName in objBench.objAna.obsLoc:

        print("drawing site: " + siteName)

        timeSeries_dic = objBench.build_site_timeSeries_dic(siteName = siteName, drawStart = drawStart, drawEnd = drawEnd)

        axs = bcp.cfig.subplots(len(timeSeries_dic), 1, sharex =  True)
        bcp.cfig.subplots_adjust(hspace = 0.1) 

        for ax,  tsName in zip(axs, timeSeries_dic):
            LinePlots_kwargs = myConfig[tsName + "_kwargs"]
            objALP = AnalysisLinePlots(timeSeries_dic[tsName], LinePlots_kwargs = LinePlots_kwargs)
            ax = objALP.plot(ax = ax)

        bcp.cfig.autofmt_xdate()
        bcp.cfig.suptitle(siteName)
        bcp.next_page()

    bcp.close()

def build_site_timeSeries_dataset(objBench, drawStart, drawEnd):

    objBench.site_timeSeries_dataset = {}

    objAna = objBench.objAna

    objBench.site_timeSeries_dataset["ds_obs"] = objAna.get_obs_timelist(drawStart, drawEnd)

    objBench.site_timeSeries_dataset["ds_initHx"] = objAna.get_initHx_timelist(drawStart, drawEnd)
    objBench.site_timeSeries_dataset["ds_initHxBck"] = objAna.get_initHxBck_timelist(drawStart, drawEnd)

    objBench.site_timeSeries_dataset["ds_priorHx"] = objAna.get_priorHx_timelist(drawStart, drawEnd)
    objBench.site_timeSeries_dataset["ds_procHx"] = objAna.get_procHx_timelist(drawStart, drawEnd)

    objBench.site_timeSeries_dataset["ds_priorHxBck"] = objAna.get_priorHxBck_timelist(drawStart, drawEnd)
    objBench.site_timeSeries_dataset["ds_procHxBck"] = objAna.get_procHxBck_timelist(drawStart, drawEnd)

    objBench.site_timeSeries_dataset["ds_finalHx"] = objAna.get_finalHx_timelist(drawStart, drawEnd)
    objBench.site_timeSeries_dataset["ds_finalHxBck"] = objAna.get_finalHxBck_timelist(drawStart, drawEnd)

    
    objBench.site_timeSeries_dataset["ds_bckPrior"] = objAna.get_bckPrior_timelist(drawStart, drawEnd)
    objBench.site_timeSeries_dataset["ds_bckProc"] = objAna.get_bckProc_timelist(drawStart, drawEnd)

    objBench.site_timeSeries_dataset["ds_sbckPrior"] = objAna.get_sbckPrior_timelist(drawStart, drawEnd)
    objBench.site_timeSeries_dataset["ds_sbckProc"] = objAna.get_sbckProc_timelist(drawStart, drawEnd)

    
    objBench.site_timeSeries_dataset["ds_priorEmiss"] = objAna.get_priorEmissPoints_timelist(drawStart, drawEnd)
    objBench.site_timeSeries_dataset["ds_priorSigma"] = objAna.get_priorSigmaPoints_timelist(drawStart, drawEnd)

    objBench.site_timeSeries_dataset["ds_posteriorEmiss"] = objAna.get_posteriorEmissPoints_timelist(drawStart, drawEnd)
    objBench.site_timeSeries_dataset["ds_posteriorSigma"] = objAna.get_posteriorSigmaPoints_timelist(drawStart, drawEnd)

    objBench.site_timeSeries_loaded = True

def build_sites_observation_dic(objBench, drawStart, drawEnd):
    if not objBench.site_timeSeries_loaded:
        objBench.build_site_timeSeries_dataset(drawStart, drawEnd)

    ds_obs = objBench.site_timeSeries_dataset["ds_obs"]
    sites_nameList = ds_obs.site.values
    sites_observation_dic = {}

    for siteName in sites_nameList:
        sites_observation_dic[siteName] = ds2ATS(ds_line = ds_obs.sel(site = siteName), label = siteName)
    
    return sites_observation_dic


def build_site_timeSeries_dic(objBench, siteName, drawStart, drawEnd):

    #myConfig = objBench.myConfig["site_timeSeries_kwargs"]
    if not objBench.site_timeSeries_loaded:
        objBench.build_site_timeSeries_dataset(drawStart, drawEnd)

    site_timeSeries_dic = {}
 
    ds_obs = objBench.site_timeSeries_dataset["ds_obs"].sel(site = siteName)
    ds_initHxBck = objBench.site_timeSeries_dataset["ds_initHxBck"].sel(site = siteName)
    ds_procHxBck = objBench.site_timeSeries_dataset["ds_procHxBck"].sel(site = siteName)
    ds_finalHxBck = objBench.site_timeSeries_dataset["ds_finalHxBck"].sel(site = siteName)
    

    site_timeSeries_dic["plot_concentration"] = {
        "obsConc": ds2ATS(ds_line = ds_obs, label = "OBS"),
        "initConc": ds2ATS(ds_line = ds_initHxBck, label = "prior"),
        #"priorConc": ds2ATS(ds_line = ds_priorHxBck, label = "Prior"),
        "postConc": ds2ATS(ds_line = ds_procHxBck, label = "post_primary"),
        "finalConc": ds2ATS(ds_line = ds_finalHxBck, label = "post_final"),
    }

    ds_priorHx = objBench.site_timeSeries_dataset["ds_priorHx"].sel(site = siteName)
    ds_procHx = objBench.site_timeSeries_dataset["ds_procHx"].sel(site = siteName)
    ds_finalHx = objBench.site_timeSeries_dataset["ds_finalHx"].sel(site = siteName)

    site_timeSeries_dic["plot_Hx"] = {
        "priorHx": ds2ATS(ds_line = ds_priorHx, label = "Hx_prior"),
        "postHx": ds2ATS(ds_line = ds_procHx, label = "Hx_post"),
        "finalHx": ds2ATS(ds_line = ds_finalHx, label = "Hx_final"),
    }

    ds_bckPrior = objBench.site_timeSeries_dataset["ds_bckPrior"].sel(site = siteName)
    ds_sbckPrior = objBench.site_timeSeries_dataset["ds_sbckPrior"].sel(site = siteName)
    ds_bckProc = objBench.site_timeSeries_dataset["ds_bckProc"].sel(site = siteName)
    ds_sbckProc = objBench.site_timeSeries_dataset["ds_sbckProc"].sel(site = siteName)

    site_timeSeries_dic["plot_bck"] = {
        "bckPrior": ds2ATS(ds_line = ds_bckPrior, ds_spread = ds_sbckPrior, label = "BCK_orig"),
        "bckProc": ds2ATS(ds_line = ds_bckProc, ds_spread = ds_sbckProc, label = "BCK_optimized"),
    }

    ds_priorEmiss = objBench.site_timeSeries_dataset["ds_priorEmiss"].sel(site = siteName)
    ds_priorSigma = objBench.site_timeSeries_dataset["ds_priorSigma"].sel(site = siteName)
    ds_posteriorEmiss = objBench.site_timeSeries_dataset["ds_posteriorEmiss"].sel(site = siteName)
    ds_posteriorSigma = objBench.site_timeSeries_dataset["ds_posteriorSigma"].sel(site = siteName)

    site_timeSeries_dic["plot_point_emiss"] = {
        "priorEmiss": ds2ATS(ds_line = ds_priorEmiss, ds_spread = ds_priorSigma, label = "Prior"),
        "posteriorEmiss": ds2ATS(ds_line = ds_posteriorEmiss, ds_spread = ds_posteriorSigma, label = "Posterior"),
    }

    return site_timeSeries_dic

