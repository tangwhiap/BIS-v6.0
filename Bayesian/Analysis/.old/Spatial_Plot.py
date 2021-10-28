#!/usr/bin/env python
# Benchmark of spatial plotting.
# Authors:
#   Wenhan TANG - 06/2021
#   ...

import numpy as np
import datetime as dtm
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
from pdb import set_trace
from .utils.Bench import BenchPlots, multicases_subplot_scheme

def draw_cbar(cs, bcp, loc = (0.2, 0.2, 0.7, 0.1), label = None, **kwargs):
    ax_cbar = bcp.vpage(loc[0], loc[1], loc[2], loc[3])
    cbar = bcp.cfig.colorbar(cs, cax = ax_cbar, **kwargs)
    if label is not None:
        cbar.set_label(label)

def get_clim(cs_list):
    vmin = np.inf
    vmax = -np.inf
    for cs in cs_list:
        vmin_, vmax_ = cs.get_clim()
        vmin = vmin_ if vmin_ < vmin else vmin
        vmax = vmax_ if vmax_ > vmax else vmax
    return vmin, vmax

def draw_spatial(time, case, bcp, CaseType = "REAL", vmin = -40, vmax = 40):

    CaseType = CaseType.upper()
    assert CaseType in ["OSSE", "REAL"]
    isOSSE = CaseType == "OSSE"
    cmap = "coolwarm"
    #global case, bcp
    prior = case.get_prior_emiss(time)
    posterior = case.get_posterior_emiss(time)
    RQ, _ = case.get_reduced_total_unc(time)

    if isOSSE:
        truth = case.get_truth_emiss(time)

    vp1 = 3 if isOSSE else 2
    vp2 = 2
    ax_prior, cs_prior = case.plot2D(prior, fig = bcp.cfig, subplot = (vp1,vp2,1), vmin = vmin, vmax = vmax, cmap = cmap)
    ax_posterior, cs_posterior = case.plot2D(posterior, fig = bcp.cfig, subplot = (vp1,vp2,2), vmin = vmin, vmax = vmax, cmap = cmap)
    ax_RQ, cs_RQ = case.plot2D(RQ, fig = bcp.cfig, subplot = (vp1,vp2,3), vmin = vmin, vmax = vmax, cmap = cmap)
    ax_pmp, cs_pmp = case.plot2D(posterior - prior, fig = bcp.cfig, subplot = (vp1,vp2,4), vmin = vmin, vmax = vmax, cmap = cmap)

    if isOSSE:
        ax_truth, cs_truth = case.plot2D(truth, fig = bcp.cfig, subplot = (vp1,vp2,5), vmin = vmin, vmax = vmax, cmap = cmap)
        ax_pmt, cs_pmt = case.plot2D(posterior - truth, fig = bcp.cfig, subplot = (vp1,vp2,6), vmin = vmin, vmax = vmax, cmap = cmap)

    ax_prior.set_title("Prior")
    ax_posterior.set_title("Posterior")
    ax_RQ.set_title("Reduced Uncertainty")
    ax_pmp.set_title("Posterior-Prior")

    if isOSSE:
        ax_truth.set_title("Truth")
        ax_pmt.set_title("Posterior-Truth")

    draw_cbar(cs = cs_pmp, bcp = bcp, loc = (0.15, 0.05, 0.7, 0.03), label = "umol/m^2/s", orientation = "horizontal")
    bcp.cfig.suptitle(time.strftime("%Y-%m-%d %H"))
    #bcp.next_page()

def opti_prog(time, case, bcp, vr_obs = None, vr_foot = None, vr_dEmiss = None, vr_RQ = None):

    fig_row = 4; fig_col = 4
    cmap = "jet"
    def things_after_page(cs_obs_list, cs_foot_list, cs_dEmiss_list, cs_RQ_list):
        assert len(cs_obs_list) > 0 and len(cs_foot_list) > 0 and len(cs_RQ_list) > 0 and len(cs_dEmiss_list) > 0
        if vr_obs is None:
            vmin_obs, vmax_obs = get_clim(cs_obs_list)
        else:
            vmin_obs, vmax_obs = vr_obs
        for ics in range(len(cs_obs_list)):
            cs_obs_list[ics].set_clim(vmin_obs, vmax_obs)
        if vr_foot is None:
            vmin_foot, vmax_foot = get_clim(cs_foot_list)
        else:
            vmin_foot, vmax_foot = vr_foot
        for ics in range(len(cs_foot_list)):
            cs_foot_list[ics].set_clim(vmin_foot, vmax_foot)
        if vr_dEmiss is None:
            vmin_dEmiss, vmax_dEmiss = get_clim(cs_dEmiss_list)
        else:
            vmin_dEmiss, vmax_dEmiss = vr_dEmiss
        for ics in range(len(cs_dEmiss_list)):
            cs_dEmiss_list[ics].set_clim(vmin_dEmiss, vmax_dEmiss)
        if vr_RQ is None:
            vmin_RQ, vmax_RQ = get_clim(cs_RQ_list)
        else:
            vmin_RQ, vmax_RQ = vr_RQ
        for ics in range(len(cs_RQ_list)):
            cs_RQ_list[ics].set_clim(vmin_RQ, vmax_RQ)
        
        cbar_y = 0.04
        cbar_length = 0.23
        cbar_height = 0.03
        cbar_center = [(1/fig_col) * (0.5 + icol) for icol in range(fig_col)]
        cbar_x = [icenter - cbar_length / 2 for icenter in cbar_center]
        #print(cbar_x)
        """ 
        cax_cbar_obs = bcp.vpage(cbar_x[0], cbar_y, cbar_length, cbar_height)
        cbar_obs = bcp.cfig.colorbar(cs_obs_list[0], cax = cax_cbar_obs)
        cbar_obs.label("ppm")
        cax_cbar_foot = bcp.vpage(cbar_x[1], cbar_y, cbar_length, cbar_height)
        cbar_foot = bcp.cfig.colorbar(cs_foot_list[0], cax = cax_cbar_foot)
        cbar_foot.label("ppm m^2 s umol^-1")
        cax_cbar_RQ = bcp.vpage(cbar_x[2], cbar_y, cbar_length, cbar_height)
        cbar_RQ = bcp.cfig.colorbar(cs_RQ_list[0], cax = cax_cbar_RQ)
        cbar_RQ.label("umol m^2 s^-1")
        """
        draw_cbar(cs = cs_obs_list[0], bcp = bcp, loc = (cbar_x[0], cbar_y, cbar_length, cbar_height), orientation = "horizontal", label = "ppm")
        draw_cbar(cs = cs_foot_list[0], bcp = bcp, loc = (cbar_x[1], cbar_y, cbar_length, cbar_height), orientation = "horizontal", label = "log10(ppm m^-2 s umol^-1)")
        draw_cbar(cs = cs_dEmiss_list[0], bcp = bcp, loc = (cbar_x[2], cbar_y, cbar_length, cbar_height), orientation = "horizontal", label = "umol m^-2 s^-1")
        draw_cbar(cs = cs_RQ_list[0], bcp = bcp, loc = (cbar_x[3], cbar_y, cbar_length, cbar_height), orientation = "horizontal", label = "log10(umol^2 m-^4 s^-2)")
        bcp.cfig.suptitle(time.strftime("%Y-%m-%d %H"), fontsize = 20)
        

    mcss = multicases_subplot_scheme(fig_row, fig_col, span_dim = "row", remain_adjust = False, auto_reverse = False)
    opti_time_list = case.opti_recepTime_4eachTime[time.strftime("%Y%m%d%H")]
    N_reduced = len(opti_time_list)
    if N_reduced == 0:
        bcp.cfig.text(0.5, 0.5, "No optimized!", fontsize = 80, horizontalalignment = "center", verticalalignment = "center")#, transform = ax.transAxes)
        bcp.cfig.suptitle(time.strftime("%Y-%m-%d %H"), fontsize = 20)
        return

    mcss.layout(N_reduced * 4)

    cs_obs_list = []
    cs_foot_list = []
    cs_RQ_list = []
    cs_dEmiss_list = []
    for i_reduced in range(N_reduced):
        recepTime = opti_time_list[i_reduced]

        # Draw obs value
        isNextPage, sub1, sub2, sub3 = mcss.get_subloc()
        if isNextPage:
            things_after_page(cs_obs_list, cs_foot_list, cs_dEmiss_list, cs_RQ_list)
            bcp.next_page()
            cs_obs_list = []
            cs_foot_list = []
            cs_RQ_list = []
            cs_dEmiss_list = []

        print("RecepTime: ", recepTime)
        #print("OBS", recepTime)
        #ax_obs = bcp.subplot(sub1, sub2, sub3)
        #print(sub1, sub2, sub3)
        ax_obs, cs_obs, _ = case.draw_OBS_scatter(time = recepTime, fig = bcp.cfig, subplot = (sub1, sub2, sub3), size = 50, draw_cbar = False)
        ax_obs.set_title("OBS " + recepTime.strftime("%Y-%m-%d %H"))
        cs_obs_list.append(cs_obs)
        # Draw footprint
        isNextPage, sub1, sub2, sub3 = mcss.get_subloc()
        footprint = case.get_footprint(time_list = [time], RecepTime_list = [recepTime])
        #print(sub1, sub2, sub3)
        ax_foot, cs_foot = case.plot2D(np.log10(footprint), fig = bcp.cfig, subplot = (sub1, sub2, sub3), cmap = cmap)
        cs_foot_list.append(cs_foot)
        ax_foot.set_title("Footprint")

      
        # Draw emiss increment
        isNextPage, sub1, sub2, sub3 = mcss.get_subloc()
        dEmiss = case.get_dEmiss(time = time, RecepTime = recepTime)[0]["dEmiss"]
        #print(sub1, sub2, sub3)
        ax_dEmiss, cs_dEmiss = case.plot2D(dEmiss, fig = bcp.cfig, subplot =  (sub1, sub2, sub3), cmap = cmap)
        cs_dEmiss_list.append(cs_dEmiss)
        ax_dEmiss.set_title("Emiss increment")

        # Draw reduced uncertainty 
        isNextPage, sub1, sub2, sub3 = mcss.get_subloc()
        RQ = case.get_reduced_unc(time = time, RecepTime = recepTime)[0]["RQ"]
        #print(sub1, sub2, sub3)
        ax_RQ, cs_RQ = case.plot2D(np.log10(RQ), fig = bcp.cfig, subplot =  (sub1, sub2, sub3), cmap = cmap)
        cs_RQ_list.append(cs_RQ)
        ax_RQ.set_title("Reduced uncertainty")

    things_after_page(cs_obs_list, cs_foot_list, cs_dEmiss_list, cs_RQ_list)
    #bcp.next_page()

def Bench_draw_spatial(Start, End, dt, case, OutDir, BenchName = "Bench_spatial", CaseType = "REAL", vmin = -40, vmax = 40, isDraw = True):
    if isDraw == False:
        return
    if isinstance(Start, str):
        Start = dtm.datetime.strptime(Start, "%Y-%m-%d_%H:%M:%S")
    if isinstance(End, str):
        End = dtm.datetime.strptime(End, "%Y-%m-%d_%H:%M:%S")
    if isinstance(dt, int):
        dt = dtm.timedelta(hours = dt)
    bcp = BenchPlots(save_name = OutDir + "/" + BenchName)
    Current = Start
    while(Current <= End):
        print("Spatial Plotting:", Current)
        draw_spatial(time = Current, case = case, bcp = bcp, CaseType = CaseType, vmin = vmin, vmax = vmax)
        bcp.next_page()
        Current += dt
    bcp.close()

def Bench_draw_optprog(Start, End, dt, case, OutDir, BenchName = "Bench_optprog", vr_obs = None, vr_foot = None, vr_dEmiss = None, vr_RQ = None, isDraw = True):
    if isDraw == False:
        return
    if isinstance(Start, str):
        Start = dtm.datetime.strptime(Start, "%Y-%m-%d_%H:%M:%S")
    if isinstance(End, str):
        End = dtm.datetime.strptime(End, "%Y-%m-%d_%H:%M:%S")
    if isinstance(dt, int):
        dt = dtm.timedelta(hours = dt)
    bcp = BenchPlots(save_name = OutDir + "/" + BenchName, figsize = (16, 12))
    Current = Start
    while(Current <= End):
        print("Optimizing Plotting:", Current)
        opti_prog(time = Current, case = case, bcp = bcp, vr_obs = vr_obs, vr_foot = vr_foot, vr_dEmiss = vr_dEmiss, vr_RQ = vr_RQ)
        bcp.next_page()
        Current += dt
    bcp.close()
