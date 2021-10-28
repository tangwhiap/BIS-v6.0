#!/usr/bin/env python
# Authors:
#   Wenhan TANG - 05/2021
#   ...

import numpy as np
import datetime as dtm
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
from pdb import set_trace

#from OSSE_Analyst import OSSE_Analyst_zero_one as Analyst
#from OSSE_Analyst import OSSE_Analyst_const_meic as Analyst
from REAL_Analyst import REAL_Analyst_Independent as Analyst
from Bench import BenchPlots, multicases_subplot_scheme

#str_StartTime = "2019-01-01_08:00:00"
#str_EndTime = "2019-01-05_02:00:00"

str_StartTime = "2020-08-01_08:00:00"
str_EndTime = "2020-08-30_08:00:00"

#str_StartTime = "2021-01-03_00:00:00"
#str_EndTime = "2021-01-05_18:00:00"
dtHrs = 6
#CaseName = "Exp_OSSE=zero_one_Lt=48_Ls=100"
#CaseName = "Exp_OSSE=zero_one_Lt=48_Ls=30"
#CaseName = "OSSE_const_meic_t"
#CaseName = "test_high"
#CaseName = "bhd_4"
CaseName = "BTH_d02_test5"
OutDir = "/home/tangwh/public_html/BIS_show/" + CaseName

StartTime = dtm.datetime.strptime(str_StartTime, "%Y-%m-%d_%H:%M:%S")
EndTime = dtm.datetime.strptime(str_EndTime, "%Y-%m-%d_%H:%M:%S")
dt = dtm.timedelta(hours = dtHrs)

if not(os.path.exists(OutDir)):
    os.makedirs(OutDir)

bcp1 = BenchPlots(save_name = OutDir + "/Bench1")
#bcp2 = BenchPlots(save_name = OutDir + "/Bench2", figsize = (12, 12))

case = Analyst(CaseName = CaseName, window = 24)#, Start = StartTime, End = EndTime)

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
        #print(cs, "aaa", vmin_, "bbb", vmax_)
        vmin = vmin_ if vmin_ < vmin else vmin
        vmax = vmax_ if vmax_ > vmax else vmax
    return vmin, vmax

def draw_etot_tlist():
    global case, bcp1
    tlist, prior_tot = case.total_emiss_timelist(emiss_type = "prior")
    _, posterior_tot = case.total_emiss_timelist(emiss_type = "posterior")
    #_, truth_tot = case.total_emiss_timelist(get_emiss_fun = "get_truth_emiss")
    ax = bcp1.subplot(1, 1, 1)
    case.plot_tlist(tlist, prior_tot, ax = ax, label = "prior", linestyle = "--", linewidth = 3, color = "k")
    case.plot_tlist(tlist, posterior_tot, ax = ax, label = "posterior")
    #case.plot_tlist(tlist, truth_tot, ax = ax, label = "truth", linestyle = "--", linewidth = 3, color = "r")
    ax.legend()
    bcp1.cfig.suptitle(case.CaseName)
    bcp1.next_page()


def draw_pptr(time):
    global case, bcp1
    vmin = -40
    vmax = 40
    prior = case.get_prior_emiss(time)
    #truth = case.get_truth_emiss(time)
    posterior = case.get_posterior_emiss(time)
    RQ, _ = case.get_reduced_total_unc(time)
    ax_prior, cs_prior = case.plot2D(prior, fig = bcp1.cfig, subplot = (2,2,1), vmin = vmin, vmax = vmax)
    #ax_truth, cs_truth = case.plot2D(truth, fig = bcp1.cfig, subplot = (3,2,2), vmin = vmin, vmax = vmax)
    ax_posterior, cs_posterior = case.plot2D(posterior, fig = bcp1.cfig, subplot = (2,2,2), vmin = vmin, vmax = vmax)
    ax_RQ, cs_RQ = case.plot2D(RQ, fig = bcp1.cfig, subplot = (2,2,3), vmin = vmin, vmax = vmax)
    ax_pmp, cs_pmp = case.plot2D(posterior - prior, fig = bcp1.cfig, subplot = (2,2,4), vmin = vmin, vmax = vmax)
    #ax_pmt, cs_pmt = case.plot2D(posterior - truth, fig = bcp1.cfig, subplot = (3,2,6), vmin = vmin, vmax = vmax)

    ax_prior.set_title("Prior")
    #ax_truth.set_title("Truth")
    ax_posterior.set_title("Posterior")
    ax_RQ.set_title("Reduced Uncertainty")
    ax_pmp.set_title("Posterior-Prior")
    #ax_pmt.set_title("Posterior-Truth")
    draw_cbar(cs = cs_pmp, bcp = bcp1, loc = (0.15, 0.05, 0.7, 0.03), label = "umol/m^2/s", orientation = "horizontal")
    bcp1.cfig.suptitle(time.strftime("%Y-%m-%d %H"))
    bcp1.next_page()

def opti_prog(time):
    global case, bcp2
    fig_row = 4; fig_col = 3
    def things_after_page(cs_obs_list, cs_foot_list, cs_RQ_list):
        assert len(cs_obs_list) > 0 and len(cs_foot_list) > 0 and len(cs_RQ_list) > 0
        vmin_obs, vmax_obs = get_clim(cs_obs_list)
        print("vmin = ", vmin_obs, "vmax = ", vmax_obs)
        for ics in range(len(cs_obs_list)):
            cs_obs_list[ics].set_clim(vmin_obs, vmax_obs)
        vmin_foot, vmax_foot = get_clim(cs_foot_list)
        for ics in range(len(cs_foot_list)):
            cs_foot_list[ics].set_clim(vmin_foot, vmax_foot)
        vmin_RQ, vmax_RQ = get_clim(cs_RQ_list)
        for ics in range(len(cs_RQ_list)):
            cs_RQ_list[ics].set_clim(vmin_RQ, vmax_RQ)
        cbar_y = 0.04
        cbar_length = 0.3
        cbar_height = 0.03
        cbar_center = [(1/fig_col) * (0.5 + icol) for icol in range(fig_col)]
        cbar_x = [icenter - cbar_length / 2 for icenter in cbar_center]
        print(cbar_x)
        """ 
        cax_cbar_obs = bcp2.vpage(cbar_x[0], cbar_y, cbar_length, cbar_height)
        cbar_obs = bcp2.cfig.colorbar(cs_obs_list[0], cax = cax_cbar_obs)
        cbar_obs.label("ppm")
        cax_cbar_foot = bcp2.vpage(cbar_x[1], cbar_y, cbar_length, cbar_height)
        cbar_foot = bcp2.cfig.colorbar(cs_foot_list[0], cax = cax_cbar_foot)
        cbar_foot.label("ppm m^2 s umol^-1")
        cax_cbar_RQ = bcp2.vpage(cbar_x[2], cbar_y, cbar_length, cbar_height)
        cbar_RQ = bcp2.cfig.colorbar(cs_RQ_list[0], cax = cax_cbar_RQ)
        cbar_RQ.label("umol m^2 s^-1")
        """
        draw_cbar(cs = cs_obs_list[0], bcp = bcp2, loc = (cbar_x[0], cbar_y, cbar_length, cbar_height), orientation = "horizontal", label = "ppm")
        draw_cbar(cs = cs_foot_list[0], bcp = bcp2, loc = (cbar_x[1], cbar_y, cbar_length, cbar_height), orientation = "horizontal", label = "log10(ppm m^2 s umol^-1)")
        draw_cbar(cs = cs_RQ_list[0], bcp = bcp2, loc = (cbar_x[2], cbar_y, cbar_length, cbar_height), orientation = "horizontal", label = "umol m^2 s^-1")
        bcp2.cfig.suptitle(time.strftime("%Y-%m-%d %H"))
        

    mcss = multicases_subplot_scheme(fig_row, fig_col, span_dim = "row", remain_adjust = False, auto_reverse = False)
    opti_time_list = case.opti_recepTime_4eachTime[time.strftime("%Y%m%d%H")]
    N_reduced = len(opti_time_list)

    mcss.layout(N_reduced * 3)

    cs_obs_list = []
    cs_foot_list = []
    cs_RQ_list = []
    for i_reduced in range(N_reduced):
        recepTime = opti_time_list[i_reduced]

        # Draw obs value
        isNextPage, sub1, sub2, sub3 = mcss.get_subloc()
        if isNextPage:
            things_after_page(cs_obs_list, cs_foot_list, cs_RQ_list)
            bcp2.next_page()
            cs_obs_list = []
            cs_foot_list = []
            cs_RQ_list = []

        print("OBS", recepTime)
        #ax_obs = bcp2.subplot(sub1, sub2, sub3)
        print(sub1, sub2, sub3)
        ax_obs, cs_obs, _ = case.draw_OBS_scatter(time = recepTime, fig = bcp2.cfig, subplot = (sub1, sub2, sub3), size = 50, draw_cbar = False)
        #set_trace()
        ax_obs.set_title("OBS " + recepTime.strftime("%Y-%m-%d %H"))
        cs_obs_list.append(cs_obs)
        # Draw footprint
        isNextPage, sub1, sub2, sub3 = mcss.get_subloc()
        footprint = case.get_footprint(time_list = [time], RecepTime_list = [recepTime])
        print(sub1, sub2, sub3)
        ax_foot, cs_foot = case.plot2D(np.log10(footprint), fig = bcp2.cfig, subplot = (sub1, sub2, sub3))
        cs_foot_list.append(cs_foot)
        ax_foot.set_title("Footprint")

       
        # Draw reduced uncertainty 
        isNextPage, sub1, sub2, sub3 = mcss.get_subloc()
        RQ = case.get_reduced_unc(time = time, RecepTime = recepTime)[0]["RQ"]
        print(sub1, sub2, sub3)
        ax_RQ, cs_RQ = case.plot2D(RQ, fig = bcp2.cfig, subplot =  (sub1, sub2, sub3))
        cs_RQ_list.append(cs_RQ)
        ax_RQ.set_title("Reduced uncertainty")

    things_after_page(cs_obs_list, cs_foot_list, cs_RQ_list)
    bcp2.next_page()
        

    #case.draw_OBS_scatter(time = time, fig = bcp2.cfig, subplot = ())
    
if __name__ == "__main__":
    draw_etot_tlist()
    Current = StartTime
    while(Current <= EndTime):
        #Current = dtm.datetime(2021,1,5,18)
        #EndTime = Current
        print(Current.strftime("%Y-%m-%d_%H:%M:%S"))
        draw_pptr(Current)
        #opti_prog(Current)
        Current += dt
    bcp1.print(); bcp1.close()
    #bcp2.print(); bcp2.close()
