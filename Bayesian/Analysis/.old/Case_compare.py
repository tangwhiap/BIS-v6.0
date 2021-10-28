#!/usr/bin/env python
# Benchplot for multicases comparing.
# Authors:
#    Wenhan TANG - 05/2021 (Original Version)
#    ...

import numpy as np
import matplotlib.pyplot as plt
import datetime as dtm
from Bench import BenchPlots, multicases_subplot_scheme
from OSSE_Analyst import *
from config_CaseCompare import CaseDict as CaseConfig
"""
CaseDict = {
               "casename": [LongName, Start, End, dt, window, n_tback], 
               "Ls200": ["Exp_OSSE=zero_one_Lt=48_Ls=200", "2021-01-03_00:00:00", "2021-01-06_00:00:00", 1, 6, 24]
           }
"""
#print(CaseDict["Ls100"])

str_Analyse_Start = "2021-01-03_00:00:00"
str_Analyse_End = "2021-01-05_18:00:00"
OutDir = "/home/tangwh/public_html/BIS_show/1"


###
Analyse_Start = dtm.datetime.strptime(str_Analyse_Start, "%Y-%m-%d_%H:%M:%S")
Analyse_End = dtm.datetime.strptime(str_Analyse_End, "%Y-%m-%d_%H:%M:%S")
Cases_dic = {}
for case in CaseConfig:
    #Cases_dic[case] = OSSE_Analyst_zero_one(**CaseConfig[case])
    Cases_dic[case] = OSSE_Analyst_const_meic(**CaseConfig[case])
Ncase = len(Cases_dic)
bcp = BenchPlots(save_name = OutDir + "/Ls_compare")


def draw_cbar(cs, loc = (0.2, 0.2, 0.7, 0.1), label = None, **kwargs):
    global bcp
    ax_cbar = bcp.vpage(loc[0], loc[1], loc[2], loc[3])
    cbar = bcp.cfig.colorbar(cs, cax = ax_cbar, **kwargs)
    if label is not None:
        cbar.set_label(label)

# Mean Posterior Comparing
vmin = 0
vmax = 50
mcss =  multicases_subplot_scheme(3, 2, remain_adjust = False)
mcss.layout(Ncase)

for case_name in Cases_dic:
    #print(case_name)
    case = Cases_dic[case_name]
    print(case.CaseDir)
    emiss = case.posterior_mean(start = Analyse_Start, end = Analyse_End)
    print(emiss.mean())
    isNextPage, sbp1, sbp2, sbp3 = mcss.get_subloc()
    if isNextPage:
        
        draw_cbar(cs, loc = (0.15, 0.05, 0.7, 0.03), label = "umol/m^2/s", orientation = "horizontal")    
        bcp.cfig.suptitle("Posterior emiss")
        bcp.next_page()
    ax, cs = case.plot2D(emiss, fig = bcp.cfig, subplot = (sbp1, sbp2, sbp3), vmin = vmin, vmax = vmax)
    ax.set_title(case_name)

draw_cbar(cs, loc = (0.15, 0.05, 0.7, 0.03), label = "umol/m^2/s", orientation = "horizontal")    
bcp.cfig.suptitle("Posterior emiss")
bcp.next_page()

# Mean reduced uncertainty
vmin = 0
vmax = 60
mcss =  multicases_subplot_scheme(3, 2, remain_adjust = False)
mcss.layout(Ncase)
for case_name in Cases_dic:
    print(case_name)
    case = Cases_dic[case_name]
    RQ_mean = case.get_reduced_total_unc_mean(Start = Analyse_Start, End = Analyse_End)
    isNextPage, sbp1, sbp2, sbp3 = mcss.get_subloc()
    if isNextPage:
        draw_cbar(cs, loc = (0.15, 0.05, 0.7, 0.03), label = "umol/m^2/s", orientation = "horizontal")
        bcp.cfig.suptitle("Reduced uncertainty")
        bcp.next_page()
    ax, cs = case.plot2D(RQ_mean, fig = bcp.cfig, subplot = (sbp1, sbp2, sbp3), vmin = vmin, vmax = vmax)
    #bcp.cfig.colorbar(cs)
    ax.set_title(case_name)
draw_cbar(cs, loc = (0.15, 0.05, 0.7, 0.03), label = "umol/m^2/s", orientation = "horizontal") 
bcp.cfig.suptitle("Reduced uncertainty")
bcp.next_page()

# Total emiss time list
ax = bcp.subplot(1, 1, 1)
for case_name in Cases_dic:
    case = Cases_dic[case_name]
    time, etot = case.total_emiss_timelist(emiss_type = "posterior", start = Analyse_Start, end = Analyse_End)
    case.plot_tlist(time, etot, ax = ax, label = case_name)
time, etot_prior = case.total_emiss_timelist(emiss_type = "prior", start = Analyse_Start, end = Analyse_End)
case.plot_tlist(time, etot_prior, ax = ax, label = "Prior", linestyle = "--", linewidth = 3, color = "k")
time, etot_truth = case.total_emiss_timelist(get_emiss_fun = "get_truth_emiss", start = Analyse_Start, end = Analyse_End)
case.plot_tlist(time, etot_truth, ax = ax, label = "Truth", linestyle = "--", linewidth = 3, color = "r")
ax.legend()
ax.set_ylabel("Total emiss (tCO2)")
bcp.next_page()

# Total reduced uncertainty time list
ax = bcp.subplot(1, 1, 1)
for case_name in Cases_dic:
    print(case_name)
    case = Cases_dic[case_name]
    time, RQtot = case.total_RQ_timelist(Start =  Analyse_Start, End = Analyse_End)
    case.plot_tlist(time, RQtot, ax = ax, label = case_name)
ax.legend()
ax.set_ylabel("Total reduced uncertainty of emiss (tCO2)")

bcp.print()
bcp.close()


    


