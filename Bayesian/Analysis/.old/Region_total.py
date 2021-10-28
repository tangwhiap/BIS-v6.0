#!/usr/bin/env python
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

#str_StartTime = "2020-08-01_08:00:00"
#str_EndTime = "2020-08-30_08:00:00"

#str_StartTime = "2021-01-03_00:00:00"
#str_EndTime = "2021-01-05_18:00:00"
#dtHrs = 6
#CaseName = "Exp_OSSE=zero_one_Lt=48_Ls=100"
#CaseName = "Exp_OSSE=zero_one_Lt=48_Ls=30"
#CaseName = "OSSE_const_meic_t"
#CaseName = "test_high"
#CaseName = "bhd_4"
CaseName = "BTH_d02_test5"
OutDir = "/home/tangwh/public_html/BIS_show/" + CaseName

#StartTime = dtm.datetime.strptime(str_StartTime, "%Y-%m-%d_%H:%M:%S")
#EndTime = dtm.datetime.strptime(str_EndTime, "%Y-%m-%d_%H:%M:%S")
#dt = dtm.timedelta(hours = dtHrs)

if not(os.path.exists(OutDir)):
    os.makedirs(OutDir)

#bcp1 = BenchPlots(save_name = OutDir + "/Bench1")
#bcp2 = BenchPlots(save_name = OutDir + "/Bench2", figsize = (12, 12))

bcp = BenchPlots(save_name = OutDir + "/emiss_total")
case = Analyst(CaseName = CaseName, window = 24)#, Start = StartTime, End = EndTime)

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

def draw_region_total_tlist(Region, ax):
    global case, bcp
    tlist, prior_total, prior_sigma_total = case.total_emiss_timelist(emiss_type = "prior", Region = Region, sigma_out = True)
    _, posterior_total, posterior_sigma_total = case.total_emiss_timelist(emiss_type = "posterior", Region = Region, sigma_out = True)
    tlist = tlist[:-1]
    prior_total = prior_total[:-1]
    prior_sigma_total = prior_sigma_total[:-1]
    posterior_total = posterior_total[:-1]
    posterior_sigma_total = posterior_sigma_total[:-1]
    ax.plot(tlist, prior_total, color = "black", label = "Prior")
    ax.fill_between(tlist, prior_total - prior_sigma_total, prior_total + prior_sigma_total, facecolor = "grey", alpha = 0.3)
    ax.plot(tlist, posterior_total, color = "red", label = "Posterior")
    ax.fill_between(tlist, posterior_total - posterior_sigma_total, posterior_total + posterior_sigma_total, facecolor = "red", alpha = 0.3)
    Total_prior = prior_total.sum()
    Sigma_prior = prior_sigma_total.sum()
    Total_posterior = posterior_total.sum()
    Sigma_posterior = posterior_sigma_total.sum()
    ax.grid(color = "grey", linestyle = "--", linewidth = 0.2)
    ax.set_title(Region + " unit: tCO2")
    #ax.set_title(Region + "\nPrior: %.4f ± %.4f Posterior:  %.4f ± %.4f" % (Total_prior, Sigma_prior, Total_posterior, Sigma_posterior))
    ax.text(0.02, 0.95, "Prior: %.2f ± %.2f" % (Total_prior, Sigma_prior), horizontalalignment = "left", verticalalignment = "top", transform = ax.transAxes)
    ax.text(0.02, 0.85, "Posterior: %.2f ± %.2f" % (Total_posterior, Sigma_posterior), horizontalalignment = "left", verticalalignment = "top", transform = ax.transAxes)
    bcp.cfig.autofmt_xdate()


if __name__ == "__main__":
    for region in ["Beijing", "Tianjin"]:
        ax = bcp.subplot(3,1,2)
        draw_region_total_tlist(region, ax)
        bcp.next_page()
    bcp.close()
    #tlist, posterior_tot = case.total_emiss_timelist(emiss_type = "posterior")
    #tlist, posterior_tot_BJ, sigma_total_BJ = case.total_emiss_timelist(emiss_type = "prior", Region = "Beijing", sigma_out = True)
    #tlist, posterior_tot_TJ = case.total_emiss_timelist(emiss_type = "posterior", Region = "Tianjin")
    #plt.plot(tlist, posterior_tot, label = "ALL")
    #plt.plot(tlist, posterior_tot_BJ, label = "BJ")
    #plt.plot(tlist, sigma_total_BJ, label = "BJ_sigma")
    #plt.plot(tlist, posterior_tot_TJ, label = "TJ")
    #plt.legend()
    #plt.show()

