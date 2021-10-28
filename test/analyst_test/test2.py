#!/usr/bin/env python

from Bayesian.app import case
from Bayesian.Analysis.Analyst.bis_analyst import Analyst

import datetime as dtm
import matplotlib.pyplot as plt

from pdb import set_trace

a = Analyst(case)

start = dtm.datetime(2019, 1, 1)

end = dtm.datetime(2019, 1, 6)

#unc, LON, LAT = a.get_reduced_uncertainty(time, recepTime)
obs = a.get_obs_timelist(start, end)

initHx = a.get_initHx_timelist(start, end)
initHxBck = a.get_initHxBck_timelist(start, end)

priorHx = a.get_priorHx_timelist(start, end)
procHx = a.get_procHx_timelist(start, end)

priorHxBck = a.get_priorHxBck_timelist(start, end)
procHxBck = a.get_procHxBck_timelist(start, end)

finalHx = a.get_finalHx_timelist(start, end)
finalHxBck = a.get_finalHxBck_timelist(start, end)

bckPrior = a.get_bckPrior_timelist(start, end)
bckProc = a.get_bckProc_timelist(start, end)

freq = "daily"
priorEmiss = a.get_priorEmissPoints_timelist(start, end, freq)
priorSigma = a.get_priorSigmaPoints_timelist(start, end, freq)

posteriorEmiss = a.get_posteriorEmissPoints_timelist(start, end, freq)
posteriorSigma = a.get_posteriorSigmaPoints_timelist(start, end, freq)

ind = -1
"""
obs.data[ind].plot(label = "OBS")
initHxBck.data[ind].plot(label = "init")
priorHxBck.data[ind].plot(label = "prior")
procHxBck.data[ind].plot(label = "proc")
finalHxBck.data[ind].plot(label = "final")
plt.legend()
plt.show()

bckPrior.data[ind].plot(label = "Prior")
bckProc.data[ind].plot(label = "Proc")
plt.legend()
plt.show()
"""

priorEmiss.data[ind].plot(label = "Prior")
posteriorEmiss.data[ind].plot(label = "Posterior")
plt.legend()
plt.show()

priorSigma.data[ind].plot(label = "Prior")
posteriorSigma.data[ind].plot(label = "Posterior")
plt.legend()
plt.show()
