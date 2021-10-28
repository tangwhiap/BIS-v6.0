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
freq = "daily"

def fun(freq):
    priorEmiss = a.get_priorEmissPoints_timelist(start, end, freq)
    priorSigma = a.get_priorSigmaPoints_timelist(start, end, freq)

    posteriorEmiss = a.get_posteriorEmissPoints_timelist(start, end, freq)
    posteriorSigma = a.get_posteriorSigmaPoints_timelist(start, end, freq)

    ind = -1

    priorEmiss.data[ind].plot(label = "Prior_" + freq)
    posteriorEmiss.data[ind].plot(label = "Posterior_" + freq)

for freq in ["hourly", "daily"]:
    print(freq)
    fun(freq)

plt.legend()
plt.show()
