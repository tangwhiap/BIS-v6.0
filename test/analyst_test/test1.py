#!/usr/bin/env python

from Bayesian.app import case
from Bayesian.Analysis.Analyst.bis_analyst import Analyst

import datetime as dtm
import matplotlib.pyplot as plt

a = Analyst(case)


time = dtm.datetime(2019, 1, 3)

recepTime = dtm.datetime(2019, 1, 10)

#unc, LON, LAT = a.get_reduced_uncertainty(time, recepTime)
unc, LON, LAT = a.get_reduced_uncertainty(time, recepTime)

plt.contourf(LON, LAT, unc, cmap = "jet")

plt.colorbar()

plt.show()

