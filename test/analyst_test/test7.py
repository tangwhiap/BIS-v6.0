#!/usr/bin/env python
from Bayesian.app import case
from Bayesian.Analysis.Analyst.bis_analyst import Analyst

import datetime as dtm
import matplotlib.pyplot as plt

from pdb import set_trace

a = Analyst(case)

start = dtm.datetime(2019, 1, 1)

end = dtm.datetime(2019, 1, 6)

a.build_region_maskout(regionName = "Beijing", shpName = ["Beijing"], LON = a.sampleLON, LAT = a.sampleLAT)
a.build_region_maskout(regionName = "Tianjin", shpName = ["Tianjin"], LON = a.sampleLON, LAT = a.sampleLAT)
"""
a.regionsMaskout["Beijing"].plot()
plt.show()
a.regionsMaskout["BT"].plot()
plt.show()
"""
timeScale = "hourly"
ds = a.get_posteriorEmissRegion_timelist(start = start, end = end, regions = ["Beijing", "Tianjin"], timeScale = timeScale, areaWeighted = False)
dss = a.get_posteriorSigmaRegion_timelist(start = start, end = end, regions = ["Beijing", "Tianjin"], timeScale = timeScale, areaWeighted = False)

time = ds.time.values
emiss_BJ = ds.sel(site = "Beijing").data.values
sigma_BJ = dss.sel(site = "Beijing").data.values
emiss_TJ = ds.sel(site = "Tianjin").data.values
sigma_TJ = dss.sel(site = "Tianjin").data.values

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(time, emiss_BJ, label = "Beijing", color = "red")
ax.plot(time, emiss_TJ, label = "Tianjin", color = "black")
ax.fill_between(time, emiss_BJ - sigma_BJ, emiss_BJ + sigma_BJ, color = "red", alpha = 0.3)
ax.fill_between(time, emiss_TJ - sigma_TJ, emiss_TJ + sigma_TJ, color = "grey", alpha = 0.3)

fig.autofmt_xdate()

ax.legend()
plt.show()
