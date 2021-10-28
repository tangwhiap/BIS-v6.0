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
ds_daily = a.get_priorEmissRegion_timelist(start = start, end = end, regions = ["Beijing", "Tianjin"], timeScale = "daily")
ds_hourly = a.get_priorEmissRegion_timelist(start = start, end = end, regions = ["Beijing", "Tianjin"], timeScale = "hourly")

print(ds_hourly.sel(site = "Tianjin").data.values[:-1].sum())
print(ds_daily.sel(site = "Tianjin").data.values.sum())

ds_hourly.sel(site = "Beijing").data.plot(label = "Beijing")
ds_hourly.sel(site = "Tianjin").data.plot(label = "Tianjin")
plt.legend()
plt.show()
