#!/usr/bin/env python

from Bayesian.emiss.NDRC import interface

import numpy as np
import xarray as xr
import datetime as dtm
import matplotlib.pyplot as plt

sectorName = "heat"

ds = xr.open_dataset("/home/tangwh/modeling/BIS_v6.0/cases/PIS_Ls3/emiss/processing/Proc_total_2019-01-28_17:00:00.nc")
LON = ds.LON.values
LAT = ds.LAT.values
time = dtm.datetime.strptime("2019-07-01_00:00:00", "%Y-%m-%d_%H:%M:%S")

emiss = interface(time, LON, LAT, sector = sectorName)

plt.pcolormesh(LON, LAT, np.log10(emiss), cmap = "coolwarm")
plt.colorbar()
plt.show()
