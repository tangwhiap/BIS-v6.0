#!/usr/bin/env python
from footprint_stilt import stilt

import datetime as dtm
import matplotlib.pyplot as plt

foot = stilt(receptors = ["IAPtower_h1", "IAPtower_h2", "XiangHe_h1"])
time = dtm.datetime(2020, 9, 15, 20)
recepTime = dtm.datetime(2020, 9, 15, 23)

footprint, LON, LAT = foot.get_footprint("IAPtower_h1", time, recepTime)

#plt.pcolormesh(LON, LAT, footprint, cmap = "jet")
#plt.show()

