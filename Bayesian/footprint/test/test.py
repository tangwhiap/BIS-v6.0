#!/usr/bin/env python
from footprint_stilt import interface

import numpy as np
import datetime as dtm
import matplotlib.pyplot as plt

receptorName = "IAPtower_h1"
#location = {
#            "lon": "116.3667",
#            "lat": "39.9667",
#            "agl": "80",
#            }

time = "202009152000"
recepTime = "202009152300"


recepTime = dtm.datetime.strptime(recepTime, "%Y%m%d%H%M")
time = dtm.datetime.strptime(time, "%Y%m%d%H%M")

footPrint, lonList, latList = interface(receptorName, location, time, recepTime)

plt.pcolormesh(lonList, latList, np.log10(footPrint))
plt.show()
