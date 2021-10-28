#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import datetime as dtm
import xarray as xr
from pdb import set_trace


str_StartDT = "2021-01-03_00:00:00"
str_EndDT = "2021-01-06_00:00:00"

dtHrs = 1
ShpDir = "/home/tangwh/datasets/china_shp/cnhimap.shp"
#InDir = "/home/tangwh/modeling/BIS_cases/test_lt=6_ls=10/output/processing"
#InDir = "/home/tangwh/modeling/BIS_cases/test1/output/processing"
#InDir = "/home/tangwh/modeling/BIS_cases/test_lt=6_ls=30/output/processing"
InDir = "/home/tangwh/modeling/BIS_cases/osse2_lt=6_ls=30/output/processing"
#InDir = "/home/tangwh/modeling/BIS_cases/Exp_OSSE=zero_one_Lt=48_Ls=1000/output/processing"
OutDir = "./output"
InPrefix = "BISOUT"
OutPrefix = "POS"

Start = dtm.datetime.strptime(str_StartDT, "%Y-%m-%d_%H:%M:%S")
End = dtm.datetime.strptime(str_EndDT, "%Y-%m-%d_%H:%M:%S")
dt = dtm.timedelta(hours = dtHrs)


def areaS(lon, lat):
    R = 6371000
    dlon = (lon[1:] - lon[:-1]).mean()
    dlat = (lat[1:] - lat[:-1]).mean()
    LON, LAT = np.meshgrid(lon, lat)
    S = (R**2) * np.deg2rad(dlon) * (np.sin(np.deg2rad(LAT + dlat/2 )) - np.sin(np.deg2rad(LAT - dlat/2 )))
    assert S.shape[0] == len(lat)
    assert S.shape[1] == len(lon)
    return S

time_tlist = [] 
emiss_total_tlist = []
truth_total_tlist = []
Current = Start
while( Current < End ):
    print(Current.strftime("%Y-%m-%d_%H:%M:%S"))
    ds = xr.open_dataset(InDir + "/" + InPrefix + "_" + Current.strftime("%Y-%m-%d_%H:%M:%S") + ".nc")
    emiss = ds["emiss"].values
    lon = ds["lon"].values
    lat = ds["lat"].values
    S = areaS(lon, lat)
    # umol/m^2/s --> tCO2
    emiss_total = (emiss * S).sum() * 3600 / 1000000 * 44 / 1000000
    truth_total = ((emiss * 0 + 1) * S).sum() * 3600 / 1000000 * 44 / 1000000
    time_tlist.append(np.datetime64(Current))
    emiss_total_tlist.append(emiss_total)
    truth_total_tlist.append(truth_total)
    #print(emiss_total)
    emiss_total_tlist
    Current += dt

time_tlist = np.array(time_tlist)
emiss_total_tlist = np.array(emiss_total_tlist)
set_trace()
fig = plt.figure()
ax = fig.add_subplot(2,1,1)
ax.plot(time_tlist, emiss_total_tlist)
ax.set_ylim([-10000, 200000])
ax.set_ylabel("tCO2/hour")
#ax.plot(time_tlist, truth_total_tlist)
ax = fig.add_subplot(2,1,2)
ax.set_ylim([-10000, 200000])
ax.plot(time_tlist, emiss_total_tlist, label = "posterior")
ax.plot(time_tlist, truth_total_tlist, label = "truth")
#ax.set_ylim([-10000, 200000])
ax.set_ylabel("tCO2/hour")
ax.legend()
plt.show()
