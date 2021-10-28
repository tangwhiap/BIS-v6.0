#!/usr/bin/env python
import numpy as np
from OSSE_Analyst import OSSE_Analyst_zero_one as case
import datetime as dtm
import matplotlib.pyplot as plt

a = case(CaseName = "Exp_OSSE=zero_one_Lt=48_Ls=100")
time = a.pos_emiss_times[26]
#RecepTime = a.receptor_times[0]
print(time)
#print(RecepTime)
#foot = a.get_footprint(Receptors_list = ["IAPtower_h1"], time_list = time, RecepTime_list = [RecepTime])
foot = a.get_footprint(Receptors_list = ["IAPtower_h1", "LuanCheng_h1", "ShangDianZi_h1"], time_list = [time])
fig = plt.figure()
ax, cs = a.plot2D(foot, fig = fig)
fig.colorbar(cs)
plt.show()


