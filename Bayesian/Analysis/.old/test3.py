#!/usr/bin/env python
import numpy as np
from OSSE_Analyst import OSSE_Analyst_zero_one as case
import datetime as dtm
import matplotlib.pyplot as plt
time = dtm.datetime(2021,1,4,0)
a = case(CaseName = "Exp_OSSE=zero_one_Lt=48_Ls=100")
#obs = a.get_OBS(time = dtm.datetime(2021,1,4), receptors_list = ["IAPtower_h1", "IAPtower_h2", "XiangHe_h1"])
obs = a.get_OBS(time = time)
print(obs)
fig = plt.figure()
ax, cs = a.draw_OBS_scatter(time = time, fig = fig, vmin = 0, vmax = 2)
cs.set_clim(0,3)
plt.show()
