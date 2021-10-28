#!/usr/bin/env python


import matplotlib.pyplot as plt
#from BIS_Analyst import BIS_Case_Analyst
from OSSE_Analyst import OSSE_Analyst_zero_one

#a=BIS_Case_Analyst("/home/tangwh/modeldata/BIS_cases/Exp_OSSE=zero_one_Lt=48_Ls=500", "2021-01-03_00:00:00", "2021-01-06_00:00:00", CaseName = "Exp_OSSE=zero_one_Lt=48_Ls=500")
a=OSSE_Analyst_zero_one("/home/tangwh/modeldata/BIS_cases/Exp_OSSE=zero_one_Lt=48_Ls=500", "2021-01-03_00:00:00", "2021-01-06_00:00:00", CaseName = "Exp_OSSE=zero_one_Lt=48_Ls=500")

emiss = a.posterior_mean()
t, y = a.total_emiss_timelist(emiss_type = "posterior")
t, tr = a.total_emiss_timelist(get_emiss_fun = "get_prior_emiss")
fig = plt.figure()
ax1, cs1 = a.plot2D(emiss, fig = fig, subplot = (2,2,1), vmin = 0, vmax = 0.5)
cbar1 = fig.colorbar(cs1)
ax2, cs2 = a.plot2D(emiss, fig = fig, subplot = (2,2,2), vmin = 0.5, vmax = 1)
cbar2 = fig.colorbar(cs2)
ax3, line1 = a.plot_tlist(t, y, fig = fig, subplot = (2,1,2), label = "posterior")
_, line2 = a.plot_tlist(t, tr, ax = ax3, label = "prior")


#cbar1 = fig.colorbar(cs1)
#cbar2 = fig.colorbar(cs2)
ax3.legend()
plt.show()
