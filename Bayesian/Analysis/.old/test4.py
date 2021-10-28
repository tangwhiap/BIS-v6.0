#!/usr/bin/env python
import numpy as np
from OSSE_Analyst import OSSE_Analyst_zero_one as case
import datetime as dtm
import matplotlib.pyplot as plt
from pdb import set_trace
time = dtm.datetime(2021,1,4,0)
Start = dtm.datetime(2021,1,3,0)
End = dtm.datetime(2021,1,5,18)
a = case(CaseName = "Exp_OSSE=zero_one_Lt=48_Ls=100")
#obs = a.get_OBS(time = dtm.datetime(2021,1,4), receptors_list = ["IAPtower_h1", "IAPtower_h2", "XiangHe_h1"])
fig = plt.figure()


#tlist, RQtot_list = a.total_RQ_timelist(Start = Start, End = End)
#a.plot_tlist(tlist, RQtot_list,fig = fig)
emiss = a.get_reduced_total_unc_mean(Start = Start, End = End)
a.plot2D(emiss, fig = fig)
plt.show()
