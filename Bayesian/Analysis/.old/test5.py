#!/usr/bin/env python
from BIS_Analyst import BIS_Case_Analyst as ana
import datetime as dtm
a=ana(CaseName = "real_test3")
v = a.emiss_point_value(dtm.datetime(2019,1,3,1), 116, 40, EmissType = "Posterior", get_variance = True)
print(v)
