#!/usr/bin/env python

"""
CaseDict = {
               "casename": [CaseDir, Start, End, dt, window, n_tback], 
               "Ls200": ["Exp_OSSE=zero_one_Lt=48_Ls=200", "2021-01-03_00:00:00", "2021-01-06_00:00:00", 1, 6, 24]
           }
"""

Start = "2021-01-03_00:00:00"
End = "2021-01-05_18:00:00"
dt = 1
window = 6
n_tback = 24
CasesDir = "/home/tangwh/modeldata/BIS_cases"
Ls_list = [5, 10, 30, 50, 100, 200, 500, 1000]
CaseDict = {}
for Ls in Ls_list:
    CaseName = "Ls" + str(Ls)
    CaseDict[CaseName] = {
        #"CaseDir": CasesDir + "/Exp_OSSE=zero_one_Lt=48_Ls=" + str(Ls),
        "CaseDir": CasesDir + "/Exp_OSSE=const_meic_Lt=48_Ls=" + str(Ls),
        "Start": Start, 
        "End": End, 
        "CaseName": CaseName,
        "dt": dt, 
        "window": window, 
        "n_tback": n_tback
    }


