#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from .indicators import H_Indicator
from ..indicators.indicator_operator import IndicatorOperator_Zero
from pdb import set_trace

def __Hxi(H, x):
    HH = H.to_operator()
    xx = x.to_operator()
    return H.to_operator() * x.to_operator()

def Hx(timeList, Xtype, recepTime, class_XIndicator, objExp):
    sumHx = IndicatorOperator_Zero()
    for time in timeList:
        H = H_Indicator(time = time, recepTime = recepTime, typeToSector = objExp.typeToSector_H, 
                fileDir = objExp.myConfig["HDir"], filePrefix = objExp.myConfig["H_Prefix"], fun_get_filename = objExp.get_HFile_name)
        x = class_XIndicator(time = time, Xtype = Xtype, objExp = objExp)
        Hxi = __Hxi(H, x)
        sumHx = sumHx + Hxi
    return sumHx


def Hx2obs(method = "sum", *args, **kwargs):
    loc = locals()
    #print("Hello")
    #obs_Indopt = Hx2obs_sum(*args, **kwargs)
    #print("Done")
    exec("obs_Indopt = Hx2obs_" + method + "(*args, **kwargs)")
    obs_Indopt = loc["obs_Indopt"]
    return obs_Indopt

def Hx2obs_sum(sumHx):
    """
        obs = H1*x1 + H2*x2 + ... + Hn * xn
    """
    return sumHx.get_sum()



def compute_Hx(recepTimeList, backtime_j2i, **Hx_kwargs):
    if isinstance(recepTimeList, list):
        Hx_list = []
        for recepTime in timeList:
            timeList = backtime_j2i(recepTime)
            Hx_list.append(Hx(timeList, recepTime = recepTime, **Hx_kwargs))
        return Hx_list
    else:
        recepTime = recepTimeList
        timeList = backtime_j2i(recepTime)
        return Hx(timeList, recepTime = recepTime, **Hx_kwargs)

