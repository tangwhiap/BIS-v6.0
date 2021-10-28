#!/usr/bin/env python

from Bayesian.background.wrfbck.main_class import BCK_wrfbck
from Bayesian.main.configure import BCKConfig

import datetime as dtm

Receptors = {
    "IAPtower_h1": ("116.3667", "39.9667", "80", 0.0),
    "IAPtower_h2": ("116.3667", "39.9667", "280", 0.0),
}

a = BCK_wrfbck(BCKConfig = BCKConfig, Receptors = Receptors)
time = dtm.datetime(2019,1,5,3)
print(a.compute_BCK(time))
print(a.compute_sigma_BCK(time))

