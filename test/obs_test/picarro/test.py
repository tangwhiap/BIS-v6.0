#!/usr/bin/env python

from Bayesian.observation.picarro.main_class import OBS_picarro
from Bayesian.main.configure import obsConfig

import datetime as dtm

myConfig = obsConfig["picarro"]
a = OBS_picarro(myConfig["directory"], myConfig["sites"], myConfig["errors"])

a.get_data_offline("2019-01-01_00:00:00", "2019-01-03_00:00:00", 1)

current = dtm.datetime(2019,1,2,4)

print("orig")
print(a.get_obs_orig(current))

print("proc")
print(a.get_obs_proc(current))
