#!/usr/bin/env python

from Bayesian.app import case
from Bayesian.Analysis.Analyst.bis_analyst import Analyst

import datetime as dtm
import matplotlib.pyplot as plt

from pdb import set_trace

a = Analyst(case)
LON, LAT = a.sampleLON, a.sampleLAT

g = a.obsPoints["IAPtower_h1"]
g.anchor(LON, LAT)
print(g.staLon)
print(g.belinear(a.sampleLON))
print(g.staLat)
print(g.belinear(a.sampleLAT))
