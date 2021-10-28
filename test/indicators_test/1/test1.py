#!/usr/bin/env python

from Bayesian.core.indicators.types_indicator import *
from Bayesian.core.indicators.sectors_indicator import *
from Bayesian.core.indicators.indicator_operator import *

import datetime as dtm

H = H_Indicator()
E = E_Indicator()
D = D_Indicator()
Sigma = Sigma_Indicator()

time1 = dtm.datetime(2019, 1, 2)
time2 = dtm.datetime(2019, 1, 2, 4)
print(D[time1, time2])
d = D[time1, time2]

H_cac = IndicatorOperator(indicator = H)
E_cac = IndicatorOperator(indicator = E)
HE_cac = H_cac * Sigma.diag_out()
