#!/usr/bin/env python

from Bayesian.core.indicators.indicator_operator import IndicatorOperator

import numpy as np

data = {}
sectorNames = ["A", "B"]
for sector in sectorNames:
    data[sector] = np.random.rand(10, 20)

obj = IndicatorOperator(data = data)
