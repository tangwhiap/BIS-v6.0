#!/usr/bin/env python

from Bayesian.core.optimize_core.compute_HQ import compute_HQj

import pandas as pd

ti_range = pd.date_range("2019-01-01", "2019-01-03", freq = "H")[:-1]
tj = ti_range[2]
recepTime = ti_range[-1]
HQj = compute_HQj(ti_range = ti_range, tj = tj, recepTime = recepTime)
print(HQj)
print("A")
print(HQj.get_matrix("A"))
print("B")
print(HQj.get_matrix("B"))
print("C")
print(HQj.get_matrix("C"))

