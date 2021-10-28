#!/usr/bin/env python
from Bayesian.Analysis.benchplots_base.main_class import BIS_BenchPlots
from Bayesian.app import case
import datetime as dtm
start = dtm.datetime(2019, 1, 1, 0)
end = dtm.datetime(2019, 1, 6, 0)

objBench = BIS_BenchPlots(".", case)
#objBench.make_cartoon(start, end, 1)
#objBench.site_timeSeries(start, end)
objBench.pointRegion_timeSeries(start, end, "daily")
#objBench.draw_fieldPlots(start, end, "daily")
