#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...

from ...benchplots_base.main_class import BIS_BenchPlots
from ...benchplots_base import benchplots_configure  as benchConfig

class BenchPlots_REAL(BIS_BenchPlots):

    def __init__(self, *args, **kwargs):

        BIS_BenchPlots.__init__(self, *args, **kwargs)
        self.myConfig = {**self.myConfig, **benchConfig.expREAL_kwargs}