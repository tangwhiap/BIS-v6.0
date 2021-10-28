#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...

from ..benchplots_real.main_class import BenchPlots_REAL
from ...benchplots_base import benchplots_configure  as benchConfig

class BenchPlots_independent(BenchPlots_REAL):

    def __init__(self, *args, **kwargs):

        BenchPlots_REAL.__init__(self, *args, **kwargs)
        self.myConfig = {**self.myConfig, **benchConfig.expREAL_independent_kwargs}
