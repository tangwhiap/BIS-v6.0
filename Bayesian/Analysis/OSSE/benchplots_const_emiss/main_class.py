#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 09/2021
#   ...

from ..benchplots_osse.main_class import BenchPlots_OSSE
from ...benchplots_base import benchplots_configure  as benchConfig

class BenchPlots_const_emiss(BenchPlots_OSSE):

    def __init__(self, *args, **kwargs):

        BenchPlots_OSSE.__init__(self, *args, **kwargs)
        self.myConfig = {**self.myConfig, **benchConfig.expOSSE_const_emiss_kwargs}
