#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from ..indicators.temp_indicator import HQ_Indicator, HQHt_Indicator
from ...utils.netcdf_io import get_ncvar, nc_write

import numpy as np

def compute_HQHt_R(objCore, objExp, recepTime):

    HQHt_Ind = HQHt_Indicator(sectors = objExp.sectors, get_HQHtFile_name = objCore.get_HQHtFile_name, recepTime = recepTime)

    HQHt = HQHt_Ind.data_sum()

    R = get_ncvar(objExp.get_RFile_offline_name(recepTime = recepTime), "data")

    nc_write(objCore.get_HQHtRFile_name(recepTime = recepTime), HQHt + R)

