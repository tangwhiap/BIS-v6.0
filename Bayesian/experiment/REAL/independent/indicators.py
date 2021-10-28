#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from  ..REAL_base import indicators as base_indicators
from ....utils.netcdf_io import get_ncdim, nc_addvar


class Sigma_Indicator(base_indicators.Sigma_Indicator):
    pass

class X_Indicator(base_indicators.X_Indicator):
    pass


class dX_Indicator(base_indicators.dX_Indicator):

    def modify_file(self, get_fileName):

        base_indicators.dX_Indicator.modify_file(self, get_fileName)

        for sector in self.typeNames:
            fileName = get_fileName(self.time, sector)
            nY, nX = get_ncdim(fileName, ["dimY", "dimX"])
            dArr = self.data[sector].reshape(nY, nX)
            nc_addvar(fileName, "emiss", dArr)

class dSigma_Indicator(base_indicators.dSigma_Indicator):

    def modify_file(self, get_fileName):

        dArr = base_indicators.dSigma_Indicator.modify_file(self, get_fileName)

        for sector in self.typeNames:
            fileName = get_fileName(self.time, sector)
            nY, nX = get_ncdim(fileName, ["dimY", "dimX"])
            #dArr = self.data[sector].reshape(nY, nX)
            dArr = dArr.reshape(nY, nX)
            nc_addvar(fileName, "emiss", dArr)
