#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from ...core.indicators import sectors_indicator as base_sectors_indicator
from ...core.indicators import types_indicator as base_types_indicator
from ...utils.netcdf_io import get_ncvar, nc_addvar, nc_write

import numpy as np
from pdb import set_trace

#import Sigma_Indicator, X_Indicator


class Sigma_Indicator(base_sectors_indicator.Sigma_Indicator):
    """
        Overload the constructor of Sigma_Indicator from indicators.sectors_indicator

        Getting the data from sigma files ("Prior" or "Proc") by objExp.get_sigma<type>File_name function,
        and complete the attribute (self.data, self.dim1Dict, self.dim2Dict, self.typeNames, self.nTypes) of the object.


    """
    Stype_list = ["Prior", "Proc"]
    def __init__(self, time, Stype, objExp):

        assert Stype in self.Stype_list
        self.typeNames = objExp.myConfig["sectors"]
        self.nTypes = len(self.typeNames)
        loc = locals()
        exec("get_fileName = objExp.get_sigma" + Stype + "File_name")
        get_fileName = loc["get_fileName"]
        self.data = {}
        self.dim1Dict = {}
        self.dim2Dict = {}
        for sector in self.typeNames:
            fileName = get_fileName(time, sector)
            data = get_ncvar(fileName, ["data"])
            self.data[sector] = data
            #set_trace()
            self.dim1Dict[sector], self.dim2Dict[sector] = data.shape 
        self.check()
    #def __init__(self, *args, **kwargs):



class X_Indicator(base_sectors_indicator.X_Indicator):
    """
        Overload the constructor of X_Indicator from indicators.sectors_indicator

        Getting the data from emission files ("prior" or "proc") by objExp.get_<type>File_name function,
        and complete the attribute (self.data, self.dim1Dict, self.dim2Dict, self.typeNames, self.nTypes) of the object.


    """
    #def __init__(self, time, get_fileName, exp_myConfig):
    Xtype_list = ["Prior", "Proc", "Post"]
    def __init__(self, time, Xtype, objExp):

        #print(Xtype, self.Xtype_list)
        assert Xtype in self.Xtype_list
        self.typeNames = objExp.myConfig["sectors"]
        self.nTypes = len(self.typeNames)
        loc = locals()
        exec("get_fileName = objExp.get_" + Xtype.lower() + "File_name")
        get_fileName = loc["get_fileName"]
        self.data = {}
        self.dim1Dict = {}
        self.dim2Dict = {}
        for sector in self.typeNames:
            fileName = get_fileName(time, sector)
            data = get_ncvar(fileName, ["data"])
            self.data[sector] = data
            #set_trace()
            self.dim1Dict[sector], self.dim2Dict[sector] = data.shape 
        self.check()


class dX_Indicator(base_sectors_indicator.dX_Indicator):

    """
        Add the function of modifying emission files.
    """
    def __init__(self, *args, **kwargs):
        base_sectors_indicator.dX_Indicator.__init__(self, *args, **kwargs)

    def modify_file(self, get_fileName):
        for sector in self.typeNames:
            fileName = get_fileName(self.time, sector)
            dArr = self.data[sector]
            nc_addvar(fileName, "data", dArr)

class dSigma_Indicator(base_sectors_indicator.dSigma_Indicator):

    """
        Add the function of modifying sigma files.
    """
    def __init__(self, *args, **kwargs):
        base_sectors_indicator.dSigma_Indicator.__init__(self, *args, **kwargs)

    def to_file(self, get_fileName, recepTime):
        for sector in self.typeNames:
            fileName = get_fileName(self.time, recepTime, sector)
            dSigma2 = self.data[sector]
            nc_write(fileName, dSigma2)

    def modify_file(self, get_fileName):
        if self.isSquare:
            for sector in self.typeNames:
                fileName = get_fileName(self.time, sector)
                sigma = get_ncvar(fileName, "data")
                dSigma2 = self.data[sector]
                sigma2 = sigma ** 2
                dSigma2 = np.where(dSigma2 > sigma2, 0, dSigma2)
                newSigma2 = sigma2 - dSigma2
                newSigma = np.sqrt(newSigma2)
                dArr = newSigma - sigma
                nc_addvar(fileName, "data", dArr)
                
        else:
            for sector in self.typeNames:
                fileName = get_fileName(self.time, sector)
                dArr = -self.data[sector]
                nc_addvar(fileName, "data", dArr)

        return dArr
#class H_Indicator(base_types_indicator.H_Indicator):
#    def __init__(self, *args, **kwargs):
