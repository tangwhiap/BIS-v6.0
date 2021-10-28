#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from .indicator import Indicator

import numpy as np
import pandas as pd
import datetime as dtm
from pdb import set_trace

class TypesIndicator(Indicator):
    def __init__(self, *args, **kwargs):
        """
            This function is used for testing.
            It couldn't be called by the final system!
        """
        print("Warning! The testing function has been called.")
        self.nType = 2
        self.typeNames = ["SiteObs", "SatObs"]
        self.typeToSector = {"SiteObs": ["A", "B"], "SatObs": ["C"]}
        self.get_sector_to_type()
        self.dim1Dict = {"SiteObs": 3, "SatObs": 100}
        self.dim2Dict = {"SiteObs": 1230, "SatObs": 1000}
        self.data = {}
        for type_ in self.typeNames:
            dim1 = self.dim1Dict[type_]
            dim2 = self.dim2Dict[type_]
            self.data[type_]  = np.random.rand(dim1, dim2)
        self.get_sector_to_type()

    #def dim_val(self, *args, **kwargs):
    #    Indicator.dim_val(self, *args, **kwargs)
    #    sectorHaved = []
    #    for type_ in self.typeToSector:
    #        for sector in self.typeToSector[type_]:
    #            assert sector not in sectorHaved
    #            sectorHaved.append(sector)

    def get_sector_to_type(self):
        self.sectorToType = {}
        for type_ in self.typeToSector:
            for sector in self.typeToSector[type_]:
                try:
                    assert sector not in self.sectorToType
                except:
                    set_trace()
                self.sectorToType[sector] = type_
            
    def to_operator(self):
        from .indicator_operator import IndicatorOperator
        return IndicatorOperator(indicator = self)

class H_Indicator(TypesIndicator):
    """
        Observation opterator (H matrix) Indicator.
    """
    def __init__(self, *args, **kwargs):
        TypesIndicator.__init__(self, *args, **kwargs)
       
    def trans(self):
        for type_ in self.typeNames:
            self.data[type_] = self.data[type_].T
        self.dim1Dict, self.dim2Dict = self.dim2Dict, self.dim1Dict
        return self




class D_Indicator(TypesIndicator):
    """
        Spatial correlation coef. (D matrix) indicator.
    """
    def __init__(self, *args, **kwargs):
        """
            This function is used for testing.
            It couldn't be called by the final system!
        """
        print("Waring! The testing function has been called.")
        self.nType = 3
        self.typeNames = ["AsD", "BsD", "CsD"]
        self.typeToSector = {"AsD": "A", "BsD": "B", "CsD": "C"}
        self.dim1Dict = {"AsD": 48, "BsD": 48, "CsD": 48}
        self.dim2Dict = {"AsD": 48, "BsD": 48, "CsD": 48}
        self.dim1TDict = {"AsD": [], "BsD": [], "CsD": []}
        self.dim2TDict = {"AsD": [], "BsD": [], "CsD": []}
        tlist = pd.date_range("2019-01-01", "2019-01-03", freq = "H")[:-1]
        for itime in tlist:
            for type_ in self.typeNames:
                self.dim1TDict[type_].append(itime.strftime("%Y%m%d%H%M%S"))
                self.dim2TDict[type_].append(itime.strftime("%Y%m%d%H%M%S"))
        self.data = {}
        for type_ in self.typeNames:
            dim1 = self.dim1Dict[type_]
            dim2 = self.dim2Dict[type_]
            self.data[type_]  = np.random.rand(dim1, dim2)
        self.get_sector_to_type()

    def __getitem__(self, key):

        assert isinstance(key, tuple)
        if len(key) == 2:
            time1, time2 = key
            typeList = self.typeNames
        elif len(key) == 3:
            time1, time2, typeList = key
        else:
            assert False, ""

        if isinstance(time1, dtm.datetime) or isinstance(time1, pd.Timestamp):
            time1 = time1.strftime("%Y%m%d%H%M%S")

        if isinstance(time2, dtm.datetime) or isinstance(time2, pd.Timestamp):
            time2 = time2.strftime("%Y%m%d%H%M%S")

        
        outDict = {}
        for type_ in typeList:
            assert time1 in self.dim1TDict[type_], "time1: " + time1 + " not found in D matrix"
            assert time2 in self.dim2TDict[type_], "time2: " + time2 + " not found in D matrix"
            dim1Index = self.dim1TDict[type_].index(time1)
            dim2Index = self.dim2TDict[type_].index(time2)
            outDict[type_] = self.data[type_][dim1Index, dim2Index]

        from .indicator_operator import IndicatorOperator
        #return outDict
        return IndicatorOperator(data = outDict, indicator_class = "types", sectorToType = self.sectorToType, typeToSector = self.typeToSector)


    #def __mul__(self, other):

class E_Indicator(TypesIndicator):
    """
        Spatial correlation coef. (E matrix) Indicator.
    """
    def __init__(self, *args, **kwargs):
        """
            This function is used for testing.
            It couldn't be called by the final system!
        """
        print("Waring! The testing function has been called.")
        self.nType = 3
        self.typeNames = ["AsE", "BsE", "CsE"]
        self.typeToSector = {"AsE": "A", "BsE": "B", "CsE": "C"}
        self.dim1Dict = {"AsE": 1230, "BsE": 1230, "CsE": 1000}
        self.dim2Dict = {"AsE": 1230, "BsE": 1230, "CsE": 1000}
        self.data = {}
        for type_ in self.typeNames:
            dim1 = self.dim1Dict[type_]
            dim2 = self.dim2Dict[type_]
            self.data[type_]  = np.random.rand(dim1, dim2)
        self.get_sector_to_type()


