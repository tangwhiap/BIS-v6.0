#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from .indicator import Indicator
from ...utils.sparse_matrix import create_diag

import numpy as np

class SectorsIndicator(Indicator):
    def __init__(self, *args, **kwargs):
        """
            This function is used for testing.
            It couldn't be called by the final system!
        """
        print("Warning! The testing function has been called.")
        self.nType = 3
        self.typeNames = ["A", "B", "C"]
        #self.sector_to_type = {"A": "SiteObs", "B": "SiteObs", "C": "SatObs"}
        self.dim1Dict = {"A": 1230, "B": 1230, "C": 1000}
        self.dim2Dict = {"A": 1, "B": 1, "C": 1}
        self.data = {}
        for type_ in self.typeNames:
            dim1 = self.dim1Dict[type_]
            dim2 = self.dim2Dict[type_]
            self.data[type_]  = np.random.rand(dim1, dim2)

    def to_operator(self):
        from .indicator_operator import IndicatorOperator
        return IndicatorOperator(indicator = self) 

            
        
class Sigma_Indicator(SectorsIndicator):

    """
        Uncertainty (Sigma matrix) indicator.
        Overload the dim_val method.
        Add diag_out method used for 
    """

    def __init__(self, *args, **kwargs):
        SectorsIndicator.__init__(self, *args, **kwargs)
    
    def dim_val(self):
        SectorsIndicator.dim_val(self)
        for sector in self.dim2Dict:
            assert self.dim2Dict[sector] == 1, "Error! Dimension #2 must be 1. (current: " + str(self.dim2Dict[sector]) + ")"
    
    def diag_out(self):
        # Used to convert sigma vector to diagonal sparse matrix.
        diagDict = {}
        for type_ in self.typeNames:
            diagDict[type_] = create_diag(self.data[type_].reshape(self.dim1Dict[type_]))
        from .indicator_operator import IndicatorOperator
        return IndicatorOperator(data = diagDict, indicator_class = "sectors")


class X_Indicator(SectorsIndicator):

    """
        Emission (Sigma matrix) indicator.
        Overload the dim_val method.
    """
    def __init__(self, *args, **kwargs):
        SectorsIndicator.__init__(self, *args, **kwargs)
    
    def dim_val(self):
        SectorsIndicator.dim_val(self)
        for sector in self.dim2Dict:
            assert self.dim2Dict[sector] == 1, "Error! Dimension #2 must be 1. (current: " + str(self.dim2Dict[sector]) + ")"

class dX_Indicator(SectorsIndicator):

    """
        Emission (Sigma matrix) indicator.
        Overload the dim_val method.
    """
    def __init__(self, objIndOp, time):
        self.typeNames = objIndOp.typeNames
        self.nType = len(self.typeNames)
        self.data = objIndOp.data
        self.time = time
        self.dim1Dict = {}
        self.dim2Dict = {}
        for sector in self.typeNames:
            self.dim1Dict[sector], self.dim2Dict[sector] = self.data[sector].shape
        self.check()

    
    def dim_val(self):
        SectorsIndicator.dim_val(self)
        for sector in self.dim2Dict:
            assert self.dim2Dict[sector] == 1, "Error! Dimension #2 must be 1. (current: " + str(self.dim2Dict[sector]) + ")"

class dSigma_Indicator(SectorsIndicator):

    """
        Emission (Sigma matrix) indicator.
        Overload the dim_val method.
    """
    def __init__(self, resData, time, isSquare):
        self.isSquare = isSquare
        self.typeNames = list(resData.keys())
        self.nType = len(self.typeNames)
        self.data = resData
        self.time = time
        self.dim1Dict = {}
        self.dim2Dict = {}
        for sector in self.typeNames:
            self.dim1Dict[sector], self.dim2Dict[sector] = self.data[sector].shape
        self.check()

    
    def dim_val(self):
        SectorsIndicator.dim_val(self)
        for sector in self.dim2Dict:
            assert self.dim2Dict[sector] == 1, "Error! Dimension #2 must be 1. (current: " + str(self.dim2Dict[sector]) + ")"

