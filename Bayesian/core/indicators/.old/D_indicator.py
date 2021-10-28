#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from .indicator import Indicator
from ...utils.sparse_matrix import SparseMatrix

import numpy as np
import pandas as pd

class D_indicator(Indicator):
    """
        Spatial correlation coef. (D matrix) indicator.
    """
    def __init__(self, *args, **kwargs):
        """
            This function is used for testing.
            It couldn't be called by the final system!
        """
        print("Waring! The testing function has been called.")
        self.Ntype = 3
        self.typeNames = ["A", "B", "C"]
        self.sector_to_type = {"A": "SiteObs", "B": "SiteObs", "C": "SatObs"}
        self.dim1Dict = {"A": 48, "B": 48, "C": 48}
        self.dim2Dict = {"A": 48, "B": 48, "C": 48}
		self.dim1TDict = {"A": [], "B": [], "C": []}
		self.dim2TDict = {"A": [], "B": [], "C": []}
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

    #def __mul__(self, other):
		
