#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from .indicator import Indicator

import numpy as np

class E_Indicator(Indicator):
    """
        Spatial correlation coef. (E matrix) Indicator.
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
        self.dim1Dict = {"A": 12300, "B": 12300, "C": 1000}
        self.dim2Dict = {"A": 12300, "B": 12300, "C": 1000}
        self.data = {}
        for type_ in self.typeNames:
            dim1 = self.dim1Dict[type_]
            dim2 = self.dim2Dict[type_]
            self.data[type_]  = np.random.rand(dim1, dim2)


