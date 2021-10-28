#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from .indicator import Indicator

import numpy as np

class H_indicator(Indicator):
    """
        Observation opterator (H matrix) Indicator.
    """
    def __init__(self, *args, **kwargs):
        """
            This function is used for testing.
            It couldn't be called by the final system!
        """
        print("Waring! The testing function has been called.")
        self.Ntype = 2
        self.typeNames = ["SiteObs", "SatObs"]
        self.dim1Dict = {"SiteObs": 3, "SatObs": 100}
        self.dim2Dict = {"SiteObs": 12300, "SatObs": 1000}
        self.data = {}
        for type_ in self.typeNames:
            dim1 = self.dim1Dict[type_]
            dim2 = self.dim2Dict[type_]
            self.data[type_]  = np.random.rand(dim1, dim2)


