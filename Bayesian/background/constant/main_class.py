#!/usr/bin/env python
# Authors:
#   Wenhan TANG - 06/2021
#   Wenhan TANG - 07/2021
#   ...

import numpy as np
import netCDF4 as nc
from ..base.main_class import BCK_base

class BCK_constant(BCK_base):

    def __init__(self, Receptors, const, errors):

        self.BCK_name = "constant"
        self.Receptors = Receptors
        self.nStation = len(self.Receptors)
        self.const = const
        self.errors = errors
        
    def compute_BCK(self, Time):

        BCK_list = np.array([self.const] * self.nStation)

        return BCK_list

    def compute_sigma_BCK(self, Time):

        errorConst = self.errors
        errorList = [errorConst] * self.nStation
        errorList = np.array(errorList)

        return errorList


