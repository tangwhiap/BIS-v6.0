#!/usr/bin/env python

#-> Base class of observation API <-#

# Authors:
#   Wenhan TANG - 07/2021
#   ...

import numpy as np

class OBS_base(object):

    makeType = "offline" # or "online"

    def __init__(self):
        """
            This function is used for testing.
        """
        print("Warning! test function (OBS_base constructor) is called.")
        self.obsName_dic = {"pr003CMAH1": 0, "pi368": 1, "pi602": 2}

    def get_obs(self, time):
        pass

    def get_error(self, time):
        pass

    def filter_min_max(self, arr, vmin = -np.inf, vmax = np.inf):
        return np.where( (arr >= vmin) & (arr <= vmax), arr, np.nan)
