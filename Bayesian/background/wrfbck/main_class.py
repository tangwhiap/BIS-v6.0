#!/usr/bin/env python
# Authors:
#   Wenhan TANG - 06/2021
#   Wenhan TANG - 07/2021
#   ...

import numpy as np
import netCDF4 as nc
import xarray as xr
from ..base.main_class import BCK_base
from .getWrfco2BCK import getBCK, make_LocDic

class BCK_wrfbck(BCK_base):
    def __init__(self, Receptors, Wrfco2Dir, Wrfco2Prefix, Domid, errors, offset = 0):

        #BCK_base.__init__(self, *args, **kwargs)
        self.BCK_name = "wrfbck"
        #self.myConfig = self.bckConfig[self.BCK_name]
        self.Wrfco2Dir = Wrfco2Dir
        self.Wrfco2Prefix = Wrfco2Prefix 
        self.Domid = Domid
        self.errors = errors
        self.offset = offset

        #self.Ls_bck = self.argv["Ls_bck"]
        self.LocDic = make_LocDic(Receptors, self.Wrfco2Dir, self.Wrfco2Prefix, self.Domid)
        self.nStation = len(Receptors)
        #self.compute_Ebck()
        
    def compute_BCK(self, Time):

        BCK_dict = getBCK(Time, self.LocDic, self.Wrfco2Dir, self.Wrfco2Prefix, self.Domid)
        BCK_list = [BCK_dict[station] for station in BCK_dict]
        BCK_list = np.array(BCK_list) + self.offset

        return BCK_list

    def compute_sigma_BCK(self, Time):
        errorConst = self.errors

        #Sigma_bck = np.diag(np.ones(Ns) * Qbck_const)
        #Qbck = np.matmul( np.matmul( Sigma_bck, self.Ebck ), Sigma_bck )

        errorList = [errorConst] * self.nStation
        errorList = np.array(errorList)

        return errorList


