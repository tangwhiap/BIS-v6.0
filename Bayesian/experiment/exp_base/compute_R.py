#!/usr/bin/env python
# Authors:
#   Wenhan TANG - 07/2021
#   ...

import numpy as np
from pdb import set_trace
from ...utils.netcdf_io import nc_write

def diagError(recepTime, error_inReceptors, **kwargs):
    if isinstance(error_inReceptors, dict):
        error_inReceptors = np.array([error_inReceptors[receptor] for receptor in error_inReceptors])
    assert isinstance(error_inReceptors, np.ndarray) and len(error_inReceptors.shape) == 1
    return np.diag(error_inReceptors)



funNameDic = {
    "diagError": diagError,
}

def compute_R(time = None, recepTime = None, funName = None, get_RFile_name = None, **kwargs):
    assert funName in funNameDic
    assert get_RFile_name is not None
    compute_sigma_kwargs = {
        "time": time,
        "recepTime": recepTime,
    }
    compute_sigma_kwargs = {**compute_sigma_kwargs, **kwargs}
    func = funNameDic[funName]
    fileName = get_RFile_name(recepTime = recepTime)
    R = func(**compute_sigma_kwargs)
    nc_write(fileName, R)



