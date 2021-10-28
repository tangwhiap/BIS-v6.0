#!/usr/bin/env python

# Indicators for optimization core.

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from  ..indicators import types_indicator as  base_types_indicators
from ...utils.netcdf_io import get_ncvar, get_nctime
from ...utils.sparse_matrix import SparseMatrix, getNCSPR_spr

import numpy as np
from pdb import set_trace

def types_indicator_constructor(self, *args, **kwargs):
    """
        Constructor for types_indicator defined in package \"optimize_core\"
        * Argument:
            time:       Time of footprint (or other thing which can represent the observation operator.
            recepTime:  Time of observation.
            typeToSector: A dictionary with the key of each type name whose value is a list of sector corresponded.
                          ( default: {"": [""]} )
            fileDir:    Directory of data file.
            fun_get_filename : A handle of a function which can return the name of file. Used to define the format of file name.
    """
    #assert "time" in kwargs, "Error! key word \"time\" not found."
    #assert "recepTime" in kwargs, "Error! key word \"recepTime\" not found."
    if "time" in kwargs:
        time = kwargs["time"]
    else:
        time = None

    if "recepTime" in kwargs:
        recepTime = kwargs["recepTime"]
    else:
        recepTime = None
    
    if "typeToSector" in kwargs:
        self.typeToSector = kwargs["typeToSector"]
    else:
        self.typeToSector = {"": [""]}

    self.typeNames = list(self.typeToSector.keys())
    self.nType = len(self.typeNames)

    self.get_sector_to_type()

    if "fileDir" in kwargs:
        self.fileDir = kwargs["fileDir"]

    if "filePrefix" in kwargs:
        self.filePrefix = kwargs["filePrefix"]

    if "fun_get_filename" in kwargs:
        self.fun_get_filename = kwargs["fun_get_filename"]
    else:
        self.fun_get_filename = None

    self.data = {}
    self.dim1Dict = {}
    self.dim2Dict = {}

    for type_ in self.typeNames:

        assert self.matrixType in ["sparse", "sparseFile", "numpy"]

        if self.matrixType == "sparse":
            self.data[type_] = getNCSPR_spr(self.get_filename(time = time, recepTime = recepTime, typeName = type_))

        if self.matrixType == "sparseFile":
            self.data[type_] = SparseMatrix(self.get_filename(time = time, recepTime = recepTime, typeName = type_)) 

            #self.data[type_].close()
        if self.matrixType == "numpy":
            self.data[type_] = get_ncvar(self.get_filename(time = time, recepTime = recepTime, typeName = type_), "data")

        dim1, dim2 = self.data[type_].shape
        self.dim1Dict[type_] = dim1
        self.dim2Dict[type_] = dim2

    

    return self

class H_Indicator(base_types_indicators.H_Indicator):

    fileDir = None
    filePrefix = "H"
    matrixType = "sparse"

    def __init__(self, *args, **kwargs):
        self = types_indicator_constructor(self, *args, **kwargs)

    def get_filename(self, time, recepTime, typeName, *args, **kwargs):
        if self.fun_get_filename is None:
            return self.fileDir + "/" + self.filePrefix + "_t" + time.strftime("%Y%m%d%H%M") + "_r" + recepTime.strftime("%Y%m%d%H%M") + "_" + typeName + ".nc"
        else:
            return self.fun_get_filename(time, recepTime, typeName, *args, **kwargs)

    def close(self):
        for sect in self.data:
            if isinstance(self.data[sect], SparseMatrix):
                self.data[sect].close()
                
class E_Indicator(base_types_indicators.E_Indicator):

    fileDir = None
    filePrefix = "E"
    matrixType = "sparseFile"

    def __init__(self, *args, **kwargs):
        self = types_indicator_constructor(self, *args, **kwargs)

    def get_filename(self, typeName, *args, **kwargs):
        if self.fun_get_filename is None:
            return self.fileDir + "/" + self.filePrefix + "_" + typeName + ".nc"
        else:
            return self.fun_get_filename(typeName, *args, **kwargs)

    def close(self):
        for sect in self.data:
            if isinstance(self.data[sect], SparseMatrix):
                self.data[sect].close()
                
class D_Indicator(base_types_indicators.D_Indicator):

    fileDir = None
    filePrefix = "D"
    matrixType = "numpy"

    def __init__(self, *args, **kwargs):
        self = types_indicator_constructor(self, *args, **kwargs)
        #time = kwargs["time"]
        #recepTime = kwargs["recepTime"]
        self.dim1TDict = {}
        self.dim2TDict = {}
        for type_ in self.typeNames:
            self.dim1TDict[type_] = [time.strftime("%Y%m%d%H%M%S") for time in get_nctime(self.get_filename(typeName = type_), ["dim1TList"])]
            self.dim2TDict[type_] = [time.strftime("%Y%m%d%H%M%S") for time in get_nctime(self.get_filename(typeName = type_), ["dim2TList"])]

    def get_filename(self, typeName, *args, **kwargs):
        if self.fun_get_filename is None:
            return self.fileDir + "/" + self.filePrefix + "_" + typeName + ".nc"
        else:
            return self.fun_get_filename(typeName, *args, **kwargs)

"""
class Sigma_Indicator(base_indicators.types_indicator.D_Indicator):

    fileDir = None
    filePrefix = "Sigma"
    isSparse = False

    def __init__(self, *args, **kwargs):
        self = types_indicator_constructor(self, *args, **kwargs)

    def get_filename(self, typeName):
        if self.fun_get_filename is None:
            return self.filePrefix + "_" + typeName + ".nc"
        else:
            return self.fun_get_filename(typeName)   
"""
