#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from .indicator import Indicator
from .sectors_indicator import SectorsIndicator
from .indicator_operator import IndicatorOperator
from ...utils.sparse_matrix import Sparse_to_nc, getNCSPR_spr
from ...utils.netcdf_io import get_ncvar, nc_write

from scipy.sparse import csr_matrix
from pdb import set_trace


class HQ_Indicator(SectorsIndicator):
    def __init__(self, *args, **kwargs):
        HQ_constructor(self, *args, **kwargs)

    def to_file(self, *args, **kwargs):
        HQ_toFile(self, *args, **kwargs)

    def trans(self):
        for sector in self.typeNames:
            self.data[sector] = self.data[sector].T
        self.dim1Dict, self.dim2Dict = self.dim2Dict, self.dim1Dict
        return self

    def sel_dim1(self, ind):
        resData = {}
        for sector in self.typeNames:
            resData[sector] = self.data[sector][ind, :]
        return HQ_Indicator(data = resData, jtime = self.jtime, recepTime = self.recepTime)

    def sel_dim2(self, ind):
        resData = {}
        for sector in self.typeNames:
            resData[sector] = self.data[sector][:, ind]
        return HQ_Indicator(data = resData, jtime = self.jtime, recepTime = self.recepTime)

    def __getitem__(self, key):
        resData = {}
        for sector in self.typeNames:
            resData[sector] = self.data[sector][key]
        return HQ_Indicator(data = resData, jtime = self.jtime, recepTime = self.recepTime)


class HdQ_Indicator(SectorsIndicator):

    def __init__(self, *args, **kwargs):
        HdQ_constructor(self, *args, **kwargs)

    def to_file(self, *args, **kwargs):
        HdQ_toFile(self, *args, **kwargs)

    def trans_to_sparse(self):
        for sector in self.typeNames:
            self.data[sector] = csr_matrix(self.data[sector])
        return self

class HQHt_Indicator(SectorsIndicator):
    def __init__(self, *args, **kwargs):
        HQHt_constructor(self, *args, **kwargs)

    def to_file(self, *args, **kwargs):
        HQHt_toFile(self, *args, **kwargs)
    
    def data_sum(self):
        return HQHt_sum(self)


class INVHQ_Indicator(SectorsIndicator):
    def __init__(self, *args, **kwargs):
        INVHQ_constructor(self, *args, **kwargs)

    def to_file(self, *args, **kwargs):
        INVHQ_toFile(self, *args, **kwargs)

class HiHQti_Indicator(SectorsIndicator):
    def __init__(self, *args, **kwargs):
        HiHQti_constructor(self, *args, **kwargs)

    def to_file(self, *args, **kwargs):
        HiHQti_toFile(self, *args, **kwargs)

#class INV_Indicator(Indicator):
#
#    def __init__(self, fileName):
#        self.INV = get_ncvar(fileName, "data")


def HQ_constructor(self, *args, **kwargs):

    assert "jtime" in kwargs
    assert "recepTime" in kwargs

    self.jtime = kwargs["jtime"]
    self.recepTime = kwargs["recepTime"]

    self.dim1Dict = {}
    self.dim2Dict = {}

    if "objIndOp" in kwargs:
        objIndOp = kwargs["objIndOp"]
        self.data = objIndOp.data
        self.typeNames = objIndOp.typeNames
        isFile = False

    elif "sectors" in kwargs and "get_HQFile_name" in kwargs and "jtime" in kwargs and "recepTime" in kwargs:
        get_HQFile_name = kwargs["get_HQFile_name"]
        self.typeNames = kwargs["sectors"]
        self.data = {}
        for sector in self.typeNames:
            self.data[sector] = getNCSPR_spr(get_HQFile_name(jtime = kwargs["jtime"], recepTime = kwargs["recepTime"], sectorName = sector))
        isFile = True

    else:
        assert "data" in kwargs# and "jtime" in kwargs, and "recepTime" in kwargs
        self.data = kwargs["data"]
        #self.jtime = kwargs["jtime"]
        #self.recepTime = kwargs["recepTime"]
        self.typeNames = list(self.data.keys())

    for sector in self.data:
        self.dim1Dict[sector], self.dim2Dict[sector] = self.data[sector].shape

def HdQ_constructor(self, *args, **kwargs):

    assert "jtime" in kwargs
    assert "recepTime" in kwargs

    self.jtime = kwargs["jtime"]
    self.recepTime = kwargs["recepTime"]

    self.dim1Dict = {}
    self.dim2Dict = {}

    if "objIndOp" in kwargs:
        objIndOp = kwargs["objIndOp"]
        self.data = objIndOp.data
        self.typeNames = objIndOp.typeNames
        isFile = False

    elif "sectors" in kwargs and "get_HdQFile_name" in kwargs and "jtime" in kwargs and "recepTime" in kwargs:
        get_HdQFile_name = kwargs["get_HdQFile_name"]
        self.typeNames = kwargs["sectors"]
        self.data = {}
        for sector in self.typeNames:
            self.data[sector] = get_ncvar(get_HdQFile_name(jtime = kwargs["jtime"], recepTime = kwargs["recepTime"], sectorName = sector), "data")
        isFile = True

    for sector in self.data:
        self.dim1Dict[sector], self.dim2Dict[sector] = self.data[sector].shape

def HQHt_constructor(self, *args, **kwargs):

    assert "recepTime" in kwargs

    self.recepTime = kwargs["recepTime"]

    if "objIndOp" in kwargs:
        objIndOp = kwargs["objIndOp"]
        self.data = objIndOp.data
        self.typeNames = objIndOp.typeNames
        for sector in self.data:
            self.data[sector] = self.data[sector].toarray()
        isFile = False

    elif "sectors" in kwargs and "get_HQHtFile_name" in kwargs and "recepTime" in kwargs:
        get_HQHtFile_name = kwargs["get_HQHtFile_name"]
        self.typeNames = kwargs["sectors"]
        self.data = {}
        for sector in self.typeNames:
            self.data[sector] = get_ncvar(get_HQHtFile_name(recepTime = kwargs["recepTime"], sectorName = sector), "data")
        isFile = True
    else:
        assert False

def INVHQ_constructor(self, *args, **kwargs):

    assert "jtime" in kwargs
    assert "recepTime" in kwargs

    self.jtime = kwargs["jtime"]
    self.recepTime = kwargs["recepTime"]

    self.dim1Dict = {}
    self.dim2Dict = {}

    if "objIndOp" in kwargs:
        objIndOp = kwargs["objIndOp"]
        self.data = objIndOp.data
        self.typeNames = objIndOp.typeNames
        isFile = False

    elif "sectors" in kwargs and "get_INVHQFile_name" in kwargs and "jtime" in kwargs and "recepTime" in kwargs:
        get_INVHQFile_name = kwargs["get_INVHQFile_name"]
        self.typeNames = kwargs["sectors"]
        self.data = {}
        for sector in self.typeNames:
            self.data[sector] = get_ncvar(get_INVHQFile_name(jtime = kwargs["jtime"], recepTime = kwargs["recepTime"], sectorName = sector), "data")
        isFile = True

    else:
        assert "data" in kwargs# and "jtime" in kwargs, and "recepTime" in kwargs
        self.data = kwargs["data"]
        self.typeNames = list(self.data.keys())

    for sector in self.data:
        self.dim1Dict[sector], self.dim2Dict[sector] = self.data[sector].shape

def HiHQti_constructor(self, *args, **kwargs):

    assert "recepTime" in kwargs
    assert "hisTime" in kwargs

    self.recepTime = kwargs["recepTime"]
    self.hisTime = kwargs["hisTime"]

    self.dim1Dict = {}
    self.dim2Dict = {}

    if "objIndOp" in kwargs:
        objIndOp = kwargs["objIndOp"]
        self.data = objIndOp.data
        self.typeNames = objIndOp.typeNames
        isFile = False

    elif "sectors" in kwargs and "get_HiHQtiFile_name" in kwargs and "recepTime" in kwargs:
        get_HiHQtiFile_name = kwargs["get_HiHQtiFile_name"]
        self.typeNames = kwargs["sectors"]
        self.data = {}
        for sector in self.typeNames:
            self.data[sector] = getNCSPR_spr(get_HiHQtiFile_name(recepTime = kwargs["recepTime"], hisTime=  kwargs["hisTime"], sectorName = sector))
        isFile = True

    else:
        assert "data" in kwargs# and "jtime" in kwargs, and "recepTime" in kwargs
        self.data = kwargs["data"]
        self.typeNames = list(self.data.keys())

    for sector in self.data:
        self.dim1Dict[sector], self.dim2Dict[sector] = self.data[sector].shape

def HQHt_sum(self):
    dataSum = 0
    for sector in self.typeNames:
        dataSum = dataSum + self.data[sector]
    return dataSum


def HQ_toFile(self, get_HQFile_name):
    for sector in self.typeNames:
        Sparse_to_nc(self.data[sector], get_HQFile_name(jtime = self.jtime, recepTime = self.recepTime, sectorName = sector))

def HdQ_toFile(self, get_HdQFile_name):
    for sector in self.typeNames:
        nc_write(outFile = get_HdQFile_name(jtime = self.jtime, recepTime = self.recepTime, sectorName = sector), arr = self.data[sector])

def HQHt_toFile(self, get_HQHtFile_name):
    for sector in self.typeNames:
        nc_write(outFile = get_HQHtFile_name(recepTime = self.recepTime, sectorName = sector), arr = self.data[sector])

def INVHQ_toFile(self, get_INVHQFile_name):
    for sector in self.typeNames:
        #Sparse_to_nc(self.data[sector], get_INVHQFile_name(jtime = self.jtime, recepTime = self.recepTime, sectorName = sector))
        nc_write(outFile = get_INVHQFile_name(jtime = self.jtime, recepTime = self.recepTime, sectorName = sector), arr = self.data[sector])

def HiHQti_toFile(self, get_HiHQtiFile_name):
    for sectorName in self.typeNames:
        Sparse_to_nc(spr = self.data[sectorName], FileDirName = get_HiHQtiFile_name(recepTime = self.recepTime, hisTime = self.hisTime, sectorName = sectorName))

