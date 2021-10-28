#!/usr/bin/env python
# Authors:
#   Wenhan TANG - 06/2021 (Original version)
#   Wenhan TANG - 07/2021 (Compatible to BIS_v6.0 later)
#   ...

import numpy as np
from scipy.sparse import csr_matrix
import netCDF4 as nc
import matplotlib.pyplot as plt
from pdb import set_trace
import warnings
import time

warnings.filterwarnings('ignore')

Spr_max = 10000

class SparseMatrix(object):

    def __init__(self, DataDir):

        self.ncf = nc.Dataset(DataDir, "r")
        self.data = self.ncf.variables["data"]
        self.indices = self.ncf.variables["indices"]
        self.indptr = self.ncf.variables["indptr"]
        self.nValid = self.ncf.dimensions["DataDim"].size
        #self.Ns = Ns
        #self.Ng = Ng
        self.ny = self.ncf.ny
        self.nx = self.ncf.nx
        self.shape = (self.ny, self.nx)
        #self.nlon = nlon
        #self.nlat = nlat
        #self.max_N = max_N
    
    def getLine(self, iline):
        ptr1 = self.indptr[iline]
        ptr2 = self.indptr[iline + 1]
        colind = self.indices[ptr1:ptr2]
        vdata = self.data[ptr1:ptr2]
        arr = np.zeros(self.nx)
        arr[colind] = vdata
        return arr
   
    def draw_relation(self, *args, **kwargs):
        nlon, nlat = kwargs["nlon"], kwargs["nlat"]
        if len(args) == 1:
            iline = args[0]
        elif len(args) == 2:
            iy, ix = args
            iline = iy * nlat + ix
        else:
            assert False
        relation = self.getLine(iline).reshape(nlat, nlon)
        plt.contourf(relation)
        plt.colorbar()
        plt.show()

    def SparseOut(self):
        return csr_matrix( (self.data[:].filled(np.nan), self.indices[:].filled(np.nan), self.indptr[:].filled(np.nan)) , shape = (self.ny, self.nx))

    def SparseLineSel(self, start, stop):
        #set_trace()
        indptr = self.indptr[start:stop + 1].filled(np.nan)
        #indptr -= indptr[0]
        ptr1 = indptr[0]
        ptr2 = indptr[-1]
        #print(" -->", ptr1, " ", ptr2, end = "")
        data = self.data[ptr1:ptr2].filled(np.nan)
        indices = self.indices[ptr1:ptr2].filled(np.nan)
        #print(" len =", len(data))
        return csr_matrix( (data, indices, indptr - indptr[0]) , shape = (len(indptr) - 1, self.nx))

    def __rmul__(self, other):
        return Sparse_mul(other, self)
    
    def close(self):
        self.ncf.close()
    
    #def __del__(self):
    #    self.close()

def create_diag(diag):
    N = len(diag)
    data = diag
    indices = np.arange(N)
    indptr = np.arange(N + 1)
    return csr_matrix( (data, indices, indptr), shape = (N, N) )

def array_to_sparse_to_nc(array, FileDirName):
    assert isinstance(array, np.ndarray)
    Sparse_to_nc(csr_matrix(array), FileDirName)

def Sparse_to_nc(spr, FileDirName):
    ncf = nc.Dataset(FileDirName, "w", format = "NETCDF4")
    ny, nx = spr.shape
    DLen = len(spr.data)
    ncf.createDimension("DataDim", None)
    ncf.createDimension("PtrDim", ny + 1)
    ncf.createVariable("data", "f", ("DataDim"))[:] = spr.data
    ncf.createVariable("indices", "i", ("DataDim"))[:] = spr.indices
    ncf.createVariable("indptr", "i", ("PtrDim"))[:] = spr.indptr
    ncf.setncattr("ny", ny)
    ncf.setncattr("nx", nx)
    ncf.close()

def getNCSPR(FileName):
    A = SparseMatrix(FileName)
    return A
    # Don't forget to A.close(), thanks~

def getNCSPR_spr(FileName):
    A = SparseMatrix(FileName)
    spr_A = A.SparseOut()
    A.close()
    return spr_A

def getNCSPR_array(FileName):
    A = SparseMatrix(FileName)
    data = A.SparseOut().toarray()
    A.close()
    return data


def Spr_zeros(ny, nx):
    return csr_matrix(([],[],np.zeros(ny + 1)), shape = (ny, nx))

def Sparse_mul(spr_A, B):#, right_Diag = None):

    assert spr_A.shape[1] == B.ny
    
    Ndiv = B.nValid // Spr_max + 1

    #print("Divided into ", Ndiv)
    MulDim = B.ny
    PerLen = MulDim // Ndiv + 1
    assert Ndiv <= MulDim, "The value of Spr_max = " + str(Spr_max) + " is too small."
    #print("PerLen = ", PerLen)
    #if right_Diag is not None:
    #    assert len(right_Diag) == MulDim
    idiv = 0

    #for idiv in range(Ndiv):
    while(True):

        start = idiv * PerLen
        stop = (idiv + 1) * PerLen

        assert start < MulDim
        if stop > MulDim:
            stop = MulDim
        #print("--> ", start, stop)
        #print("Processing: ", idiv, " (", start, ".", stop, ") ",end = '')
        temp = spr_A[:, start:stop] * B.SparseLineSel(start = start, stop = stop)
        if idiv == 0:
            MulOut = temp
        else:
            MulOut += temp
        if stop >= MulDim:
            break
        idiv += 1

    return MulOut
