#!/usr/bin/env python


import numpy as np
from pdb import set_trace

def areaS(LON, LAT):
    # LON, LAT unit: degree
    # S unit: km^2
    R = 6371
    dlon = (LON[0, 1:] - LON[0, :-1]).mean()
    dlat = (LAT[1:, 0] - LAT[:-1, 0]).mean()
    S = (R**2) * np.deg2rad(dlon) * (np.sin(np.deg2rad(LAT + dlat/2 )) - np.sin(np.deg2rad(LAT - dlat/2 )))
    assert S.shape == LON.shape
    assert S.shape == LAT.shape
    return S


def resample_points(pData, pLon, pLat, LON, LAT):


    lon_s = LON.min()
    lon_e = LON.max()
    lat_s = LAT.min()
    lat_e = LAT.max()
    nLat, nLon = LON.shape
    dlon = (lon_e - lon_s) / (nLon - 1)
    dlat = (lat_e - lat_s) / (nLat - 1)

    indexWithin = (pLon >= lon_s) & (pLon <= lon_e) & (pLat >= lat_s) & (pLat <= lat_e)
    indexWithout = ~indexWithin
    if np.sum(indexWithout) > 0:
        print("Warning! There are " + str(np.sum(indexWithout)) + " point source is out of range")
    
    pLon = pLon[indexWithin]
    pLat = pLat[indexWithin]
    pData = pData[indexWithin]
    
    ixList = (pLon - lon_s + dlon/2)/dlon
    iyList = (pLat - lat_s + dlat/2)/dlat
    ixList = ixList.astype(np.int)
    iyList = iyList.astype(np.int)

    emissField = np.zeros_like(LON)
    for ix, iy, iData in zip(ixList, iyList, pData):
        emissField[iy, ix] += iData

    return emissField


