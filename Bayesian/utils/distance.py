#!/usr/bin/env python
# Authors: Wenhan TANG - 06/2021

import numpy as np

from pdb import set_trace

def radians(x):
    return x * np.pi / 180.

def DIS_point_to_points(plon, plat, pslon, pslat, iy = None, ix = None):
    R = 6371 #unit = KM
    Npoints = len(pslon)
    assert Npoints == len(pslat)
    rplon = radians(plon)
    rplat = radians(plat)
    rpslon = radians(pslon)
    rpslat = radians(pslat)
    #if iy == 56 and ix == 0:
    #   set_trace()
    cosl = np.cos(rplat) * np.cos(rpslat) * np.cos(rplon - rpslon) + np.sin(rplat) * np.sin(rpslat) 
    cosl = np.where(cosl > 1, 1, cosl)
    cosl = np.where(cosl < -1, -1, cosl)
    return R * np.arccos(cosl)

def areaS(LON, LAT):
    # LON, LAT unit: degree
    # S unit: m^2
    R = 6371000
    dlon = (LON[0, 1:] - LON[0, :-1]).mean()
    dlat = (LAT[1:, 0] - LAT[:-1, 0]).mean()
    S = (R**2) * np.deg2rad(dlon) * (np.sin(np.deg2rad(LAT + dlat/2 )) - np.sin(np.deg2rad(LAT - dlat/2 )))
    assert S.shape == LON.shape
    assert S.shape == LAT.shape
    return S
