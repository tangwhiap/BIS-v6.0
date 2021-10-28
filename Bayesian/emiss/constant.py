#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

hasInterpolated = True
interpolate_method = None

def interface(time, *args, LON = None, LAT = None, const = 0, **kwargs):

    assert LON is not None and LAT is not None
    assert LON.shape == LAT.shape
    emiss = LON * 0 + const
    return emiss, LON, LAT
