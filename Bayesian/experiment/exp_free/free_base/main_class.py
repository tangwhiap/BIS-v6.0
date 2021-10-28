#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 09/2021
#   ...

from ...get_freeBaseClass import get_freeBaseClass
from ....main.configure import EXPTYPE, BASEEXP

baseExpClass = get_freeBaseClass(EXPTYPE, BASEEXP)

class ExpFree(baseExpClass):
    pass

"""
def get_expFreeClass(EXPNAME, EXPTYPE, BASEEXP):

    baseExpClass = get_freeBaseClass(EXPNAME, EXPTYPE, BASEEXP)

    class ExpFree(baseExpClass):
        pass

    return baseExpClass
"""
