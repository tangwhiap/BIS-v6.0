#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 09/2021
#   ...

#from ..main.configure import ISFREE, EXPTYPE, BASEEXP

def get_freeBaseClass(EXPTYPE, BASEEXP):

    loc = locals()
    exec("from ." + EXPTYPE + "." + BASEEXP + ".main_class import Exp" + EXPTYPE + "_" + BASEEXP + " as exp")
    baseExpClass = loc["exp"]

    return baseExpClass

