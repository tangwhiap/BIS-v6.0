#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 09/2021
#   ...

def get_expClass(EXPNAME, EXPTYPE, ISFREE = False, BASEEXP = None):

    if ISFREE:
        loc = locals()
        exec("from .exp_free." + EXPNAME + ".main_class import get_" + EXPNAME + " as get_expClass")
        get_expClass = loc["get_expClass"]
        expClass = get_expClass(EXPNAME, EXPTYPE, BASEEXP)

    else:
        loc = locals()
        exec("from ." + EXPTYPE + "." + EXPNAME + ".main_class import Exp" + EXPTYPE + "_" + EXPNAME + " as exp")
        expClass = loc["exp"]

    return expClass
