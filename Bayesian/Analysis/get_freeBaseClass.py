#!/usr/bin/env python


# Authors:
#   Wenhan TANG - 09/2021
#   ...

def get_freeBaseClass(EXPTYPE, BASEEXP):

    loc = locals()
    exec("from ." + EXPTYPE + ".benchplots_" + BASEEXP + ".main_class import BenchPlots_" + BASEEXP + " as exp")

    baseExpClass = loc["exp"]
    return baseExpClass


