#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...


def get_benchClass(expName, expType, isFree = False):
    assert expType.upper() in ["REAL", "OSSE"]
    if isFree:
        loc = locals()
        cmd = "from .FREE.benchplots_" + expName + ".main_class import BenchPlots_" + expName + " as BenchClass"
        exec(cmd)
        BenchClass = loc["BenchClass"]
    else:
        loc = locals()
        cmd = "from ." + expType.upper() + ".benchplots_" + expName + ".main_class import BenchPlots_" + expName + " as BenchClass"
        exec(cmd)
        BenchClass = loc["BenchClass"]
    return BenchClass
