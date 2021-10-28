#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...

import types

def module2dic(module):
    subList = dir(module)
    dic = {}
    for subName in subList:
        if subName[:2] == "__" and subName[-2:] == "__":
            continue
        loc = locals()
        exec("isModule = isinstance(module." + subName + ", types.ModuleType)")
        isModule = loc["isModule"]
        if isModule:
            print(subName)
            continue
        loc = locals()
        exec("sub = module." + subName)
        sub = loc["sub"]
        dic[subName] = sub
    return dic
