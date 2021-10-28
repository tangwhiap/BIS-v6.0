#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

import numpy as np
import datetime as dtm

def optimize_iteration(iterName):
    loc = locals()
    exec("Iteration = Iteration_" + iterName)
    Iteration = loc["Iteration"]
    return Iteration

class Iteration(object):
    def __init__(self):
        self.iterTimeList = []
        self.obsTimeList = []

    def __iter__(self):
        return self.__next__()

    def __next__(self):
        return self.__generator()

    def __generator(self):
        for i in self.iterList:
            yield i

class Iteration_Window(Iteration):

    def __init__(self, ExpCore, *args, **kwargs):
        self.iterTimeList = self.get_iterTimeList(ExpCore, *args, **kwargs)
        self.obsTimeList = self.get_obsTimeList()

        self.iterList = self.iterTimeList

    def get_iterTimeList(self, ExpCore, start = None, end = None, dt_obs = None, nBack = None, **kwargs):

        if start is None:
            start = ExpCore.start

        if end is None:
            end = ExpCore.end

        if dt_obs is None:
            #dt_obs = ExpCore.myConfig["dt_obs"]
            dt_obs = 1

        if nBack is None:
            nBack = ExpCore.myConfig["nBack"]

        iterTimeList = []
        time = start + dtm.timedelta(hours = nBack)

        while(time <= end):
            iterTimeList.append(time)
            time += dtm.timedelta(hours = dt_obs)

        return iterTimeList

    def get_obsTimeList(self):
        return self.iterTimeList

class Iteration_specific_hour(Iteration):

    def __init__(self, ExpCore, *args, **kwargs):
        self.iterTimeList = self.get_iterTimeList(ExpCore, *args, **kwargs)
        self.obsTimeList = self.get_obsTimeList()

        self.iterList = self.iterTimeList

    def get_iterTimeList(self, ExpCore, start = None, end = None, nBack = None, hourList = None, UTC = 0):

        

        assert hourList is not None

        if start is None:
            start = ExpCore.start

        if end is None:
            end = ExpCore.end

        #if dt_obs is None:
        #    dt_obs = ExpCore.myConfig["dt_obs"]

        if nBack is None:
            nBack = ExpCore.myConfig["nBack"]

        hourList_utc0 = (np.array(hourList).astype("int") - UTC).tolist()
        iterTimeList = []
        time = start + dtm.timedelta(hours = nBack)

        while(time <= end):
            if time.hour  in hourList_utc0:
                iterTimeList.append(time)

            time += dtm.timedelta(hours = 1)

        return iterTimeList

    def get_obsTimeList(self):
        return self.iterTimeList

class Iteration_3DVarWindow(Iteration_Window):
    pass

class Iteration_4DVarWindow(Iteration_Window):
    pass

class Iteration_3DVar_specific_hour(Iteration_specific_hour):
    pass

class Iteration_4DVar_specific_hour(Iteration_specific_hour):
    pass

#-- 4DVarWindow is as same as 3DVarWindow in the current version --#
#Iteration_4DVarWindow = Iteration_3DVarWindow

#-- 4DVar_specific_hour is as same as 3DVar_specific_hour in the current version --#
#Iteration_4DVar_specific_hour = Iteration_3DVar_specific_hour

