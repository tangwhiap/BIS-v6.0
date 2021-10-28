#!/usr/bin/env python

#-> 3D-Var bayesian inversion algorithm <-#

# Authors:
#    Wenhan TANG - 07/2021
#    ...

from ..core_optimize.main_class import CoreOptimize
from ..core_optimize.iteration import optimize_iteration

import datetime as dtm
from pdb import set_trace

class Core3DVar(CoreOptimize):

    def __init__(self, *args, **kwargs):

        CoreOptimize.__init__(self, *args, **kwargs)
        self.methodName = "3DVar"
        coreConfig = kwargs["coreConfig"]
        self.myConfig = {**self.myConfig, **coreConfig["3DVar"]}
        #self.obsTimeList = self.get_obsTimeList()
        #self.backTimeList = self.backtime_j2i(self.recepTime)
        Iteration = optimize_iteration(self.myConfig["iteration"])
        self.objIter = Iteration(ExpCore = self, **self.myConfig["iteration_kwargs"])
        #for time in self.objIter:
        #    print(time)

    def backtime_j2i(self, timeJ):
        nBack = self.myConfig["nBack"]
        dt = dtm.timedelta(hours = self.dtHrs)
        start = timeJ - dt * nBack
        end = timeJ - dt * 1
        backTimeList = []
        time = start 
        while(time <= end):
            backTimeList.append(time)
            time += dt
        return backTimeList

    def get_obsTimeList(self, start = None, end = None, window = None, nBack = None):
        if start is None:
            start = self.start
        if end is None:
            end = self.end
        if window is None:
            window = self.myConfig["window"]
        if nBack is None:
            nBack = self.myConfig["nBack"]
        obsTimeList = []
        time = start + dtm.timedelta(hours = nBack)
        while(time <= end):
            obsTimeList.append(time)
            time += dtm.timedelta(hours = window)
        return obsTimeList

    def compute_HQ(self, recepTime, objExp):
        tj_range = self.backtime_j2i(recepTime)
        CoreOptimize.compute_HQ(self, objExp = objExp, recepTime = recepTime, tj_range = tj_range, backtime_j2i = self.backtime_j2i)

    def compute_HQHt(self, recepTime, objExp):
        CoreOptimize.compute_HQHt(self, objExp = objExp, recepTime = recepTime, backtime_j2i = self.backtime_j2i)

    def inversion(self, recepTime, objExp):
        tj_range = self.backtime_j2i(recepTime)
        CoreOptimize.inversion(self, objExp = objExp, recepTime = recepTime, tj_range = tj_range)

    def optimize(self, recepTime, objExp):
        #CoreOptimize.optimize(self, recepTime = recepTime, objExp = objExp)
        self.make_workdir(recepTime)
        self.compute_HQ(recepTime = recepTime, objExp = objExp)
        self.compute_HQHt(recepTime = recepTime, objExp = objExp)
        self.compute_HQHt_R(recepTime = recepTime, objExp = objExp)
        self.compute_d(recepTime = recepTime, objExp = objExp)
        if objExp.hasBCK:
            self.background_optimize(recepTime = recepTime, objExp = objExp, optTimeList = [recepTime])
            self.compute_d(recepTime = recepTime, objExp = objExp, bckPrior = False)
        #self.compute_INV(recepTime = recepTime, objExp = objExp)
        self.inversion(recepTime = recepTime, objExp = objExp)
        
