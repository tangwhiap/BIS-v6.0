#!/usr/bin/env python

#-> 4D-Var bayesian inversion algorithm <-#

# Authors:
#    Wenhan TANG - 09/2021
#    ...

from ..core_optimize.main_class import CoreOptimize
from ..core_optimize.iteration import optimize_iteration

import datetime as dtm
from pdb import set_trace

class Core4DVar(CoreOptimize):

    def __init__(self, *args, **kwargs):

        CoreOptimize.__init__(self, *args, **kwargs)
        self.methodName = "4DVar"
        coreConfig = kwargs["coreConfig"]
        self.myConfig = {**self.myConfig, **coreConfig["4DVar"]}
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

    def relevant_trange(self, timeJ):
        relevant_length = self.myConfig["relevant_length"]
        unitT = dtm.timedelta(hours = 1)
        start = timeJ - relevant_length * unitT
        end = timeJ + relevant_length * unitT
        start = start if start >= self.start else self.start
        end = end if end <= self.end else self.end
        dt = dtm.timedelta(hours = self.dtHrs)
        trange = []
        time = start
        while(time <= end):
            trange.append(time)
            time += dt
        #print("timeJ = ", timeJ)
        #print("trange = ", trange)
        return trange

    def compute_HQ(self, recepTime, objExp):
        #tj_range = self.backtime_j2i(recepTime)
        tj_range = self.relevant_trange(recepTime)
        if self.myConfig["isCompute_dQ"]:
            CoreOptimize.compute_HQ(self, objExp = objExp, recepTime = recepTime, tj_range = tj_range, backtime_j2i = self.backtime_j2i, Stype = "Prior")
        else:
            CoreOptimize.compute_HQ(self, objExp = objExp, recepTime = recepTime, tj_range = tj_range, backtime_j2i = self.backtime_j2i, Stype = "Proc")

    def compute_HQHt(self, recepTime, objExp):
        CoreOptimize.compute_HQHt(self, objExp = objExp, recepTime = recepTime, backtime_j2i = self.backtime_j2i)

    def compute_HdQ(self, recepTime, objExp):
        tj_range = self.relevant_trange(recepTime)
        CoreOptimize.compute_HdQ(self, objExp = objExp, recepTime = recepTime, tj_range = tj_range, backtime_j2i = self.backtime_j2i)

    def inversion(self, recepTime, objExp):
        #tj_range = self.backtime_j2i(recepTime)
        tj_range = self.relevant_trange(recepTime)
        CoreOptimize.inversion(self, objExp = objExp, recepTime = recepTime, tj_range = tj_range)

    def optimize(self, recepTime, objExp):
        #CoreOptimize.optimize(self, recepTime = recepTime, objExp = objExp)
        self.make_workdir(recepTime)
        if self.myConfig["isCompute_dQ"]:
            self.compute_HdQ(recepTime = recepTime, objExp = objExp)
        self.compute_HQ(recepTime = recepTime, objExp = objExp)
        self.compute_HQHt(recepTime = recepTime, objExp = objExp)
        self.compute_HQHt_R(recepTime = recepTime, objExp = objExp)
        self.compute_d(recepTime = recepTime, objExp = objExp)
        if objExp.hasBCK:
            self.background_optimize(recepTime = recepTime, objExp = objExp, optTimeList = [recepTime])
            self.compute_d(recepTime = recepTime, objExp = objExp, bckPrior = False)
        #self.compute_INV(recepTime = recepTime, objExp = objExp)
        self.inversion(recepTime = recepTime, objExp = objExp)
        
