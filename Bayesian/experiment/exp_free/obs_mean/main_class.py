#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 09/2021
#   ...


from ..free_base.main_class import ExpFree
from ....core.core_optimize.iteration import Iteration_specific_hour
from ....utils.sparse_matrix import array_to_sparse_to_nc

import numpy as np
import datetime as dtm
from pdb import set_trace

myName = "obs_mean"
class ExpFree_obs_mean(ExpFree):

    def __init__(self, *args, **kwargs):

        config = kwargs["configure"]
        UTC = config["ExpConfig"][myName]["obs_mean_UTC"]
        self.obs_mean_start_hour = config["ExpConfig"][myName]["mean_start_hour"] - UTC
        self.obs_mean_end_hour = config["ExpConfig"][myName]["mean_end_hour"] - UTC

        if self.obs_mean_start_hour < 0:
            self.obs_mean_start_hour += 24
        if self.obs_mean_start_hour > 23:
            self.obs_mean_start_hour -= 24

        if self.obs_mean_end_hour < 0:
            self.obs_mean_end_hour += 24
        if self.obs_mean_end_hour > 23:
            self.obs_mean_end_hour -= 24

        if self.obs_mean_start_hour <= self.obs_mean_end_hour:
            self.mean_length = self.obs_mean_end_hour - self.obs_mean_start_hour

        if self.obs_mean_start_hour > self.obs_mean_end_hour:
            self.mean_length = 24 - (self.obs_mean_start_hour - self.obs_mean_end_hour)

        #coreName = config["ExpConfig"]["core"]
        backwardTime = config["CoreConfig"]["optimize"]["nBack"]
        kwargs["configure"]["CoreConfig"]["optimize"]["nBack"] = backwardTime + self.mean_length

        print("Change nBack from " + str(backwardTime) + " to " + str(kwargs["configure"]["CoreConfig"]["optimize"]["nBack"]) + ".")
        ExpFree.__init__(self, *args, **kwargs)
        self.myConfig = dict(**self.myConfig, **config["ExpConfig"][myName])
        

        assert isinstance(self.objCore.objIter, Iteration_specific_hour), "\nSetting Error: objIter must be an object belongs to Iteration_specific_hour class."

        assert self.check_iterTimeList(self.objCore.objIter.iterTimeList), "\nSetting Error: objIter's hourList must have only one number which is the end hour of observation mean."

        print("OBS mean time range from " + str(self.obs_mean_start_hour) + " o'clock (UTC) to " + str(self.obs_mean_end_hour) + " o'clock (UTC).")


    def check_iterTimeList(self, timeList):
        isPass = True
        for time in timeList:
            if time.hour != self.obs_mean_end_hour:
                isPass = False
                break
        return isPass

    def process_H(self, kwargs):

        assert "time" in kwargs
        assert "recepTime" in kwargs
        assert "Htype" in kwargs
        assert "objFoot" in kwargs
        assert "getHFileName" in kwargs

        time = kwargs["time"]
        recepTime = kwargs["recepTime"]
        Htype = kwargs["Htype"]
        objFoot = kwargs["objFoot"]
        get_HFile_name = kwargs["getHFileName"]

        obs_mean_start_hour = self.obs_mean_start_hour
        obs_mean_end_hour = self.obs_mean_end_hour
        mean_length = self.mean_length

        assert recepTime.hour == obs_mean_end_hour

        recepTime_range_start = recepTime - dtm.timedelta(hours = mean_length)
        assert recepTime_range_start.hour == obs_mean_start_hour

        recepTime_range = []
        iRecepTime = recepTime_range_start
        while(iRecepTime <= recepTime):
            recepTime_range.append(iRecepTime)
            iRecepTime += self.dt

        arrayH = []
        for receptor in objFoot.receptors:
            foot_sum = 0
            for iRecepTime in recepTime_range:
                foot, _, _ = objFoot.get_footprint(time = time, recepTime = recepTime, receptorName = receptor)
                #print("Warning!! footprint not found. time = " + time.strftime("%Y-%m-%d %H:%M:%S") + ", recepTime = " + recepTime.strftime("%Y-%m-%d %H:%M:%S") + ", receptor = " + receptor) 
                foot_sum = foot_sum + foot
            foot_mean = foot_sum / len(recepTime_range)
            arrayH.append(foot_mean.flatten())
        arrayH = np.array(arrayH)
        array_to_sparse_to_nc(arrayH, get_HFile_name(time = time, recepTime = recepTime, typeName = Htype))





