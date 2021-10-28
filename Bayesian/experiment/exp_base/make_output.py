#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...

from ...utils.distance import areaS
from ...utils.netcdf_io import output_write

import numpy as np
import datetime as dtm


from pdb import set_trace

def __isHourly_output(time):
    if time.minute == 0 and time.second == 0:
        return True
    return False

def __isDaily_output(time):
    if __isHourly_output(time):
        if time.hour == 0:
            return True
    return False

def __isWeekly_output(time):
    if __isDaily_output(time):
        if time.weekday() == 0:
            return True
    return False

def __isMonthly_output(time):
    if __isDaily_output(time):
        if time.day == 1:
            return True
    return False

def __is_output(time, start, timeScale):

    if timeScale.lower() == "hourly":
        is_beginTime = __isHourly_output(time)
        is_longerThan_oneStep = (time - start) >= dtm.timedelta(hours = 1)
        return is_beginTime, is_longerThan_oneStep

    if timeScale.lower() == "daily":
        is_beginTime = __isDaily_output(time)
        is_longerThan_oneStep = (time - start) >= dtm.timedelta(hours = 24)
        return is_beginTime, is_longerThan_oneStep

    if timeScale.lower() == "weekly":
        is_beginTime = __isWeekly_output(time)
        is_longerThan_oneStep = (time - start) >= dtm.timedelta(days = 7)
        return is_beginTime, is_longerThan_oneStep

    if timeScale.lower() == "monthly":
        is_beginTime = __isMonthly_output(time)
        is_longerThan_oneStep = (time - start) >= dtm.timedelta(days = 28)
        return is_beginTime, is_longerThan_oneStep

    if timeScale.lower() == "all":
        return False, False

def make_output(objExp):

    make_output_dic = objExp.make_output_dic
    output_timeScales = objExp.myConfig["output_timeScales"]
    isHourly = "hourly" in output_timeScales
    isDaily = "daily" in output_timeScales
    isWeekly = "weekly" in output_timeScales
    isMonthly = "monthly" in output_timeScales
    isAll = "all" in output_timeScales

    outputType = objExp.myConfig["outputType"].lower()
    assert outputType in ["sum", "mean"] 
    isSum = outputType == "sum"
    isMean = outputType == "mean"

    """
    sumData = {}
    for variableName in objExp.make_output_dic:
        sumData[variableName] = {}
        for timeScale in output_timeScales:
            sumData[variableName][timeScale] = 0
    """

    for sector in objExp.sectors:
        timeScale_timeList = {}
        time_dic = {}
        sumData = {}
        numData = {}
        for timeScale in output_timeScales:
            sumData[timeScale] = {}
            numData[timeScale] = {}
            for variableName in make_output_dic:
                sumData[timeScale][variableName] = 0
                numData[timeScale][variableName] = 0

            timeScale_timeList[timeScale] = []
            time_dic[timeScale] = objExp.Start

        current = objExp.Start
        while(current <= objExp.End):
            print("postproc: " + current.strftime("%Y-%m-%d %H:%M:%S")) 
            variable_thisTime_dic = {}
            LONLAT_thisTime_dic = {}
            for variableName in make_output_dic:
                array, LON, LAT = make_output_dic[variableName]["inputFun"](time = current, sectorName = sector)
                variable_thisTime_dic[variableName] = make_output_dic[variableName]["transFun"](array, LON, LAT)
                LONLAT_thisTime_dic[variableName] = LON, LAT

            for timeScale in output_timeScales:
                for variableName in make_output_dic:
                    sumData[timeScale][variableName] += variable_thisTime_dic[variableName]
                    numData[timeScale][variableName] += 1

            nextTime = current + objExp.dt
            
            for timeScale in output_timeScales:

                is_beginTime, is_longerThan_oneStep =  __is_output(nextTime, start = objExp.Start, timeScale = timeScale)

                if is_beginTime:
                    if is_longerThan_oneStep:
                        for variableName in make_output_dic:
                            outArr = make_output_dic[variableName]["outFun"](sumData[timeScale][variableName])
                            if isMean:
                                outArr = outArr / numData[timeScale][variableName]
                            outFile = make_output_dic[variableName]["outNameFun_" + timeScale](time = time_dic[timeScale], sectorName = sector)
                            #if timeScale == "daily":
                            #    set_trace()
                            output_write(outFile, outArr, LONLAT_thisTime_dic[variableName][0], LONLAT_thisTime_dic[variableName][1])
                            timeScale_timeList[timeScale].append(time_dic[timeScale])
                    
                    for variableName in make_output_dic:
                        sumData[timeScale][variableName] = 0
                        numData[timeScale][variableName] = 0

                    time_dic[timeScale] = nextTime

            current += objExp.dt

        if isAll:
            for variableName in make_output_dic:
                outArr = make_output_dic[variableName]["outFun"](sumData["all"][variableName])
                outFile = make_output_dic[variableName]["outNameFun_all"](sectorName = sector)
                output_write(outFile, outArr, LONLAT_thisTime_dic[variableName][0], LONLAT_thisTime_dic[variableName][1])

                        



