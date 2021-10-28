#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...


from ...utils.distance import areaS
from ...utils.netcdf_io import output_write

import numpy as np
import datetime as dtm

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


def make_output(objExp):

    output_timeScales = objExp.myConfig["output_timeScales"]
    isHourly = "hourly" in output_timeScales
    isDaily = "daily" in output_timeScales
    isWeekly = "weekly" in output_timeScales
    isMonthly = "monthly" in output_timeScales
    isAll = "all" in output_timeScales

    for sector in objExp.sectors:
        timeScale_timeList = {}
        if isHourly:
            emiss_priorHourly = 0
            sigma2_priorHourly = 0
            emiss_postHourly = 0
            sigma2_postHourly = 0
            timeScale_timeList["hourly_timeList"] = []
            time_hourly = objExp.Start
        if isDaily:
            emiss_priorDaily = 0
            sigma2_priorDaily = 0
            emiss_postDaily = 0
            sigma2_postDaily = 0
            timeScale_timeList["daily_timeList"] = []
            time_daily = objExp.Start
        if isWeekly:
            emiss_priorWeekly = 0
            sigma2_priorWeekly = 0
            emiss_postWeekly = 0
            sigma2_postWeekly = 0
            timeScale_timeList["weekly_timeList"] = []
            time_weekly = objExp.Start
        if isMonthly:
            emiss_priorMonthly = 0
            sigma2_priorMonthly = 0
            emiss_postMonthly = 0
            sigma2_postMonthly = 0
            timeScale_timeList["monthly_timeList"] = []
            time_monthly = objExp.Start
        if isAll:
            emiss_priorAll = 0
            sigma2_priorAll = 0
            emiss_postAll = 0
            sigma2_postAll = 0
            time_all = objExp.Start

        matrixS = None
        current = objExp.Start
        while(current <= objExp.End):
            print("postproc: " + current.strftime("%Y-%m-%d %H:%M:%S"))
            priorEmiss, LON, LAT = objExp.get_prior_emiss(time = current, sectorName = sector)
            priorSigma, _, _ = objExp.get_prior_sigma(time = current, sectorName = sector)
            postEmiss, _, _ = objExp.get_posterior_emiss(time = current, sectorName = sector)
            postSigma, _, _ = objExp.get_posterior_sigma(time = current, sectorName = sector)

            if matrixS is None:
                matrixS = areaS(LON, LAT)

            # emiss: umol/m^2/s --> tCO2/cell/hour
            priorEmissTrans = priorEmiss * matrixS * 3600 / 1.0e6 * 44 / 1.0e6
            postEmissTrans = postEmiss * matrixS * 3600 / 1.0e6 * 44 / 1.0e6

            priorSigma2Trans = (priorSigma* matrixS * 3600 / 1.0e6 * 44 / 1.0e6) ** 2
            postSigma2Trans = (postSigma* matrixS * 3600 / 1.0e6 * 44 / 1.0e6) ** 2

            if isHourly:
                emiss_priorHourly = emiss_priorHourly + priorEmissTrans
                emiss_postHourly = emiss_postHourly + postEmissTrans
                sigma2_priorHourly = sigma2_priorHourly + priorSigma2Trans
                sigma2_postHourly = sigma2_postHourly + postSigma2Trans

            if isDaily:
                emiss_priorDaily = emiss_priorDaily + priorEmissTrans
                sigma2_priorDaily = sigma2_priorDaily + priorSigma2Trans
                emiss_postDaily = emiss_postDaily + postEmissTrans
                sigma2_postDaily = sigma2_postDaily + postSigma2Trans

            if isWeekly:
                emiss_priorWeekly = emiss_priorWeekly + priorEmissTrans
                sigma2_priorWeekly = sigma2_priorWeekly + priorSigma2Trans
                emiss_postWeekly = emiss_postWeekly + postEmissTrans
                sigma2_postWeekly = sigma2_postWeekly + postSigma2Trans

            if isMonthly:
                emiss_priorMonthly = emiss_priorMonthly + priorEmissTrans
                sigma2_priorMonthly = sigma2_priorMonthly + priorSigma2Trans
                emiss_postMonthly = emiss_postMonthly + postEmissTrans
                sigma2_postMonthly = sigma2_postMonthly + postSigma2Trans

            if isAll:
                emiss_priorAll = emiss_priorAll + priorEmissTrans
                sigma2_priorAll = sigma2_priorAll + priorSigma2Trans
                emiss_postAll = emiss_postAll + postEmissTrans
                sigma2_postAll = sigma2_postAll + postSigma2Trans

            nextTime = current + objExp.dt

            if isHourly and __isHourly_output(nextTime):

                if (nextTime - objExp.Start) >= dtm.timedelta(hours = 1):

                    outFile = objExp.get_priorHourlyFile_name(time = time_hourly, sectorName = sector)
                    output_write(outFile, emiss_priorHourly, LON, LAT)

                    outFile = objExp.get_spriorHourlyFile_name(time = time_hourly, sectorName = sector)
                    output_write(outFile, np.sqrt(sigma2_priorHourly), LON, LAT)

                    outFile = objExp.get_postHourlyFile_name(time = time_hourly, sectorName = sector)
                    output_write(outFile, emiss_postHourly, LON, LAT)

                    outFile = objExp.get_spostHourlyFile_name(time = time_hourly, sectorName = sector)
                    output_write(outFile, np.sqrt(sigma2_postHourly), LON, LAT)

                    timeScale_timeList["hourly_timeList"].append(time_hourly)

                emiss_priorHourly = 0
                sigma2_priorHourly = 0
                emiss_postHourly = 0
                sigma2_postHourly = 0
                time_hourly = nextTime

            if isDaily and __isDaily_output(nextTime):

                if (nextTime - objExp.Start) >= dtm.timedelta(days = 1):

                    outFile = objExp.get_priorDailyFile_name(time = time_daily, sectorName = sector)
                    output_write(outFile, emiss_priorDaily, LON, LAT)

                    outFile = objExp.get_spriorDailyFile_name(time = time_daily, sectorName = sector)
                    output_write(outFile, np.sqrt(sigma2_priorDaily), LON, LAT)

                    outFile = objExp.get_postDailyFile_name(time = time_daily, sectorName = sector)
                    output_write(outFile, emiss_postDaily, LON, LAT)

                    outFile = objExp.get_spostDailyFile_name(time = time_daily, sectorName = sector)
                    output_write(outFile, np.sqrt(sigma2_postDaily), LON, LAT)

                    timeScale_timeList["daily_timeList"].append(time_daily)

                emiss_priorDaily = 0
                sigma2_priorDaily = 0
                emiss_postDaily = 0
                sigma2_postDaily = 0
                time_daily = nextTime

            if isWeekly and __isWeekly_output(nextTime):

                if (nextTime - objExp.Start) >= (dtm.timedelta(days = 7)):

                    outFile = objExp.get_priorWeeklyFile_name(time = time_weekly, sectorName = sector)
                    output_write(outFile, emiss_priorWeekly, LON, LAT)

                    outFile = objExp.get_spriorWeeklyFile_name(time = time_weekly, sectorName = sector)
                    output_write(outFile, np.sqrt(sigma2_priorWeekly), LON, LAT)

                    outFile = objExp.get_postWeeklyFile_name(time = time_weekly, sectorName = sector)
                    output_write(outFile, emiss_postWeekly, LON, LAT)

                    outFile = objExp.get_spostWeeklyFile_name(time = time_weekly, sectorName = sector)
                    output_write(outFile, np.sqrt(sigma2_postWeekly), LON, LAT)

                    timeScale_timeList["weekly_timeList"].append(time_weekly)

                emiss_priorWeekly = 0
                sigma2_priorWeekly = 0
                emiss_postWeekly = 0
                sigma2_postWeekly = 0
                time_weekly = nextTime

            if isMonthly and __isMonthly_output(nextTime):

                if (nextTime - objExp.Start) >= (dtm.timedelta(days = 28)):

                    outFile = objExp.get_priorMonthlyFile_name(time = time_monthly, sectorName = sector)
                    output_write(outFile, emiss_priorMonthly, LON, LAT)

                    outFile = objExp.get_spriorMonthlyFile_name(time = time_monthly, sectorName = sector)
                    output_write(outFile, np.sqrt(sigma2_priorMonthly), LON, LAT)

                    outFile = objExp.get_postMonthlyFile_name(time = time_monthly, sectorName = sector)
                    output_write(outFile, emiss_postMonthly, LON, LAT)

                    outFile = objExp.get_spostMonthlyFile_name(time = time_monthly, sectorName = sector)
                    output_write(outFile, np.sqrt(sigma2_postMonthly), LON, LAT)

                    timeScale_timeList["monthly_timeList"].append(time_monthly)

                emiss_priorMonthly = 0
                sigma2_priorMonthly = 0
                emiss_postMonthly = 0
                sigma2_postMonthly = 0
                time_monthly = nextTime

            current += objExp.dt

        if isAll:
            
            outFile = objExp.get_priorAllFile_name(sectorName = sector)
            output_write(outFile, emiss_priorAll, LON, LAT)

            outFile = objExp.get_spriorAllFile_name(sectorName = sector)
            output_write(outFile, np.sqrt(sigma2_priorAll), LON, LAT)

            outFile = objExp.get_postAllFile_name(sectorName = sector)
            output_write(outFile, emiss_postAll, LON, LAT)

            outFile = objExp.get_spostAllFile_name(sectorName = sector)
            output_write(outFile, np.sqrt(sigma2_postAll), LON, LAT)

    return timeScale_timeList