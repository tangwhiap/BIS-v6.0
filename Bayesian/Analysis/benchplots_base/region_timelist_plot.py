#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...

from .draw_object_class import AnalysisTimeSeries, AnalysisLinePlots, ds2ATS
from ..Analyst.local_point import LocPoint

from ..utils.Bench import BenchPlots_png as BenchPlots
from ..utils.Bench import multicases_subplot_scheme as mcss

from pdb import set_trace



def pointRegion_timeSeries(objBench, drawStart, drawEnd, timeScale):

    
    myConfig = objBench.myConfig["point_region_timeSeries_kwargs"]
    pointRegionList = objBench.myConfig["pointRegionList"]
    pointRegion_kwargs = objBench.myConfig["point_region_kwargs"]
    point_vrange = objBench.myConfig["point_vrange"]
    region_vrange = objBench.myConfig["region_vrange"]

    #objBench.build_point_region_timeSeries_dataset(pointRegionList = pointRegionList, drawStart = drawStart, drawEnd = drawEnd, timeScale = timeScale)
    pointRegion_timeSeries_dic = objBench.build_pointRegion_timeSeries_dic(pointRegionList = pointRegionList, drawStart = drawStart, drawEnd = drawEnd, timeScale = timeScale)

    point_timeSeries_dic = {}
    region_timeSeries_dic = {}

    for pointRegionName in pointRegion_timeSeries_dic:

        if pointRegion_kwargs[pointRegionName]["type"].lower() == "point":
            point_timeSeries_dic[pointRegionName] = pointRegion_timeSeries_dic[pointRegionName]

        elif pointRegion_kwargs[pointRegionName]["type"].lower() == "region":
            region_timeSeries_dic[pointRegionName] = pointRegion_timeSeries_dic[pointRegionName]

        else:
            assert False



    bcp = BenchPlots(figsize = myConfig["figsize"], save_name = objBench.benchDir + "/" + myConfig["BenchName"] + "_" + timeScale)
    #objMcss = mcss(Nrow_max = 4, Ncol_max = 1, span_dim = "row", auto_reverse = False, remain_adjust = False)
    #objMcss.layout(len(region_timeSeries_dic))

    Nrow = 4
    hspace = 0.2 

    axs = bcp.subplots(nrows = Nrow, ncols = 1, sharex = True)
    bcp.cfig.subplots_adjust(hspace = hspace) 
    axi = 0

    for regionName in region_timeSeries_dic:

        print("Drawing -> " + regionName)

        if axi >= Nrow:
            bcp.cfig.autofmt_xdate()
            bcp.cfig.suptitle("Region time series") 
            bcp.next_page()
            axs = bcp.subplots(nrows = Nrow, ncols = 1, sharex = True)
            bcp.cfig.subplots_adjust(hspace = hspace) 
            axi = 0


        ax = axs[axi]
        #tS_dic = region_timeSeries_dic[regionName]
        timeUnit_dic = objBench.myConfig["timeUnit_dic"] #if objBench.isSum else "umol/m^2/s"
        LinePlots_kwargs = {**myConfig[regionName + "_kwargs"], "vrange": region_vrange[timeScale if objBench.isSum else "mean"], "unit": ("tCO2/" + timeUnit_dic[timeScale]) if objBench.isSum else "umol/m^2/s", "title": regionName + " city region"}
        objALP = AnalysisLinePlots(region_timeSeries_dic[regionName], LinePlots_kwargs = LinePlots_kwargs)
        ax = objALP.plot(ax = ax)
        if objBench.isSum:
            resultList = objALP.compute_total({"Prior": "priorEmiss", "Posterior": "posteriorEmiss"}, scaleFactor = 0.000001, unit = "MtCO2", sumOrMean = "sum")
        if objBench.isMean:
            resultList = objALP.compute_total({"Prior": "priorEmiss", "Posterior": "posteriorEmiss", "Truth": "truthEmiss"}, scaleFactor = 1.0, unit = "umol/m^2/s", sumOrMean = "mean")
        ax = objALP.show_results(resultList, ax = ax)
        axi += 1

    bcp.cfig.autofmt_xdate()
    bcp.cfig.suptitle("Region time series") 
    bcp.next_page()

    ##################

    axs = bcp.subplots(nrows = Nrow, ncols = 1, sharex = True)
    bcp.cfig.subplots_adjust(hspace = hspace) 
    axi = 0

    for pointName in point_timeSeries_dic:

        print("Drawing -> " + pointName)

        if axi >= Nrow:
            bcp.cfig.autofmt_xdate()
            bcp.cfig.suptitle("Points time series") 
            bcp.next_page()
            axs = bcp.subplots(nrows = Nrow, ncols = 1, sharex = True)
            bcp.cfig.subplots_adjust(hspace = hspace) 
            axi = 0


        ax = axs[axi]
        #tS_dic = point_timeSeries_dic[pointName]
        timeUnit_dic = objBench.myConfig["timeUnit_dic"] #if objBench.isSum else "umol/m^2/s"
        LinePlots_kwargs = {**myConfig[pointName + "_kwargs"], "vrange": point_vrange[timeScale if objBench.isSum else "mean"], "unit": "tCO2/km^2/" + timeUnit_dic[timeScale] if objBench.isSum else "umol/m^2/s", "title": pointName + " location"}
        objALP = AnalysisLinePlots(point_timeSeries_dic[pointName], LinePlots_kwargs = LinePlots_kwargs)
        ax = objALP.plot(ax = ax)
        axi += 1

    bcp.cfig.autofmt_xdate()
    bcp.cfig.suptitle("Points time series")
    bcp.next_page()

    
    bcp.close()

def build_pointRegion_timeSeries_dataset(objBench, pointRegionList, drawStart, drawEnd, timeScale, inherit = False):

    print("Loading point/region time series data for timeScale: " + timeScale)
    pointRegion_timeSeries_datasets = {}
    objAna = objBench.objAna
    areaS = objAna.sampleAreaS
    areaSLON = objAna.sampleLON
    areaSLAT = objAna.sampleLAT
    pointRegionConfig = objBench.myConfig["point_region_kwargs"]
    points_dic = {}
    region_dic = {}
    for pointRegion in pointRegionList:
        config = pointRegionConfig[pointRegion]
        if config["type"].lower() == "point":
            objLoc = LocPoint(staName = pointRegion, staLon = config["location"][0], staLat = config["location"][1])
            points_dic[pointRegion] = objLoc
        elif config["type"].lower() == "region":
            region_dic[pointRegion] = config["shpName"]
        else:
            assert False, config["type"] + " is not one of \"point\" or \"region\"."
    
    ds_priorEmiss_points = objAna.get_priorEmissPoints_timelist(start = drawStart, end = drawEnd, points_dic = points_dic, timeScale = timeScale, interpMethod = "nearest")
    ds_priorSigma_points = objAna.get_priorSigmaPoints_timelist(start = drawStart, end = drawEnd, points_dic = points_dic, timeScale = timeScale, interpMethod = "nearest")
    ds_posteriorEmiss_points = objAna.get_posteriorEmissPoints_timelist(start = drawStart, end = drawEnd, points_dic = points_dic, timeScale = timeScale, interpMethod = "nearest")
    ds_posteriorSigma_points = objAna.get_posteriorSigmaPoints_timelist(start = drawStart, end = drawEnd, points_dic = points_dic, timeScale = timeScale, interpMethod = "nearest")

    for regionName in region_dic:
        objAna.build_region_maskout(regionName = regionName, shpName = region_dic[regionName])
    regionsList = list(region_dic.keys())

    sumOrMean = "sum" if objBench.isSum else "mean"
    ds_priorEmiss_regions = objAna.get_priorEmissRegion_timelist(start = drawStart, end = drawEnd, regions = regionsList, timeScale = timeScale, areaWeighted = False, sumOrMean = sumOrMean)
    ds_priorSigma_regions = objAna.get_priorSigmaRegion_timelist(start = drawStart, end = drawEnd, regions = regionsList, timeScale = timeScale, areaWeighted = False, sumOrMean = sumOrMean)
    ds_posteriorEmiss_regions = objAna.get_posteriorEmissRegion_timelist(start = drawStart, end = drawEnd, regions = regionsList, timeScale = timeScale, areaWeighted = False, sumOrMean = sumOrMean)
    ds_posteriorSigma_regions = objAna.get_posteriorSigmaRegion_timelist(start = drawStart, end = drawEnd, regions = regionsList, timeScale = timeScale, areaWeighted = False, sumOrMean = sumOrMean)

    pointRegion_timeSeries_datasets["ds_priorEmiss_points"] = ds_priorEmiss_points
    pointRegion_timeSeries_datasets["ds_priorSigma_points"] = ds_priorSigma_points
    pointRegion_timeSeries_datasets["ds_posteriorEmiss_points"] = ds_posteriorEmiss_points
    pointRegion_timeSeries_datasets["ds_posteriorSigma_points"] = ds_posteriorSigma_points
    pointRegion_timeSeries_datasets["areaS_points"] = objAna.get_arr_points(areaS, LON = areaSLON, LAT = areaSLAT, points_dic = points_dic, interpMethod = "nearest")

    pointRegion_timeSeries_datasets["ds_priorEmiss_regions"] = ds_priorEmiss_regions
    pointRegion_timeSeries_datasets["ds_priorSigma_regions"] = ds_priorSigma_regions
    pointRegion_timeSeries_datasets["ds_posteriorEmiss_regions"] = ds_posteriorEmiss_regions
    pointRegion_timeSeries_datasets["ds_posteriorSigma_regions"] = ds_posteriorSigma_regions

    objBench.pointRegion_timeSeries_datasets_dic[timeScale] = pointRegion_timeSeries_datasets

    objBench.pointRegion_timeSeries_loaded.append(timeScale)

    if inherit:
        return objAna, points_dic, region_dic, regionsList

def build_pointRegion_timeSeries_dic(objBench, pointRegionList, drawStart, drawEnd, timeScale, inherit = False):

    if timeScale not in objBench.pointRegion_timeSeries_datasets_dic:
        objBench.build_pointRegion_timeSeries_dataset(pointRegionList = pointRegionList, drawStart = drawStart, drawEnd = drawEnd, timeScale = timeScale)
    
    pointRegion_timeSeries_dic = {}
    objAna = objBench.objAna
    pointRegion_timeSeries_datasets = objBench.pointRegion_timeSeries_datasets_dic[timeScale]
    areaS = pointRegion_timeSeries_datasets["areaS_points"]
    areaS = areaS.reshape(len(areaS), 1)
    pointRegionConfig = objBench.myConfig["point_region_kwargs"]

    for pointRegionName in pointRegionList:

        config = pointRegionConfig[pointRegionName]
        typeName = config["type"].lower()

        ds_priorEmiss = pointRegion_timeSeries_datasets["ds_priorEmiss_" + typeName + "s"].copy()
        ds_priorSigma = pointRegion_timeSeries_datasets["ds_priorSigma_" + typeName + "s"].copy()
        ds_posteriorEmiss = pointRegion_timeSeries_datasets["ds_posteriorEmiss_" + typeName + "s"].copy()
        ds_posteriorSigma = pointRegion_timeSeries_datasets["ds_posteriorSigma_" + typeName + "s"].copy()

        if typeName == "point":
            if objBench.isSum:
                ds_priorEmiss.data.values = ds_priorEmiss.data.values / areaS * 1.0e6 
                ds_priorSigma.data.values = ds_priorSigma.data.values / areaS * 1.0e6 
                ds_posteriorEmiss.data.values = ds_posteriorEmiss.data.values / areaS * 1.0e6 
                ds_posteriorSigma.data.values = ds_posteriorSigma.data.values / areaS * 1.0e6 

        ds_priorEmiss = ds_priorEmiss.sel(site = pointRegionName)
        ds_priorSigma = ds_priorSigma.sel(site = pointRegionName)
        ds_posteriorEmiss = ds_posteriorEmiss.sel(site = pointRegionName)
        ds_posteriorSigma = ds_posteriorSigma.sel(site = pointRegionName)

        

        pointRegion_timeSeries_dic[pointRegionName] = {
            "priorEmiss": ds2ATS(ds_line = ds_priorEmiss, ds_spread = ds_priorSigma, label = "Prior"),
            "posteriorEmiss": ds2ATS(ds_line = ds_posteriorEmiss, ds_spread = ds_posteriorSigma, label = "Posterior"),
        }

    if inherit:
        return pointRegion_timeSeries_dic, pointRegionList
    else:
        return pointRegion_timeSeries_dic


