#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...

from ..benchplots_independent.main_class import BenchPlots_independent
from ...benchplots_base import benchplots_configure  as benchConfig
from ...Analyst.analyst_independent_priorFnet import Analyst
from ...benchplots_base.draw_object_class import ds2ATS

class BenchPlots_independent_priorFnet(BenchPlots_independent):

    def __init__(self, *args, **kwargs):

        objCase = kwargs["objCase"]
        objAna = Analyst(objCase)

        BenchPlots_independent.__init__(self, *args, objAna = objAna, **kwargs)
        #self.myConfig = {**self.myConfig, **benchConfig.expREAL_independent_priorFnet_kwargs}
        expKwargs = benchConfig.expREAL_independent_priorFnet_kwargs
        
        fieldPlots_vminmax_append = expKwargs["fieldPlots_vminmax"]
        
        for timeScale in fieldPlots_vminmax_append:
            self.myConfig["fieldPlots_vminmax"][timeScale] = {

                **self.myConfig["fieldPlots_vminmax"][timeScale],
                **fieldPlots_vminmax_append[timeScale],

            }
        
        fieldPlots_kwargs_append = expKwargs["fieldPlots_kwargs"]

        self.myConfig["fieldPlots_kwargs"] = {

            **self.myConfig["fieldPlots_kwargs"],
            **fieldPlots_kwargs_append,

        }
        
        movie_vminmax_append = expKwargs["movie_vminmax"]
        self.myConfig["movie_vminmax"] = {

            **self.myConfig["movie_vminmax"],
            **movie_vminmax_append,

        }

        movie_kwargs_append = expKwargs["movie_kwargs"]        
        self.myConfig["movie_kwargs"] = {

            **self.myConfig["movie_kwargs"],
            **movie_kwargs_append,

        }
        ffe_cfta_Plots_kwargs = expKwargs["ffe_cfta_Plots_kwargs"]

        for pointRegion in self.myConfig["pointRegionList"]:
            self.myConfig["point_region_timeSeries_kwargs"][pointRegion + "_kwargs"] = {

                **self.myConfig["point_region_timeSeries_kwargs"][pointRegion + "_kwargs"],
                **ffe_cfta_Plots_kwargs,

            }

    def build_movie_field_dic(self, time): 

        objAna = self.objAna
        field_dic = BenchPlots_independent.build_movie_field_dic(self, time)
        field_dic["ffeEmiss"] = objAna.get_ffe_emiss(time)
        field_dic["cftaEmiss"] = objAna.get_cfta_emiss(time)

        return field_dic

    def build_field_dic(self, time, timeScale):

        objAna = self.objAna
        field_dic = BenchPlots_independent.build_field_dic(self, time, timeScale)

        ffeEmiss, LON, LAT = objAna.get_ffe_emiss(time, timeScale = timeScale)
        cftaEmiss, _, _ = objAna.get_cfta_emiss(time, timeScale = timeScale)

        #ffeSigma, _, _ = objAna.get_ffe_sigma(time, timeScale = timeScale)
        #cftaSigma, _, _ = objAna.get_cfta_sigma(time, timeScale = timeScale)

        if self.isSum:
            areaS = objAna.sampleAreaS
            
            ffeEmiss *= (0.001 * areaS / 1.0e6)
            cftaEmiss *= (0.001 * areaS / 1.0e6)

        field_dic["ffeEmiss"] = ffeEmiss, LON, LAT
        field_dic["cftaEmiss"] = cftaEmiss, LON, LAT

        return field_dic

    def build_pointRegion_timeSeries_dataset(self, pointRegionList, drawStart, drawEnd, timeScale):

        objAna, points_dic, region_dic, regionsList = BenchPlots_independent.build_pointRegion_timeSeries_dataset(self, pointRegionList, drawStart, drawEnd, timeScale, inherit = True)
        prTs_datasets_dic = self.pointRegion_timeSeries_datasets_dic[timeScale]

        ds_ffeEmiss_points = objAna.get_ffeEmissPoints_timelist(start = drawStart, end = drawEnd, points_dic = points_dic, timeScale = timeScale, interpMethod = "nearest")
        ds_ffeSigma_points = objAna.get_ffeSigmaPoints_timelist(start = drawStart, end = drawEnd, points_dic = points_dic, timeScale = timeScale, interpMethod = "nearest")
        ds_cftaEmiss_points = objAna.get_cftaEmissPoints_timelist(start = drawStart, end = drawEnd, points_dic = points_dic, timeScale = timeScale, interpMethod = "nearest")
        ds_cftaSigma_points = objAna.get_cftaSigmaPoints_timelist(start = drawStart, end = drawEnd, points_dic = points_dic, timeScale = timeScale, interpMethod = "nearest")

        sumOrMean = "sum" if self.isSum else "mean"
        ds_ffeEmiss_regions = objAna.get_ffeEmissRegion_timelist(start = drawStart, end = drawEnd, regions = regionsList, timeScale = timeScale, areaWeighted = False, sumOrMean = sumOrMean)
        ds_ffeSigma_regions = objAna.get_ffeSigmaRegion_timelist(start = drawStart, end = drawEnd, regions = regionsList, timeScale = timeScale, areaWeighted = False, sumOrMean = sumOrMean)
        ds_cftaEmiss_regions = objAna.get_cftaEmissRegion_timelist(start = drawStart, end = drawEnd, regions = regionsList, timeScale = timeScale, areaWeighted = False, sumOrMean = sumOrMean)
        ds_cftaSigma_regions = objAna.get_cftaSigmaRegion_timelist(start = drawStart, end = drawEnd, regions = regionsList, timeScale = timeScale, areaWeighted = False, sumOrMean = sumOrMean)

        prTs_datasets_dic["ds_ffeEmiss_points"] = ds_ffeEmiss_points
        prTs_datasets_dic["ds_ffeSigma_points"] = ds_ffeSigma_points
        prTs_datasets_dic["ds_cftaEmiss_points"] = ds_cftaEmiss_points
        prTs_datasets_dic["ds_cftaSigma_points"] = ds_cftaSigma_points

        prTs_datasets_dic["ds_ffeEmiss_regions"] = ds_ffeEmiss_regions
        prTs_datasets_dic["ds_ffeSigma_regions"] = ds_ffeSigma_regions
        prTs_datasets_dic["ds_cftaEmiss_regions"] = ds_cftaEmiss_regions
        prTs_datasets_dic["ds_cftaSigma_regions"] = ds_cftaSigma_regions

        self.pointRegion_timeSeries_datasets_dic[timeScale] = prTs_datasets_dic

    def build_pointRegion_timeSeries_dic(self, pointRegionList, drawStart, drawEnd, timeScale):

        pointRegion_timeSeries_dic, pointRegionList = BenchPlots_independent.build_pointRegion_timeSeries_dic(self, pointRegionList, drawStart, drawEnd, timeScale, inherit = True)
        pointRegion_timeSeries_datasets = self.pointRegion_timeSeries_datasets_dic[timeScale]
        areaS = pointRegion_timeSeries_datasets["areaS_points"]
        areaS = areaS.reshape(len(areaS), 1)
        pointRegionConfig = self.myConfig["point_region_kwargs"]

        for pointRegionName in pointRegionList:

            config = pointRegionConfig[pointRegionName]
            typeName = config["type"].lower()

            ds_ffeEmiss = pointRegion_timeSeries_datasets["ds_ffeEmiss_" + typeName + "s"].copy()
            ds_ffeSigma = pointRegion_timeSeries_datasets["ds_ffeSigma_" + typeName + "s"].copy()
            ds_cftaEmiss = pointRegion_timeSeries_datasets["ds_cftaEmiss_" + typeName + "s"].copy()
            ds_cftaSigma = pointRegion_timeSeries_datasets["ds_cftaSigma_" + typeName + "s"].copy()

            if typeName == "point":
                if self.isSum:
                    ds_ffeEmiss.data.values = ds_ffeEmiss.data.values / areaS * 1.0e6 
                    ds_ffeSigma.data.values = ds_ffeSigma.data.values / areaS * 1.0e6 
                    ds_cftaEmiss.data.values = ds_cftaEmiss.data.values / areaS * 1.0e6 
                    ds_cftaSigma.data.values = ds_cftaSigma.data.values / areaS * 1.0e6 

            ds_ffeEmiss = ds_ffeEmiss.sel(site = pointRegionName)
            ds_ffeSigma = ds_ffeSigma.sel(site = pointRegionName)
            ds_cftaEmiss = ds_cftaEmiss.sel(site = pointRegionName)
            ds_cftaSigma = ds_cftaSigma.sel(site = pointRegionName)

            pointRegion_timeSeries_dic[pointRegionName] = {
                **pointRegion_timeSeries_dic[pointRegionName], 
                "ffeEmiss": ds2ATS(ds_line = ds_ffeEmiss, ds_spread = ds_ffeSigma, label = "ffe"),
                "cftaEmiss": ds2ATS(ds_line = ds_cftaEmiss, ds_spread = ds_cftaSigma, label = "cfta"),
            }

        return pointRegion_timeSeries_dic





    

