#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 09/2021
#   ...


from ...benchplots_base.main_class import BIS_BenchPlots
from ...benchplots_base import benchplots_configure  as benchConfig
from ...benchplots_base.draw_object_class import ds2ATS

class BenchPlots_OSSE(BIS_BenchPlots):

    def __init__(self, *args, **kwargs):

        BIS_BenchPlots.__init__(self, *args, **kwargs)
        #self.myConfig = {**self.myConfig, **benchConfig.expOSSE_kwargs}
        expKwargs = benchConfig.expOSSE_kwargs

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

        self.myConfig["site_timeSeries_kwargs"]["obs_timeSeries_kwargs"]["vrange"] = expKwargs["site_plots_vrange"]["plot_obs"]
        self.myConfig["site_timeSeries_kwargs"]["plot_concentration_kwargs"]["vrange"] = expKwargs["site_plots_vrange"]["plot_concentration"]
        self.myConfig["site_timeSeries_kwargs"]["plot_Hx_kwargs"]["vrange"] = expKwargs["site_plots_vrange"]["plot_Hx"]
        self.myConfig["site_timeSeries_kwargs"]["plot_bck_kwargs"]["vrange"] = expKwargs["site_plots_vrange"]["plot_bck"]


    def build_movie_field_dic(self, time): 

        objAna = self.objAna
        field_dic = BIS_BenchPlots.build_movie_field_dic(self, time)
        field_dic["truthEmiss"] = objAna.get_truth_emiss(time)
        LON, LAT = field_dic["truthEmiss"][1], field_dic["truthEmiss"][2]
        diff = field_dic["posteriorEmiss"][0] - field_dic["truthEmiss"][0]
        field_dic["posteriorMtruth"] = diff, LON, LAT

        return field_dic

    def build_field_dic(self, time, timeScale):

        objAna = self.objAna
        field_dic = BIS_BenchPlots.build_field_dic(self, time, timeScale)

        truthEmiss, LON, LAT = objAna.get_truth_emiss(time, timeScale = timeScale)
        diffEmiss = field_dic["posteriorEmiss"][0] - truthEmiss

        if self.isSum:
            areaS = objAna.sampleAreaS
            truthEmiss *= (0.001 * areaS / 1.0e6)
            diffEmiss *= (0.001 * areaS / 1.0e6)

        field_dic["truthEmiss"] = truthEmiss, LON, LAT
        field_dic["posteriorMtruth"] = diffEmiss, LON, LAT

        return field_dic

    def build_site_timeSeries_dataset(objBench, drawStart, drawEnd):
        BIS_BenchPlots.build_site_timeSeries_dataset(objBench, drawStart, drawEnd)
        objBench.site_timeSeries_dataset["ds_truthEmiss"] = objBench.objAna.get_truthEmissPoints_timelist(drawStart, drawEnd)

    def build_site_timeSeries_dic(objBench, siteName, drawStart, drawEnd):
        site_timeSeries_dic = BIS_BenchPlots.build_site_timeSeries_dic(objBench, siteName, drawStart, drawEnd)
        ds_truthEmiss = objBench.site_timeSeries_dataset["ds_truthEmiss"].sel(site = siteName)
        site_timeSeries_dic["plot_point_emiss"]["truthEmiss"] = ds2ATS(ds_line = ds_truthEmiss, label = "Truth")
        return site_timeSeries_dic

    def build_pointRegion_timeSeries_dataset(self, pointRegionList, drawStart, drawEnd, timeScale, inherit = False):

        objAna, points_dic, region_dic, regionsList = BIS_BenchPlots.build_pointRegion_timeSeries_dataset(self, pointRegionList, drawStart, drawEnd, timeScale, inherit = True)
        prTs_datasets_dic = self.pointRegion_timeSeries_datasets_dic[timeScale]

        ds_truthEmiss_points = objAna.get_truthEmissPoints_timelist(start = drawStart, end = drawEnd, points_dic = points_dic, timeScale = timeScale, interpMethod = "nearest")

        sumOrMean = "sum" if self.isSum else "mean"
        ds_truthEmiss_regions = objAna.get_truthEmissRegion_timelist(start = drawStart, end = drawEnd, regions = regionsList, timeScale = timeScale, areaWeighted = False, sumOrMean = sumOrMean)

        prTs_datasets_dic["ds_truthEmiss_points"] = ds_truthEmiss_points

        prTs_datasets_dic["ds_truthEmiss_regions"] = ds_truthEmiss_regions

        self.pointRegion_timeSeries_datasets_dic[timeScale] = prTs_datasets_dic
        
        if inherit:
            return objAna, points_dic, region_dic, regionsList


    def build_pointRegion_timeSeries_dic(self, pointRegionList, drawStart, drawEnd, timeScale, inherit = False):

        pointRegion_timeSeries_dic, pointRegionList = BIS_BenchPlots.build_pointRegion_timeSeries_dic(self, pointRegionList, drawStart, drawEnd, timeScale, inherit = True)
        pointRegion_timeSeries_datasets = self.pointRegion_timeSeries_datasets_dic[timeScale]
        areaS = pointRegion_timeSeries_datasets["areaS_points"]
        areaS = areaS.reshape(len(areaS), 1)
        pointRegionConfig = self.myConfig["point_region_kwargs"]

        for pointRegionName in pointRegionList:

            config = pointRegionConfig[pointRegionName]
            typeName = config["type"].lower()

            ds_truthEmiss = pointRegion_timeSeries_datasets["ds_truthEmiss_" + typeName + "s"].copy()

            if typeName == "point":
                if self.isSum:
                    ds_truthEmiss.data.values = ds_truthEmiss.data.values / areaS * 1.0e6 

            ds_truthEmiss = ds_truthEmiss.sel(site = pointRegionName)

            pointRegion_timeSeries_dic[pointRegionName] = {
                **pointRegion_timeSeries_dic[pointRegionName], 
                "truthEmiss": ds2ATS(ds_line = ds_truthEmiss, label = "Truth"),
            }

        if inherit:
            return pointRegion_timeSeries_dic, pointRegionList
        else:
            return pointRegion_timeSeries_dic

