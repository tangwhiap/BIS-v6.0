#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...

from ..Analyst.bis_analyst import Analyst
from ..Analyst.local_point import LocPoint
from .spatial_plot import bis_spatial_plot, draw_fieldPlots_snap, draw_fieldPlots, draw_movie_snap, make_cartoon, build_movie_field_dic, build_field_dic
from .draw_object_class import ds2ATS
from .site_timelist_plot import site_timeSeries, build_site_timeSeries_dataset, build_site_timeSeries_dic, build_sites_observation_dic
from .region_timelist_plot import pointRegion_timeSeries, build_pointRegion_timeSeries_dataset, build_pointRegion_timeSeries_dic
from ...utils.module2dic import module2dic

from . import benchplots_configure as benchConfig

from pdb import set_trace

class BIS_BenchPlots(object):

    def __init__(self, benchDir, objCase = None, objAna = None):
        if objAna is None:
            assert objCase is not None
            self.objAna = Analyst(objCase)
        else:
            self.objAna = objAna

        self.benchDir = benchDir
        self.myConfig = module2dic(benchConfig)

        self.nProc = self.objAna.objCase.objExp.objCore.nProc
        self.site_timeSeries_loaded = False
        self.pointRegion_timeSeries_datasets_dic = {}
        self.pointRegion_timeSeries_loaded = []
        self.isSum = objCase.objExp.myConfig["outputType"].lower() == "sum"
        self.isMean = objCase.objExp.myConfig["outputType"].lower() == "mean"

    def build_movie_field_dic(self, *args, **kwargs):
        return build_movie_field_dic(self, *args, **kwargs)

    def build_field_dic(self, *args, **kwargs):
        return build_field_dic(self, *args, **kwargs)

    def build_site_timeSeries_dic(self, *args, **kwargs):
        return build_site_timeSeries_dic(self, *args, **kwargs)

    #def spatial_plot(self, *args, **kwargs):
    #    return spatial_plot(self, *args, **kwargs)

    def bis_spatial_plot(self, *args, **kwargs):
        return bis_spatial_plot(self, *args, **kwargs)

    def draw_fieldPlots_snap(self, *args, **kwargs):
        return draw_fieldPlots_snap(self, *args, **kwargs)

    def draw_fieldPlots(self, *args, **kwargs):
        return draw_fieldPlots(self, *args, **kwargs)

    def draw_movie_snap(self, *args, **kwargs):
        return draw_movie_snap(self, *args, **kwargs)

    def make_cartoon(self, *args, **kwargs):
        return make_cartoon(self, *args, **kwargs)

    def build_site_timeSeries_dataset(self, *args, **kwargs):
        return build_site_timeSeries_dataset(self, *args, **kwargs)

    def build_site_timeSeries_dic(self, *args, **kwargs):
        return build_site_timeSeries_dic(self, *args, **kwargs)

    def build_sites_observation_dic(self, *args, **kwargs):
        return build_sites_observation_dic(self, *args, **kwargs)

    def build_pointRegion_timeSeries_dataset(self, *args, **kwargs):
        return build_pointRegion_timeSeries_dataset(self, *args, **kwargs)

    def build_pointRegion_timeSeries_dic(self, *args, **kwargs):
        return build_pointRegion_timeSeries_dic(self, *args, **kwargs)

    def site_timeSeries(self, *args, **kwargs):
        return site_timeSeries(self, *args, **kwargs)

    def pointRegion_timeSeries(self, *args, **kwargs):
        return pointRegion_timeSeries(self, *args, **kwargs)

    def run(self, *args, **kwargs):
        main(self, *args, **kwargs)

def main(objBench, drawList = None, spatial_timeScaleList = None, prTs_timeScaleList = None, spatialStart = None, spatialEnd = None, movieStart = None, movieEnd = None, movieDtHrs = None, siteTsStart = None, siteTsEnd = None, prTsStart = None, prTsEnd = None):

    if drawList is None:
        drawList = ["movie", "spatial_plots", "site_timeSeries", "pointRegion_timeSeries"]

    if spatial_timeScaleList is None:
        spatial_timeScaleList = ["daily", "weekly"]

    if prTs_timeScaleList is None:
        prTs_timeScaleList = ["hourly", "daily", "weekly"] 

    expStart = objBench.objAna.objCase.objExp.Start
    expEnd = objBench.objAna.objCase.objExp.End
    expDtHrs = objBench.objAna.objCase.objExp.dtHrs

    if spatialStart is None:
        spatialStart = expStart

    if spatialEnd is None:
        spatialEnd = expEnd

    if movieStart is None:
        movieStart = expStart

    if movieEnd is None:
        movieEnd = expEnd

    if movieDtHrs is None:
        movieDtHrs = expDtHrs

    if siteTsStart is None:
        siteTsStart = expStart

    if siteTsEnd is None:
        siteTsEnd = expEnd

    if prTsStart is None:
        prTsStart = expStart

    if prTsEnd is None:
        prTsEnd = expEnd

    if "movie" in drawList:
        objBench.make_cartoon(movieStart = movieStart, movieEnd = movieEnd, movieDtHrs = movieDtHrs)

    if "spatial_plots" in drawList:
        for timeScale in spatial_timeScaleList:
            objBench.draw_fieldPlots(drawStart = spatialStart, drawEnd = spatialEnd, timeScale = timeScale)

    if "site_timeSeries" in drawList:
        objBench.site_timeSeries(drawStart = siteTsStart, drawEnd = siteTsEnd)

    if "pointRegion_timeSeries" in drawList:
        for timeScale in prTs_timeScaleList:
            objBench.pointRegion_timeSeries(drawStart = prTsStart, drawEnd = prTsEnd, timeScale = timeScale)

