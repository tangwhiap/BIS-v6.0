#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

from . import analyst_independent
from ...utils.netcdf_io import get_ncvar

class Analyst(analyst_independent.Analyst):

    def __init__(self, *args, **kwargs):
        analyst_independent.Analyst.__init__(self, *args, **kwargs)

        self.get_ffePriorFile_name = self.objCase.objExp.get_ffePriorFile_name
        self.get_cftaPriorFile_name = self.objCase.objExp.get_cftaPriorFile_name
        self.get_ffeSigmaPriorFile_name = self.objCase.objExp.get_ffeSigmaPriorFile_name
        self.get_cftaSigmaPriorFile_name = self.objCase.objExp.get_cftaSigmaPriorFile_name

        self.get_ffeHourlyFile_name = self.objCase.objExp.get_ffeHourlyFile_name
        self.get_sffeHourlyFile_name = self.objCase.objExp.get_sffeHourlyFile_name
        self.get_cftaHourlyFile_name = self.objCase.objExp.get_cftaHourlyFile_name
        self.get_scftaHourlyFile_name = self.objCase.objExp.get_scftaHourlyFile_name

        self.get_ffeDailyFile_name = self.objCase.objExp.get_ffeDailyFile_name
        self.get_sffeDailyFile_name = self.objCase.objExp.get_sffeDailyFile_name
        self.get_cftaDailyFile_name = self.objCase.objExp.get_cftaDailyFile_name
        self.get_scftaDailyFile_name = self.objCase.objExp.get_scftaDailyFile_name

        self.get_ffeWeeklyFile_name = self.objCase.objExp.get_ffeWeeklyFile_name
        self.get_sffeWeeklyFile_name = self.objCase.objExp.get_sffeWeeklyFile_name
        self.get_cftaWeeklyFile_name = self.objCase.objExp.get_cftaWeeklyFile_name
        self.get_scftaWeeklyFile_name = self.objCase.objExp.get_scftaWeeklyFile_name


        self.get_ffeMonthlyFile_name = self.objCase.objExp.get_ffeMonthlyFile_name
        self.get_sffeMonthlyFile_name = self.objCase.objExp.get_sffeMonthlyFile_name
        self.get_cftaMonthlyFile_name = self.objCase.objExp.get_cftaMonthlyFile_name
        self.get_scftaMonthlyFile_name = self.objCase.objExp.get_scftaMonthlyFile_name

        self.get_ffeAllFile_name = self.objCase.objExp.get_ffeAllFile_name
        self.get_sffeAllFile_name = self.objCase.objExp.get_sffeAllFile_name
        self.get_cftaAllFile_name = self.objCase.objExp.get_cftaAllFile_name
        self.get_scftaAllFile_name = self.objCase.objExp.get_scftaAllFile_name

    def get_ffe_emiss(self, time, sectorName = None, timeScale = "orig"):

        if sectorName is None:
            sectorName = self.sampleSector

        if timeScale.lower() == "orig":
            fileName = self.get_ffePriorFile_name(time, sectorName)

        elif timeScale.lower() == "hourly":
            fileName = self.get_ffeHourlyFile_name(time, sectorName)

        elif timeScale.lower() == "daily":
            fileName = self.get_ffeDailyFile_name(time, sectorName)

        elif timeScale.lower() == "weekly":
            fileName = self.get_ffeWeeklyFile_name(time, sectorName)

        elif timeScale.lower() == "monthly":
            fileName = self.get_ffeMonthlyFile_name(time, sectorName)

        elif timeScale.lower() == "all":
            fileName = self.get_ffeAllFile_name(time, sectorName)

        else:
            assert False

        emiss, LON, LAT = get_ncvar(fileName, ["emiss", "LON", "LAT"])

        return emiss, LON, LAT

    def get_cfta_emiss(self, time, sectorName = None, timeScale = "orig"):

        if sectorName is None:
            sectorName = self.sampleSector

        if timeScale.lower() == "orig":
            fileName = self.get_cftaPriorFile_name(time, sectorName)

        elif timeScale.lower() == "hourly":
            fileName = self.get_cftaHourlyFile_name(time, sectorName)

        elif timeScale.lower() == "daily":
            fileName = self.get_cftaDailyFile_name(time, sectorName)

        elif timeScale.lower() == "weekly":
            fileName = self.get_cftaWeeklyFile_name(time, sectorName)

        elif timeScale.lower() == "monthly":
            fileName = self.get_cftaMonthlyFile_name(time, sectorName)

        elif timeScale.lower() == "all":
            fileName = self.get_cftaAllFile_name(time, sectorName)

        else:
            assert False

        emiss, LON, LAT = get_ncvar(fileName, ["emiss", "LON", "LAT"])

        return emiss, LON, LAT


    def get_ffe_sigma(self, time, sectorName = None, timeScale = "orig"):

        if sectorName is None:
            sectorName = self.sampleSector

        if timeScale.lower() == "orig":
            fileName = self.get_ffeSigmaPriorFile_name(time, sectorName)

        elif timeScale.lower() == "hourly":
            fileName = self.get_sffeHourlyFile_name(time, sectorName)

        elif timeScale.lower() == "daily":
            fileName = self.get_sffeDailyFile_name(time, sectorName)

        elif timeScale.lower() == "weekly":
            fileName = self.get_sffeWeeklyFile_name(time, sectorName)

        elif timeScale.lower() == "monthly":
            fileName = self.get_sffeMonthlyFile_name(time, sectorName)

        elif timeScale.lower() == "all":
            fileName = self.get_sffeAllFile_name(time, sectorName)

        else:
            assert False

        emiss, LON, LAT = get_ncvar(fileName, ["emiss", "LON", "LAT"])

        return emiss, LON, LAT

    def get_cfta_sigma(self, time, sectorName = None, timeScale = "orig"):

        if sectorName is None:
            sectorName = self.sampleSector

        if timeScale.lower() == "orig":
            fileName = self.get_cftaSigmaPriorFile_name(time, sectorName)

        elif timeScale.lower() == "hourly":
            fileName = self.get_scftaHourlyFile_name(time, sectorName)

        elif timeScale.lower() == "daily":
            fileName = self.get_scftaDailyFile_name(time, sectorName)

        elif timeScale.lower() == "weekly":
            fileName = self.get_scftaWeeklyFile_name(time, sectorName)

        elif timeScale.lower() == "monthly":
            fileName = self.get_scftaMonthlyFile_name(time, sectorName)

        elif timeScale.lower() == "all":
            fileName = self.get_scftaAllFile_name(time, sectorName)

        else:
            assert False

        emiss, LON, LAT = get_ncvar(fileName, ["emiss", "LON", "LAT"])

        return emiss, LON, LAT

    def get_ffeEmiss_points(self, time, sectorName = None, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        emiss, LON, LAT = self.get_ffe_emiss(time, sectorName = sectorName, timeScale = timeScale)
        if points_dic is None:
            points_dic = self.obsPoints
        return self.get_arr_points(emiss, LON, LAT, points_dic, interpMethod = interpMethod)

    def get_ffeSigma_points(self, time, sectorName = None, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        emiss, LON, LAT = self.get_ffe_sigma(time, sectorName = sectorName, timeScale = timeScale)
        if points_dic is None:
            points_dic = self.obsPoints
        return self.get_arr_points(emiss, LON, LAT, points_dic, interpMethod = interpMethod)

    def get_cftaEmiss_points(self, time, sectorName = None, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        emiss, LON, LAT = self.get_cfta_emiss(time, sectorName = sectorName, timeScale = timeScale)
        if points_dic is None:
            points_dic = self.obsPoints
        return self.get_arr_points(emiss, LON, LAT, points_dic, interpMethod = interpMethod)

    def get_cftaSigma_points(self, time, sectorName = None, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        emiss, LON, LAT = self.get_cfta_sigma(time, sectorName = sectorName, timeScale = timeScale)
        if points_dic is None:
            points_dic = self.obsPoints
        return self.get_arr_points(emiss, LON, LAT, points_dic, interpMethod = interpMethod)

    def get_ffeEmiss_region(self, time, regions, sectorName = None, timeScale = "orig",  areaWeighted = True):
        #emiss, LON, LAT = self.get_ffe_emiss(time, sectorName = sectorName, timeScale = timeScale)
        regionsValue = []
        for region in regions:
            objRegion = self.regionsMaskout[region]
            emissValue = self.spatial_dimension_total(func = self.get_ffe_emiss, maskout = objRegion.mask_array, time = time, sectorName = sectorName, timeScale = timeScale,  areaWeighted = areaWeighted)
            regionsValue.append(emissValue)
        return regionsValue

    def get_ffeSigma_region(self, time, regions, sectorName = None, timeScale = "orig",  areaWeighted = True):
        #emiss, LON, LAT = self.get_ffe_emiss(time, sectorName = sectorName, timeScale = timeScale)
        regionsValue = []
        for region in regions:
            objRegion = self.regionsMaskout[region]
            emissValue = self.spatial_dimension_total2(func = self.get_ffe_sigma, maskout = objRegion.mask_array, time = time, sectorName = sectorName, timeScale = timeScale,  areaWeighted = areaWeighted)
            regionsValue.append(emissValue)
        return regionsValue

    def get_cftaEmiss_region(self, time, regions, sectorName = None, timeScale = "orig",  areaWeighted = True):
        #emiss, LON, LAT = self.get_cfta_emiss(time, sectorName = sectorName, timeScale = timeScale)
        regionsValue = []
        for region in regions:
            objRegion = self.regionsMaskout[region]
            emissValue = self.spatial_dimension_total(func = self.get_cfta_emiss, maskout = objRegion.mask_array, time = time, sectorName = sectorName, timeScale = timeScale,  areaWeighted = areaWeighted)
            regionsValue.append(emissValue)
        return regionsValue

    def get_cftaSigma_region(self, time, regions, sectorName = None, timeScale = "orig",  areaWeighted = True):
        #emiss, LON, LAT = self.get_cfta_emiss(time, sectorName = sectorName, timeScale = timeScale)
        regionsValue = []
        for region in regions:
            objRegion = self.regionsMaskout[region]
            emissValue = self.spatial_dimension_total2(func = self.get_cfta_sigma, maskout = objRegion.mask_array, time = time, sectorName = sectorName, timeScale = timeScale,  areaWeighted = areaWeighted)
            regionsValue.append(emissValue)
        return regionsValue



    def get_ffeEmissPoints_timelist(self, start, end, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        if points_dic is None:
            site_nameList = None
        else:
            site_nameList = list(points_dic.keys())
        getFun = self.get_ffeEmiss_points
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = site_nameList, getFun = getFun, getFun_kwargs = {"timeScale": timeScale, "interpMethod": interpMethod, "points_dic": points_dic})

    def get_ffeSigmaPoints_timelist(self, start, end, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        if points_dic is None:
            site_nameList = None
        else:
            site_nameList = list(points_dic.keys())
        getFun = self.get_ffeSigma_points
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = site_nameList, getFun = getFun, getFun_kwargs = {"timeScale": timeScale, "interpMethod": interpMethod, "points_dic": points_dic})

    def get_cftaEmissPoints_timelist(self, start, end, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        if points_dic is None:
            site_nameList = None
        else:
            site_nameList = list(points_dic.keys())
        getFun = self.get_cftaEmiss_points
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = site_nameList, getFun = getFun, getFun_kwargs = {"timeScale": timeScale, "interpMethod": interpMethod, "points_dic": points_dic})

    def get_cftaSigmaPoints_timelist(self, start, end, points_dic = None, timeScale = "orig", interpMethod = "belinear"):
        if points_dic is None:
            site_nameList = None
        else:
            site_nameList = list(points_dic.keys())
        getFun = self.get_cftaSigma_points
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = site_nameList, getFun = getFun, getFun_kwargs = {"timeScale": timeScale, "interpMethod": interpMethod, "points_dic": points_dic})

    def get_ffeEmissRegion_timelist(self, start, end, regions, timeScale = "orig", areaWeighted = True):
        getFun = self.get_ffeEmiss_region
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = regions, getFun = getFun, getFun_kwargs = {"regions": regions, "timeScale": timeScale, "areaWeighted": areaWeighted})
     
    def get_ffeSigmaRegion_timelist(self, start, end, regions, timeScale = "orig", areaWeighted = True):
        getFun = self.get_ffeSigma_region
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = regions, getFun = getFun, getFun_kwargs = {"regions": regions, "timeScale": timeScale, "areaWeighted": areaWeighted})

    def get_cftaEmissRegion_timelist(self, start, end, regions, timeScale = "orig", areaWeighted = True):
        getFun = self.get_cftaEmiss_region
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = regions, getFun = getFun, getFun_kwargs = {"regions": regions, "timeScale": timeScale, "areaWeighted": areaWeighted})
     
    def get_cftaSigmaRegion_timelist(self, start, end, regions, timeScale = "orig", areaWeighted = True):
        getFun = self.get_cftaSigma_region
        timeList_type = timeScale
        return self.get_recepTimelist(start = start, end = end, timeList_type = timeList_type, site_nameList = regions, getFun = getFun, getFun_kwargs = {"regions": regions, "timeScale": timeScale, "areaWeighted": areaWeighted})
     