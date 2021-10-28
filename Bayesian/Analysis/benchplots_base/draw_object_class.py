#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader

def alpha_vary_cmap(colormap):
    cmap = plt.get_cmap(colormap)
    my_cmap = cmap(np.arange(cmap.N))
    my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
    my_cmap = ListedColormap(my_cmap)
    return my_cmap

def spatial_plot(Arr, LON, LAT, fig, gs, shpDir, vmin = None, vmax = None, cmap = "jet", cmap_va = False, draw_gridline = True, draw_lonlat_label = True, pcolormesh_kwargs = {}):

    def alpha_vary_cmap(colormap):
        cmap = plt.get_cmap(colormap)
        my_cmap = cmap(np.arange(cmap.N))
        my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
        my_cmap = ListedColormap(my_cmap)
        return my_cmap

    assert Arr.shape == LON.shape and LON.shape == LAT.shape

    proj = ccrs.PlateCarree()
    geo_ax = fig.add_subplot(gs, projection = proj)
    lon_s = LON.min(); lon_e = LON.max() 
    lat_s = LAT.min(); lat_e = LAT.max()

    # Atloc setting
    #lon_s = 115.2
    #lon_e = 117.6
    #lat_s = 39.4
    #lat_e = 41.1
    ################

    geo_ax.set_extent([lon_s, lon_e, lat_s, lat_e], crs = proj)

    if draw_gridline:
        gl = geo_ax.gridlines(crs = proj, linestyle = "--", alpha = 0.5, draw_labels = True)
        gl.top_labels = False
        gl.right_labels = False
        if not draw_lonlat_label:
            gl.left_lables = False
            gl.bottom_lables = False

    geo_ax.add_geometries(Reader(shpDir).geometries(), proj, facecolor = "none", edgecolor = "k", linewidth = 1)

    if cmap_va:
        cmap = alpha_vary_cmap(cmap)
    cs = geo_ax.pcolormesh(LON, LAT, Arr, cmap = cmap, vmin = vmin, vmax = vmax, **pcolormesh_kwargs)
    return geo_ax, cs

def stations_plot(lonList, latList, fig, gs, shpDir, dataList = None, spatialRange = None, vmin = None, vmax = None, cmap = "jet", draw_gridline = True, draw_lonlat_label = True, scatter_kwargs = {}):

    assert len(lonList) == len(latList)
    if spatialRange is None:
        spatialRange = {
            "lon_s": 109, "lon_e": 122,
            "lat_s": 35, "lat_e": 43.2,
        }
    proj = ccrs.PlateCarree()
    geo_ax = fig.add_subplot(gs, projection = proj)
    #lon_s = LON.min(); lon_e = LON.max() 
    #lat_s = LAT.min(); lat_e = LAT.max()

    geo_ax.set_extent([spatialRange["lon_s"], spatialRange["lon_e"], spatialRange["lat_s"], spatialRange["lat_e"]], crs = proj)

    if draw_gridline:
        gl = geo_ax.gridlines(crs = proj, linestyle = "--", alpha = 0.5, draw_labels = True)
        gl.top_labels = False
        gl.right_labels = False
        if not draw_lonlat_label:
            gl.left_lables = False
            gl.bottom_lables = False

    geo_ax.add_geometries(Reader(shpDir).geometries(), proj, facecolor = "none", edgecolor = "k", linewidth = 1)

    #if cmap_va:
    #    cmap = alpha_vary_cmap(cmap)
    #cs = geo_ax.pcolormesh(LON, LAT, Arr, cmap = cmap, vmin = vmin, vmax = vmax, **pcolormesh_kwargs)
    if dataList is None:
        cs = geo_ax.scatter(lonList, latList, **scatter_kwargs)
    else:
        assert len(dataList) == len(lonList)
        scatter_kwargs["c"] = dataList
        cs = geo_ax.scatter(lonList, latList, vmin = vmin, vmax = vmax, cmap = cmap, **scatter_kwargs)
    return geo_ax, cs

class AnalysisField(object):

    def __init__(self, array, LON, LAT, title = "array", unit = "unitless", vmin = None, vmax = None):#, gsLoc = [0, 0]):
        self.array = array
        self.LON = LON
        self.LAT = LAT
        self.vmin = vmin
        self.vmax = vmax
        #self.gsLoc = gsLoc
        self.title = title
        self.unit = unit
        #assert len(self.gsLoc) == 2

    def spatial_plot(self, fig, gs, shpDir, spatial_kwargs = {}, pcolormesh_kwargs = {}, cbar_kwargs = {}):
        """
        spatial_kwargs = {
            cmap = "jet",
            cmap_va = False, 
            draw_gridline = True, 
            draw_lonlat_label = True,
        }
        """
        assert self.array.shape == self.LON.shape and self.LON.shape == self.LAT.shape
        #gs = gses[self.gsLoc[0], self.gsLoc[1]]
        geo_ax, cs = spatial_plot(self.array, self.LON, self.LAT, fig = fig, gs = gs, shpDir = shpDir, vmin = self.vmin, vmax = self.vmax, pcolormesh_kwargs = pcolormesh_kwargs, **spatial_kwargs)
        geo_ax.set_title(self.title)
        cbar = fig.colorbar(cs, ax = geo_ax, **cbar_kwargs)
        cbar.set_label(self.unit)
        return fig, geo_ax, cs, cbar

class AnalysisStations(object):

    def __init__(self, lonList, latList, dataList = None, title = "array", unit = None, vmin = None, vmax = None):#, gsLoc = [0, 0]):
        self.dataList = dataList
        self.lonList = lonList 
        self.latList = latList
        self.vmin = vmin
        self.vmax = vmax
        #self.gsLoc = gsLoc
        self.title = title
        self.unit = unit
        #assert len(self.gsLoc) == 2

    def stations_plot(self, fig, gs, shpDir, spatialRange = None, stationsPlot_kwargs = {}, scatter_kwargs = {}, cbar_kwargs = {}):
        geo_ax, cs = stations_plot(lonList = self.lonList, latList = self.latList, fig = fig, gs = gs, shpDir = shpDir, dataList = self.dataList, spatialRange = spatialRange, vmin = self.vmin, vmax = self.vmax, scatter_kwargs = scatter_kwargs, **stationsPlot_kwargs)
        geo_ax.set_title(self.title)
        if self.dataList is not None:
            cbar = fig.colorbar(cs, ax = geo_ax, **cbar_kwargs)
            cbar.set_label(self.unit)
        else:
            cbar = None
        return fig, geo_ax, cs, cbar



class AnalysisTimeSeries(object):

    def __init__(self, time, data, label, spread = None):
        self.time = time
        self.data = data
        self.label = label
        self.spread = spread
        assert len(self.time) == len(self.data)
        self.hasSpread = self.spread is not None
        if self.hasSpread:
            assert len(self.time) == len(self.spread)

    def plot(self, ax, plot_kwargs = {}, spread_kwargs = {}):
        line = ax.plot(self.time, self.data, label = self.label, **plot_kwargs)
        if self.hasSpread:
            fb = ax.fill_between(self.time, self.data - self.spread, self.data + self.spread, **spread_kwargs)
        else:
            fb = None
        return line, fb

    def compute_sum(self):
        total = np.nansum(self.data)
        if self.hasSpread:
            total_spread = np.sqrt(np.nansum(self.spread ** 2))
        else:
            total_spread = None
        return total, total_spread

    def compute_mean(self):
        mean_data = np.nanmean(self.data)
        if self.hasSpread:
            mean_spread = np.sqrt(np.nansum(self.spread **2)) / np.sum(~np.isnan(self.spread))
        else:
            mean_spread = None
        return mean_data, mean_spread

class AnalysisLinePlots(object):

    def __init__(self, timeSeries_dic, LinePlots_kwargs = {}):
        self.timeSeries_dic = timeSeries_dic
        self.LinePlots_kwargs = LinePlots_kwargs
        for tsName in self.timeSeries_dic:
            if tsName + "_kwargs" not in self.LinePlots_kwargs:
                self.LinePlots_kwargs[tsName + "_kwargs"] = {"plot_kwargs": {}, "spread_kwargs": {}}

        if "vrange" not in self.LinePlots_kwargs:
            self.LinePlots_kwargs["vrange"] = (None, None)

        if "grid_kwargs" not in self.LinePlots_kwargs:
            self.LinePlots_kwargs["grid_kwargs"] = {}
        
        if "legend_kwargs" not in self.LinePlots_kwargs:
            self.LinePlots_kwargs["legend_kwargs"] = {}

        if "unit" not in self.LinePlots_kwargs:
            self.LinePlots_kwargs["unit"] = ""


    def plot(self, ax):
        for tsName in self.timeSeries_dic:
            self.timeSeries_dic[tsName].plot(ax, **self.LinePlots_kwargs[tsName + "_kwargs"])

        if "yLine" in self.LinePlots_kwargs:
            yLines = self.LinePlots_kwargs["yLine"]
            for yLine in yLines:
                yValue, yLine_kwargs = yLine
                ax = self.plot_yline(ax = ax, yValue = yValue, **yLine_kwargs)

        ax.set_ylim(self.LinePlots_kwargs["vrange"])
        ax.set_ylabel(self.LinePlots_kwargs["unit"])
        ax.legend(**self.LinePlots_kwargs["legend_kwargs"])
        ax.grid(**self.LinePlots_kwargs["grid_kwargs"])

        if "title" in self.LinePlots_kwargs:
            ax.set_title(self.LinePlots_kwargs["title"])

        return ax

    def compute_total(self, timeSeries_nameDic, scaleFactor = 1, unit = "", sumOrMean = "sum"):
        sumOrMean = sumOrMean.lower()
        assert sumOrMean in ["sum", "mean"]
        resultList = []
        for labelName in timeSeries_nameDic:
            tsName = timeSeries_nameDic[labelName]
            if sumOrMean == "sum":
                total, total_spread = self.timeSeries_dic[tsName].compute_sum()
            elif sumOrMean == "mean":
                total, total_spread = self.timeSeries_dic[tsName].compute_mean()
            hasSpread = self.timeSeries_dic[tsName].hasSpread
            total = total * scaleFactor
            if hasSpread:
                total_spread = total_spread * scaleFactor
            result = labelName + ": %.4f" % (total)
            if hasSpread:
                result = result + " ± " + "%.4f" % (total_spread)
            result = result + " " + unit
            resultList.append(result)
        return resultList

    def show_results(self, resultList, ax):
        text_sx = 0.05
        text_sy = 0.95
        text_yint = 0.1
        text_ix = text_sx
        text_iy = text_sy
        for result in resultList:
            ax.text(text_ix, text_iy, result, transform = ax.transAxes, horizontalalignment = "left", verticalalignment = "top")
            text_iy = text_iy - text_yint
        return ax







        

            
    
    def plot_yline(self, ax, yValue, **kwargs):
        t1, t2 = self.get_tlim()
        ax.plot([t1, t2], [yValue, yValue], **kwargs)
        #ax.set_xlim([x1, x2])
        return ax
    
    def get_tlim(self):
        t1List = []
        t2List = []
        for tsName in self.timeSeries_dic:
            objTS = self.timeSeries_dic[tsName]
            t1List.append(np.array(objTS.time).min())
            t2List.append(np.array(objTS.time).max())
        t1 = np.array(t1List).min()
        t2 = np.array(t2List).max()
        return t1, t2



def ds2ATS(ds_line, ds_spread = None, label = ""):
    time, data = ds_line["time"].values, ds_line["data"].values
    if ds_spread is not None:
        spread = ds_spread["data"].values
    else:
        spread = None
    return AnalysisTimeSeries(time = time, data = data, label = label, spread = spread)



