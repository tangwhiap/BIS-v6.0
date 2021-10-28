#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 08/2021
#   ...


###########################
#      Basic setting      #
###########################

#--- Shape file directory and name ---#
#shpDir = "/home/tangwh/datasets/china_shp/cnhimap.shp"
shpDir = "/home/tangwh/datasets/china_shp/new/data/CHN_adm1.shp"

#--- define time unit ---#
timeUnit_dic = {"hourly": "hour", "daily": "day", "weekly": "week", "monthly": "month"}

#--- spatial movie yrange ---#
movie_vminmax = {
    "priorEmiss": (0, 80), # prior flux yrange. (unit: umol/m^2/s)
    "posteriorEmiss": (0, 80), # posterior flux yrange. (unit: umol/m^2/s)
    "reducedUncertainty": (0, 1000), # reduced uncertainty yrange. (unit: (umol/m^2/s)^2)
    "posteriorMprior": (-50, 50), # difference between posterior and prior flux. (unit: umol/m^2/s)
}

#--- movie keywords ---#
movie_kwargs = {

    "figsize": (10, 8), # movie figure size.
    "movieName": "spatial_movie", # mp4 file name.
    "dpi": 200, # dpi value.
    "pcolormesh_kwargs": {},  # pyplot.pcolormesh keywords
    "cbar_kwargs": {}, # pyplot.colorbar keywords.
    "spatial_kwargs": # keywords for spatial plotting function.
    {
        "cmap": "coolwarm", # colormap.
        "cmap_va": False, 
        "draw_gridline": True,  # drawing grid lines.
        "draw_lonlat_label": True, # drawing longitude & latitude label.
    },
    "priorEmiss_kwargs": # other keywords for prior flux plotting.
    {
        "unit": "umol/m^2/s",  # unit label.
        "title": "prior", # title.
    },
    "posteriorEmiss_kwargs": # other keywords for posterior flux plotting.
    {
        "unit": "umol/m^2/s", # unit label.
        "title": "posterior", # title.
    },
    "reducedUncertainty_kwargs": # other keywords for reduced uncertainty plotting.
    {
        "unit": "(umol/m^2/s)^2", # unit label.
        "title": "reduced uncertainty", # title.
    },
    "posteriorMprior_kwargs": # other keywords for posterior - prior plotting.
    {
        "unit": "umol/m^2/s", # unit label.
        "title": "posterior - prior", # title.
    }, 
}

#--- 2D field plots yrange settings for each time scale  ---#
fieldPlots_vminmax = {
    "mean":
    {
        "priorEmiss": (0, 80), # prior flux yrange. (unit: umol/m^2/s)
        "posteriorEmiss": (0,80), # posterior flux yrange. (unit: umol/m^2/s)
        "reducedUncertainty": (0, 1000), # reduced uncertainty yrange. (unit: (umol/m^2/s)^2)
        "posteriorMprior": (-50, 50), # difference between posterior and prior flux. (unit: umol/m^2/s)
    },

    "daily":
    {
        #-- flux unit: KtCO2/km^2/day --#
        #-- reduced uncertainty unit: (KtCO2/km^2/day)^2 --#
        "priorEmiss": (-5, 25),
        "posteriorEmiss": (-5, 25),
        "reducedUncertainty": (0, 0.01),
        "posteriorMprior": (-20, 20),
    },
    "weekly":
    {
        #-- flux unit: KtCO2/km^2/week --#
        #-- reduced uncertainty unit: (KtCO2/km^2/week)^2 --#
        "priorEmiss": (-5 * 7, 25 * 7),
        "posteriorEmiss": (-5 * 7, 25 * 7),
        "reducedUncertainty": (0 * 7, 0.01 * 7),
        "posteriorMprior": (-20 * 7, 20 * 7),
    },
    "monthly":
    {
        #-- flux unit: KtCO2/km^2/month --#
        #-- reduced uncertainty unit: (KtCO2/km^2/month)^2 --#
        "priorEmiss": (None, None),
        "posteriorEmiss": (None, None),
        "reducedUncertainty": (None, None),
        "posteriorMprior": (None, None),
    },
    "all":
    {
        "priorEmiss": (None, None),
        "posteriorEmiss": (None, None),
        "reducedUncertainty": (None, None),
        "posteriorMprior": (None, None),
    },
}

#--- 2D field plots keywords ---#
fieldPlots_kwargs = {

    "figsize": (10, 8), # figure size
    #"PlotsDir": benchDir,
    "PlotsName_Prefix": "spatial_temp", # pdf file name prefix ( pdf file name : <prefix>_<timeScale>.pdf )
    #"dpi": 200,
    "pcolormesh_kwargs": {},  # pyplot.pcolormesh keywords.
    "cbar_kwargs": {}, # pyplot.colorbar keywords.
    "spatial_kwargs": # BIS_BenchPlots.bis_spatial_plot keywords.
    {
        "cmap": "coolwarm", # colormap
        "cmap_va": False, 
        "draw_gridline": True, # draw grid lines.
        "draw_lonlat_label": True,# draw longitude & latitude label.
    },
    "priorEmiss_kwargs": # other keywords for prior flux.
    {
        #"vminmax": (None, None),
        "title": "prior", # title
    },
    "posteriorEmiss_kwargs": # other keywords for posterior flux.
    {
        #"vminmax": (None, None),
        "title": "posterior", # title
    },
    "reducedUncertainty_kwargs": # other keywords for reduced uncertainty.
    {
        #"vminmax": (None, None), 
        "title": "reduced uncertainty", # title
    },
    "posteriorMprior_kwargs":  # other keywords for posterior - prior.
    {
        #"vminmax": (None, None), 
        "title": "posterior - prior", # title
    }, 
}

#--- common plots setting for prior & poterior flux time series. ---#
emissPlots_kwargs = {
    "priorEmiss_kwargs": # prior flux time series plots keywords.
    {
        "plot_kwargs": {"color": "black"}, # pyplot.plot keywords.
        "spread_kwargs": {"facecolor": "grey", "alpha": 0.3}, # pyplot.fill_between keywords.
    },
    "posteriorEmiss_kwargs": # posterior flux time series plots keywords.
    {
        "plot_kwargs": {"color": "red"}, # pyplot.plot keywords.
        "spread_kwargs": {"facecolor": "red", "alpha": 0.3}, # pyplot.fill_between keywords.

    },
    "truthEmiss_kwargs": # truth flux time series plots keywords. ( for OSSE experiments ) 
    {
        "plot_kwargs": {"color": "blue", "linestyle": "--"}, # pyplot.plot keywords.
    },
}

#--- site time series plots keywords. ---#
site_timeSeries_kwargs = {

    "BenchName": "site_timeSeries_temp", # pdf file name.
    "figsize": (8, 11), # figure size.

    "station_local_kwargs":
    {
        "stationsPlot_kwargs":
        {
            "draw_gridline": True, # draw grid lines.
            "draw_lonlat_label": True,# draw longitude & latitude label.
        },
        "scatter_kwargs": 
        {
            "s": 1,
            "c": "red",  
        },
        "title": "",
    },

    "obs_timeSeries_kwargs": # obervations time series plots keywords.
    {
        "vrange": (350, 600), # yrange. (unit: ppm)
        "unit": "ppm", # unit label.
    },

    "plot_concentration_kwargs":  # site CO2 concentrations plots keywords.
    {
        #-- site observation data --#
        "obsConc_kwargs": 
        {
            "plot_kwargs": {"color": "blue", "marker": "o"} #pyplot.plot keywords.
        },
        #-- site inital concentration --#
        "initConc_kwargs": 
        {
            "plot_kwargs": {"color": "black"} #pyplot.plot keywords.
        },
        #-- site prior concentration --#
        "priorConc_kwargs": 
        {
            "plot_kwargs": {} #pyplot.plot keywords.
        },
        #-- site posterior (primary) concentration --#
        "postConc_kwargs": 
        {
            "plot_kwargs": {"color": "green"} #pyplot.plot keywords.
        },
        #-- site final concentration --#
        "finalConc_kwargs": 
        {
            "plot_kwargs": {"color": "red"} #pyplot.plot keywords.
        },
        "vrange": (350, 600), # vrange.
        "unit": "ppm", # unit label.
    },
    #-- Hx value plots --#
    "plot_Hx_kwargs": {
        #-- Hx prior keywords --#
        "priorHx_kwargs": 
        {
            "plot_kwargs": {"color": "black"} #pyplot.plot keywords.
        },
        #-- Hx posterior keywords --#
        "postHx_kwargs": 
        {
            "plot_kwargs": {"color": "green"} #pyplot.plot keywords.
        },
        #-- Hx final keywords --#
        "finalHx_kwargs": 
        {
            "plot_kwargs": {"color": "red"} #pyplot.plot keywords.
        },
        #-- horizontal line --#
        "yLine": 
        [
            (0, {"color": "grey", "linestyle": "--"}), # (y_value, pyplot.plot keywords)
        ],
        "vrange": (-10, 50), # vrange.
        "unit": "ppm", # unit label.
    },
    #-- background plots keywords --#
    "plot_bck_kwargs": {
        #-- prior background plots keywords --#
        "bckPrior_kwargs": 
        {
            "plot_kwargs": {"color": "black"}, # pyplot.plot keywords.
            "spread_kwargs": {"facecolor": "grey", "alpha": 0.3}, # pyplot.fill_between keywords.
        },
        #"sbckPrior_kwargs": {"facecolor": "grey", "alpha": 0.3},

        #-- posterior background plots keywords --#
        "bckProc_kwargs": 
        {
            "plot_kwargs": {"color": "red"},
            "spread_kwargs": {"facecolor": "red", "alpha": 0.3},
        },
        #"sbckProc_kwargs": {"facecolor": "red", "alpha": 0.3},
        "vrange": (350, 550),
        "unit": "ppm",
    },
    "plot_point_emiss_kwargs": {

        **emissPlots_kwargs,
        "vrange": (None, None),
        "unit": "umol/m^2/s",

    },
}

#--- define points/regions list. ---#
pointRegionList = ["IAPtower", "XiangHe", "Beijing", "Tianjin", "TWH-home"]
pointRegionList += ["pi424", "pi405", "pi431", "pi440", "pi445", "pi449", "pi471", "pi629", "pi637", "pi672", "pi679"]

#--- vrange of regional CO2 total flux. (unit: tCO2/timeUnit) ---#
region_vrange = {
    "mean": (-1.0, 2.0),
    "hourly": (0, 30000),
    #"hourly": (2500, 4000),
    "daily": (0 * 24, 30000 * 24),
    #"daily": (None, None),
    "weekly": (0 * 168, 30000 * 168),
}

#--- vrange of CO2 total flux in a point. (unit: tCO2/km^2/timeUnit) ---#
point_vrange = {
    #"mean": (-10, 100),
    "mean": (0, 30),
    "hourly": (-10, 20),
    "daily": (-10 * 24,  20 * 24),
    "weekly": (-10 * 168, 20 * 168),
}

#--- keywords of point_region_timeSeries plots function. ---#
point_region_timeSeries_kwargs = {

    "BenchName": "point_region_timeSeries",
    "figsize": (8, 11),
}

for prName in pointRegionList:
    point_region_timeSeries_kwargs[prName + "_kwargs"] = {
        **emissPlots_kwargs,
        "vrange": (None, None),
    }

#--- define the point / region location configure. ---#
point_region_kwargs = {

    "IAPtower": 
    {
        "type": "point", # one of "point" & "region".
        "location": (116.3667, 39.9667), # (longitude, latitude)
    },
    "XiangHe": 
    {
        "type": "point",
        "location": (116.9578, 39.7833),
    },
    "TWH-home":
    {
        "type": "point",
        "location": (116.33, 39.73),
    },
    "pi424":
    {
        "type": "point",
        "location": (117.1114, 40.1425),
    },
    "pi405":
    {
        "type": "point",
        "location": (116.2501, 39.7241),
    },
    "pi431":
    {
        "type": "point",
        "location": (116.1880, 39.7322),
    },
    "pi440":
    {
        "type": "point",
        "location": (116.6882, 39.8787),
    },
    "pi445":
    {
        "type": "point",
        "location": (116.3377, 39.9630),
    },
    "pi449":
    {
        "type": "point",
        "location": (116.4576, 40.2344),
    },
    "pi471":
    {
        "type": "point",
        "location": (115.7811, 39.5383),
    },
    "pi629":
    {
        "type": "point",
        "location": (116.6329, 40.3187),
    },
    "pi637":
    {
        "type": "point",
        "location": (116.5078, 39.7966),
    },
    "pi672":
    {
        "type": "point", 
        "location": (115.9797, 39.6419),
    },
    "pi679":
    {
        "type": "point",
        "location": (116.1343, 39.9042),
    },

    "Beijing":
    {
        "type": "region",
        "shpName": ["Beijing"], # region name in shape file.  (the name and directory of shape file is defined in variable "shpDir"
    },
    "Tianjin":
    {
        "type": "region",
        "shpName": ["Tianjin"],
    },
    "Heibei":
    {
        "type": "region",
        "shpName": ["Heibei"],
    },
    "BTH":
    {
        "type": "region",
        "shpName": ["Beijing", "Tianjin", "Heibei"],
    },

}

###############################################
#      special setting of each experiment     #
###############################################

#-- REAL experiment base --#
#-- tree (child --> parent) --# 
# BenchPlots_REAL -->  BIS_BenchPlots
expREAL_kwargs = {

}

expOSSE_kwargs = {
    "fieldPlots_vminmax": 
    {
        "mean":
        {
            "truthEmiss": (0, 80),
            "posteriorMtruth": (-50, 50)
        }, 
        "daily":
        {
            "truthEmiss": (0, 25),
            "posteriorMtruth": (-1, 1),
        },
        "weekly":
        {
            "truthEmiss": (0 * 7, 25 * 7),
            "posteriorMtruth": (-1 * 7, 1 * 7),
        },
        "monthly":
        {
            "truthEmiss": (None, None),
            "posteriorMtruth": (None, None),
        },
        "all":
        {
            "truthEmiss": (None, None),
            "posteriorMtruth": (None, None),
        },
    },

    "fieldPlots_kwargs":
    {
        "truthEmiss_kwargs":
        {
            "title": "Truth emiss",
        },
        "posteriorMtruth_kwargs":
        {
            "title": "Posterior - Truth",
        },
    },
    
    "movie_vminmax":
    {
        "truthEmiss": (0, 80),
        "posteriorMtruth": (-50, 50),
    },

    "movie_kwargs":
    {
        "truthEmiss_kwargs":
        {
            "unit": "umol/m^2/s",
            "title": "Truth",
        },

        "posteriorMtruth_kwargs":
        {
            "unit": "umol/m^2/s",
            "title": "Posterior - Truth",
        },
    },

    "site_plots_vrange":
    {
        "plot_obs": (-10, 150),
        "plot_concentration": (-10, 400),
        "plot_Hx": (-10, 150),
        "plot_bck": (-5, 5),
    },
    #"ffe_cfta_Plots_kwargs":
    #{
    #    "ffeEmiss_kwargs": 
    #    {
    #        "plot_kwargs": {"color": "brown"},
    #        "spread_kwargs": {"facecolor": "brown", "alpha": 0.3},
    #    },
    #    "cftaEmiss_kwargs": 
    #    {
    #        "plot_kwargs": {"color": "green"},
    #        "spread_kwargs": {"facecolor": "green", "alpha": 0.3},
    #    },
    #},

}


#-- OSSE experiment: const_emiss --#
#-- tree (child --> parent) --#
# BenchPlots_const_emiss --> BenchPlots_OSSE --> BIS_BenchPlots
expOSSE_const_emiss_kwargs = {

}


#-- REAL experiment: independent --#
#-- tree (child --> parent) --# 
# BenchPlots_independent --> BenchPlots_REAL --> BIS_BenchPlots
expREAL_independent_kwargs = {

}

#-- REAL experiment: independent_priorFnet --#
#-- tree (child --> parent) --# 
# BenchPlots_independent_priorFnet --> BenchPlots_independent --> BenchPlots_REAL --> BIS_BenchPlots
expREAL_independent_priorFnet_kwargs = {
    "fieldPlots_vminmax": 
    {
        "mean":
        {
            "ffeEmiss": (0, 80),
            "cftaEmiss": (-3, 3)
        }, 
        "daily":
        {
            "ffeEmiss": (0, 25),
            "cftaEmiss": (-1, 1),
        },
        "weekly":
        {
            "ffeEmiss": (0 * 7, 25 * 7),
            "cftaEmiss": (-1 * 7, 1 * 7),
        },
        "monthly":
        {
            "ffeEmiss": (None, None),
            "cftaEmiss": (None, None),
        },
        "all":
        {
            "ffeEmiss": (None, None),
            "cftaEmiss": (None, None),
        },
    },

    "fieldPlots_kwargs":
    {
        "ffeEmiss_kwargs":
        {
            "title": "Fossil fuel",
        },
        "cftaEmiss_kwargs":
        {
            "title": "Carbon flux to atmosphere",
        },
    },
    
    "movie_vminmax":
    {
        "ffeEmiss": (0, 80),
        "cftaEmiss": (-3, 3),
    },

    "movie_kwargs":
    {
        "ffeEmiss_kwargs":
        {
            "unit": "umol/m^2/s",
            "title": "FFE",
        },

        "cftaEmiss_kwargs":
        {
            "unit": "umol/m^2/s",
            "title": "Fta",
        },
    },
    "ffe_cfta_Plots_kwargs":
    {
        "ffeEmiss_kwargs": 
        {
            "plot_kwargs": {"color": "brown"},
            "spread_kwargs": {"facecolor": "brown", "alpha": 0.3},
        },
        "cftaEmiss_kwargs": 
        {
            "plot_kwargs": {"color": "green"},
            "spread_kwargs": {"facecolor": "green", "alpha": 0.3},
        },

    }

}
