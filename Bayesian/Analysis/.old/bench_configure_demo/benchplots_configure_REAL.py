#!/usr/bin/env python

shpDir = "/home/tangwh/datasets/china_shp/cnhimap.shp"


movie_kwargs = {

    "figsize": (10, 8),
    "snapDir": "/home/tangwh/modeling/BIS_v6.0/Bayesian/Analysis/test",
    "dpi": 200,
    "pcolormesh_kwargs": {}, 
    "cbar_kwargs": {},
    "spatial_kwargs":
    {
        "cmap": "coolwarm",
        "cmap_va": False, 
        "draw_gridline": True, 
        "draw_lonlat_label": True,
    }
    "priorEmiss_kwargs": 
    {
        "vminmax": (0, 60),
        "unit": "umol/m^2/s",
        "title": "prior",
    },
    "posteriorEmiss_kwargs": 
    {
        "vminmax": (0, 60),
        "unit": "umol/m^2/s",
        "title": "posterior",
    },
    "reducedUncertainty_kwargs": 
    {
        "vminmax": (0, 100), 
        "unit": "(umol/m^2/s)^2",
        "title": "reduced uncertainty",
    },
    "posteriorMprior_kwargs": 
    {
        "vminmax": (-30, 30), 
        "unit": "umol/m^2/s"
        "title": "posterior - prior"
    }, 
}

point_timeSeries_kwargs = {

    "BenchName": "point_timeSeries",
    "figsize": (8, 11),

    "plot_concentration_kwargs": {
        "obsConc_kwargs": {"color": "blue", "marker": "o"},
        "initConc_kwargs": {"color": "black"},
        "priorConc_kwargs": {},
        "postConc_kwargs": {"color": "green"},
        "finalConc_kwargs": {"color": "red"},
        "vrange": (350, 600),
    },
    "plot_Hx_kwargs": {
        "priorHx_kwargs": {"color": "black"},
        "postHx_kwargs": {"color": "green"},
        "finalHx_kwargs": {"color": "red"},
        "vrange": (-10, 50),
    },
    "plot_bck_kwargs": {
        "bckPrior_kwargs": {"color": "black"},
        "sbckPrior_kwargs": {"facecolor": "grey", "alpha": 0.3},
        "bckProc_kwargs": {"color": "red"},
        "sbckProc_kwargs": {"facecolor": "red", "alpha": 0.3},
        "vrange": (350, 550),
    },
    "plot_point_emiss_kwargs": {
        "priorEmiss_kwargs": {"color": "black"},
        "priorSigma_kwargs": {"facecolor": "grey", "alpha": 0.3},
        "posteriorEmiss_kwargs": {"color": "red"},
        "posteriorSigma_kwargs": {"facecolor": "red", "alpha": 0.3},
        "truthEmiss_kwargs": {"color": "blue", "linestyle": "--"},
        "vrange": (None, None),
    },
}