#!/usr/bin/env python

ShpDir = "/home/tangwh/datasets/china_shp/cnhimap.shp"

movieDir = "/home/tangwh/modeling/BIS_v6.0/test/bench_test/test_dir"

movie_kwargs = {
    "figsize": (8, 11),
	"movieDir": movieDir,
	"movieName": "spatial_movie",
    "snapDir": movieDir + "/.movie_snap",
    "dpi": 200,
    "vminmax": {
        "priorEmiss": (0, 10),
        "posteriorEmiss": (0, 10),
        "reducedUncertainty": (0, 10),
        "posteriorMprior": (-10, 10),
        "truthEmiss": (0, 10),
        "priorMtruth": (-10, 10),
    },
    "cmap": "coolwarm",
    "pcolormesh_kwargs": {}, 
    "cbar_kwargs": {},

}
