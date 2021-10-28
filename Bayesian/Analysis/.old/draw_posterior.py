#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import datetime as dtm
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from matplotlib.colors import ListedColormap
from pdb import set_trace


str_StartDT = "2021-01-03_00:00:00"
str_EndDT = "2021-01-06_00:00:00"
dtHrs = 1
ShpDir = "/home/tangwh/datasets/china_shp/cnhimap.shp"
#InDir = "/home/tangwh/modeling/BIS_cases/test_lt=6_ls=10/output/processing"
#InDir = "/home/tangwh/modeling/BIS_cases/test1/output/processing"
#InDir = "/home/tangwh/modeling/BIS_cases/test_lt=6_ls=30/output/processing"
#InDir = "/home/tangwh/modeling/BIS_cases/osse2_lt=6_ls=30/output/processing"
InDir = "/home/tangwh/modeling/BIS_cases/Exp_OSSE=zero_one_Lt=48_Ls=1000/output/processing"
OutDir = "./output"
InPrefix = "BISOUT"
OutPrefix = "POS"

Start = dtm.datetime.strptime(str_StartDT, "%Y-%m-%d_%H:%M:%S")
End = dtm.datetime.strptime(str_EndDT, "%Y-%m-%d_%H:%M:%S")
dt = dtm.timedelta(hours = dtHrs)

def alpha_vary_cmap(colormap):
    cmap = plt.get_cmap(colormap)
    my_cmap = cmap(np.arange(cmap.N))
    my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
    my_cmap = ListedColormap(my_cmap)
    return my_cmap


Current = Start
emiss_total = 0.0
Ntotal = 0

while( Current <= End ):

    print(Current.strftime("%Y-%m-%d_%H:%M:%S"))
    ds = xr.open_dataset(InDir + "/" + InPrefix + "_" + Current.strftime("%Y-%m-%d_%H:%M:%S") + ".nc")
    emiss = ds["emiss"].values
    emiss_total += emiss
    Ntotal += 1
    lon = ds["lon"].values
    lat = ds["lat"].values
    lon_s = lon[0]
    lon_e = lon[-1]
    lat_s = lat[0]
    lat_e = lat[-1]
    fig = plt.figure()
    proj = ccrs.PlateCarree()
    ax = fig.add_subplot(1, 1, 1, projection = proj)
    ax.add_geometries( Reader(ShpDir).geometries(), proj, facecolor = "none", edgecolor = "k", linewidth = 1)
    ax.set_extent([lon_s, lon_e, lat_s, lat_e], crs = proj)
    gl = ax.gridlines(crs = proj, linestyle = "--", alpha = 0.5, draw_labels = True) 
    gl.top_labels = False
    gl.right_labels = False
    #cmap = alpha_vary_cmap("jet")
    cmap = "coolwarm"
    cs = ax.pcolormesh(lon, lat, emiss - 1, cmap = cmap, vmin = -1, vmax = 1)
    cbar = fig.colorbar(cs, orientation = "horizontal")
    cbar.set_label("umol/m^2/s")
    ax.set_title("Posterior-Truth " + Current.strftime("%Y-%m-%d_%H:%M:%S"))
    #set_trace()
    plt.savefig(OutDir + "/" + OutPrefix + "_" + Current.strftime("%Y-%m-%d_%H:%M:%S") + ".png")
    #plt.close(fig)
    Current += dt

emiss_mean = emiss_total / Ntotal
fig = plt.figure()
proj = ccrs.PlateCarree()
ax = fig.add_subplot(1, 1, 1, projection = proj)
ax.add_geometries( Reader(ShpDir).geometries(), proj, facecolor = "none", edgecolor = "k", linewidth = 1)
ax.set_extent([lon_s, lon_e, lat_s, lat_e], crs = proj)
gl = ax.gridlines(crs = proj, linestyle = "--", alpha = 0.5, draw_labels = True)
gl.top_labels = False
gl.right_labels = False
#cmap = alpha_vary_cmap("jet")
cmap = "coolwarm"
cs = ax.pcolormesh(lon, lat, emiss_mean - 1, cmap = cmap, vmin = -1, vmax = 1)
cbar = fig.colorbar(cs, orientation = "horizontal")
cbar.set_label("umol/m^2/s")
ax.set_title("Posterior-Truth mean")
plt.savefig(OutDir + "/" + OutPrefix + "_mean.png")
