#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Check EIS and Cloud Cover

Created on Tue Nov  4 11:25:13 2025

@author: gliu

"""

import sys

import time
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import xarray as xr
import sys
import tqdm
import glob 
import scipy as sp
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
import cartopy

from scipy.io import loadmat
import matplotlib as mpl

import climlab

#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% User Edits


expnames            = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
expnames_long       = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"]

ensoid_name         = "nino34"
standardize         = False
regrid_1x1          = True

figpath             = "/home/niu4/gliu8/figures/bydate/2025-11-04/"
proc.makedir(figpath)


#%% Load regridded EIS

datpath     = "/home/niu4/gliu8/projects/scrap/regrid_1x1/"

ds9km       = xr.open_dataset(datpath + "TCo1279-DART-1950_eis_regrid1x1.nc").load()
ds5km       = xr.open_dataset(datpath + "TCo2559-DART-1950C_eis_regrid1x1.nc").load()

dsall       = [ds9km,ds5km]
expids      = [2,4]

#%% Other Functions

def init_globalmap(nrow=1,ncol=1,figsize=(12,8)):
    proj            = ccrs.Robinson(central_longitude=-180)
    bbox            = [-180,180,-90,90]
    fig,ax          = plt.subplots(nrow,ncol,subplot_kw={'projection':proj},figsize=figsize,constrained_layout=True)
    
    multiax = True
    if (type(ax) == mpl.axes._axes.Axes) or (type(ax) == cartopy.mpl.geoaxes.GeoAxes):
        ax = [ax,]
        multiax = False
    
    for a in ax:
        a.coastlines(zorder=10,lw=0.75,transform=proj)
        a.gridlines(ls ='dotted',draw_labels=True)
        
    if multiax is False:
        ax = ax[0]
    return fig,ax

#%% Calcualtions

# Compute Time Mean
dsmean = [ds.mean('time').eis for ds in dsall]

#%% Plot Time Mean



#%% Visualize the time-mean eis

proj    = ccrs.PlateCarree()
cints   = np.arange(-8,8.5,0.5)


fig,axs = init_globalmap(2,1)

for ii in range(2):
    ax      = axs[ii]
    plotvar = dsmean[ii]
    ex      = expids[ii]
    
    pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,vmin=cints[0],vmax=cints[-1],cmap='cmo.balance')
    cl      = ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=cints,
                          colors="k",linewidths=0.75)
    
    
    ax.set_title(expnames_long[ex])

# fig,ax  = ut.init_
# plotvar = eismean
# pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,vmin=cints[0],vmax=cints[-1],cmap='cmo.balance')
# cl      = ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=cints,
#                       colors="k",linewidths=0.75)

cb      = viz.hcbar(pcm,ax=ax,fraction=0.045,pad=0.1)
cb.set_label("Mean Estimated Inversion Strength")
ax.clabel(cl,levels=cints[::2],fontsize=8)

figname = "%sMeanEIS_Comparison_Global.png" % figpath
plt.savefig(figname,dpi=150,bbox_inches='tight')
plt.show()

#%% Look at Differences in mean state

diff = dsmean[1] - dsmean[0] # 5km - 9km


fig,ax  = init_globalmap(1,1)
plotvar = diff
pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,vmin=cints[0],vmax=cints[-1],cmap='cmo.balance')
cl      = ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=cints,
                      colors="k",linewidths=0.75)

cb      = viz.hcbar(pcm,ax=ax,fraction=0.045,pad=0.1)
cb.set_label("Difference (5km - 9km) (EIS)")
ax.clabel(cl,levels=cints[::2],fontsize=8)

figname = "%sMeanEIS_Comparison_9km_5km.png" % figpath
plt.savefig(figname,dpi=150,bbox_inches='tight')
plt.show()

#%% Load in ENSO Indices

bbox_tropics    = [0        , 360      , -30,30] # Tropics (from Ceppi and Fueglistaler 2021)
dstrop          = [proc.sel_region_xr(ds,bbox_tropics) for ds in dsall]
dstropmean      = [proc.area_avg_cosweight(ds) for ds in dstrop]
dsgmean         = [proc.area_avg_cosweight(ds) for ds in dsall]

#%% Plot the Tropical and Global Mean EIS

fig,axs = plt.subplots(2,1,constrained_layout=True,figsize=(12.5,6))

ax = axs[0]
ax.set_title("Tropical Mean EIS")
for ii in range(2):
    ex = expids[ii]
    plotvar = dstropmean[ii].eis
    ax.plot(plotvar.time,plotvar,label=expnames_long[ex])

ax = axs[1]
ax.set_title("Global Mean EIS")
for ii in range(2):
    ex = expids[ii]
    plotvar = dsgmean[ii].eis
    ax.plot(plotvar.time,plotvar,label=expnames_long[ex])
plt.show()

ax.legend()
    
figname = "%sTropGlobMeanEIS_Timeseries.png" % figpath
plt.savefig(figname,dpi=150,bbox_inches='tight')

#%% Also Plot the mean seasonal cycle

mons3   = proc.get_monstr()
fig,axs = viz.init_monplot(2,1,constrained_layout=True,figsize=(8,6))

ax = axs[0]
ax.set_title("Tropical Mean EIS")
for ii in range(2):
    ex = expids[ii]
    plotvar = dstropmean[ii].eis.groupby('time.month').mean('time')
    ax.plot(mons3,plotvar,label=expnames_long[ex])

ax = axs[1]
ax.set_title("Global Mean EIS")
for ii in range(2):
    ex = expids[ii]
    plotvar = dsgmean[ii].eis.groupby('time.month').mean('time')
    ax.plot(mons3,plotvar,label=expnames_long[ex])
plt.show()

ax.legend()
    
figname = "%sTropGlobMeanEIS_Scycle.png" % figpath
plt.savefig(figname,dpi=150,bbox_inches='tight')