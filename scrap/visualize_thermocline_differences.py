#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Using D20 and Dmax calculated from "check_thermocline",
visualize some of the thermocline differences




This script runs locally on Astraeus
Copied Upper section from explore_enso

Created on Fri Sep 12 13:47:42 2025

@author: gliu
"""


from tqdm import tqdm
from scipy import signal

import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr

import sys
import cmocean
import time
import glob

#%% Import Packages

amvpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/" # amv module
scmpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/"

sys.path.append(amvpath)
sys.path.append(scmpath)

from amv import proc,viz
import scm
import amv.loaders as dl
import cvd_utils as cvd

#%% Indicate Paths

figpath = "/Users/gliu/Downloads/02_Research/01_Projects/05_SMIO/02_Figures/20250909/"
datpath = "/Users/gliu/Downloads/02_Research/01_Projects/07_ENSO/01_Data/TP_Crop/"
proc.makedir(figpath)

#%% Make some necessary functions/plotting variables

proj   = ccrs.PlateCarree(central_longitude=180)
projd  = ccrs.PlateCarree()
mons3  = proc.get_monstr()

def init_tp_map():
    bbplot = [120, 290, -20, 20]
    proj   = ccrs.PlateCarree(central_longitude=180)
    projd  = ccrs.PlateCarree()
    fig,ax = plt.subplots(1,1,figsize=(12.5,4.5),subplot_kw={'projection':proj})
    #ax     = viz.add_coast_grid(ax,bbox=bbplot,fill_color='k',proj=ccrs.PlateCarree())
    return fig,ax

#%% Load Data and establish functions

expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090"]
vnames          = ["D20","Dmaxgrad"]

#%% Load the Data
nexps = 4
nvars = 2

dall  = []
for vv in range(nvars):
    
    byexp = []
    for ex in range(nexps):
        
        ncname  = "%s%s_%s.nc" % (datpath,expnames[ex],vnames[vv])
        ds      = xr.open_dataset(ncname).load()
        byexp.append(ds)
    
    dall.append(byexp)

d20s,dmaxes = dall
#%% 


d20scycle   = [ds.groupby('time.season').mean('time') for ds in d20s]
dmaxscycle  = [ds.groupby('time.season').mean('time') for ds in dmaxes]


#%% Plot Differences Across Each Season

ex  = -2
sid = 0

cints = np.arange(-50,55,5)

for ex in range(4):
    for sid in range(4):
        
        
        fig,ax  = init_tp_map()
        ax.coastlines()
        plotvar = d20scycle[ex].isel(season=sid).nz1 - dmaxscycle[ex].isel(season=sid).nz1
        pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,vmin=-50,vmax=50,
                                cmap='cmo.balance',transform=projd)
        
        # pcm    = ax.contourf(plotvar.lon,plotvar.lat,plotvar,levels=cints,
        #                         cmap='cmo.balance',transform=projd)
        cl = ax.contour(plotvar.lon,plotvar.lat,plotvar,levels=cints,
                                colors="k",linewidths=0.75,transform=projd)
        ax.clabel(cl)
        
        
        fig.colorbar(pcm,ax=ax,fraction=0.010)
        ax.set_title("D20 - Dmax (%s, %s) [meters]" % (plotvar.season.data,expnames_long[ex]))
        
        
        figname = "%sD20_v_DMax_%s_%s.png" % (figpath,expnames[ex],plotvar.season.data.item())
        plt.savefig(figname,dpi=150,bbox_inches='tight')
        print(figname)


#%% Look at the seasonal cycle of Thermocline at a location









