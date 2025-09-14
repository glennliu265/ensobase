#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 10:31:12 2025

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


#%% Load Data and establish functions

expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090"]
vname           = "sst"
dsall           = [xr.open_dataset("%s%s_%s.nc" % (datpath,ex,vname)).load() for ex in expnames]

def rename_time_dart(ds):
    if 'time_counter' in list(ds.coords):
        print("Renaming [time_counter] to [time]")
        ds = ds.rename({'time_counter':'time'})
    return ds

def movmean(ds,win):
    return np.convolve(ds.data,np.ones(win)/win,mode='same')

def remake_da(ds,dsref,name='sst'):
    coords = {'time':dsref.time}
    ds_new = xr.DataArray(ds,coords=coords,dims=coords,name=name)
    return ds_new

#%%

dsall         = [rename_time_dart(ds) for ds in dsall]
sst_anom      = [proc.xrdeseason(ds) for ds in dsall]

#%% Calculate ENSO Index
apply_movmean = False

bbox_nino34   = [-170+360,-120+360,-5,5]
dsall_reg     = [proc.sel_region_xr(ds,bbox_nino34) for ds in dsall]
nino_sim      = [proc.area_avg_cosweight(ds.sst) for ds in dsall_reg]

nino_sim      = [rename_time_dart(ds) for ds in nino_sim]

nino_ds       = [proc.xrdeseason(ds) for ds in nino_sim] # Remove scycle (anomalize)
#nino_ds       = [proc.xrdetrend(ds) for ds in nino_ds] # REmove linear trend
nino_ds       = [proc.xrdetrend_1d(ds,2) for ds in nino_ds] # REmove linear trend

if apply_movmean:
    nino34sim = [movmean(ds,5) for ds in nino_ds] # 5-month running mean
else:
    nino34sim = nino_ds

nino34_norm = [ds/np.nanstd(ds) for ds in nino34sim]
nino34_ds   = [remake_da(nino34_norm[ii],sst_anom[ii]) for ii in range(len(expnames))]

#%% Look at composites
bbplot = proc.get_bbox(dsall[0])

def init_tp_map():
    bbplot = [120, 290, -20, 20]
    proj   = ccrs.PlateCarree(central_longitude=180)
    projd  = ccrs.PlateCarree()
    fig,ax = plt.subplots(1,1,figsize=(12.5,4.5),subplot_kw={'projection':proj})
    ax     = viz.add_coast_grid(ax,bbox=bbplot,fill_color='k',proj=ccrs.PlateCarree())
    return fig,ax

ex = 0
#monsel = [0,1,2]

for ex in range(4):
    
    idin    = nino34_norm[ex]
    projd   = ccrs.PlateCarree()
    
    # Nino Composite
    fig,ax  = init_tp_map()
    ninoid  = idin >= 1
    plotsst = sst_anom[ex].sst.data[ninoid,:,:].mean(0)
    lon     = sst_anom[ex].lon
    lat     = sst_anom[ex].lat
    pcm     = ax.pcolormesh(lon,lat,plotsst,vmin=-2,vmax=2,
                            cmap='cmo.balance',
                            transform=projd)
    ax.set_title("%s El Nino Composite (n=%i/%i)" % (expnames_long[ex],ninoid.sum(),len(ninoid)))
    
    fig.colorbar(pcm,ax=ax,fraction=0.01,pad=0.01)
    figname = "%s%s_Nino_Composite.png" % (figpath,expnames[ex])
    plt.savefig(figname,dpi=150,bbox_inches='tight')
    
    
    # Nina Composite
    fig,ax  = init_tp_map()
    ninoid  = idin <= -1
    plotsst = sst_anom[ex].sst.data[ninoid,:,:].mean(0)
    lon     = sst_anom[ex].lon
    lat     = sst_anom[ex].lat
    pcm     = ax.pcolormesh(lon,lat,plotsst,vmin=-2,vmax=2,
                            cmap='cmo.balance',
                            transform=projd)
    ax.set_title("%s La Nina Composite (n=%i/%i)" % (expnames_long[ex],ninoid.sum(),len(ninoid)))
    
    fig.colorbar(pcm,ax=ax,fraction=0.01,pad=0.01)
    figname = "%s%s_Nina_Composite.png" % (figpath,expnames[ex])
    plt.savefig(figname,dpi=150,bbox_inches='tight')


#%% Plot Nino Indices

for ex in range(4):
    
    fig,ax  = plt.subplots(1,1,figsize=(8,3),constrained_layout=True)
    plot_ts = nino34_norm[ex]
    plotx   = nino_sim[ex].time
    stdev   = np.nanstd(nino_sim[ex])
    ax.plot(plotx,plot_ts)
    ax.set_ylim([-3.5,3.5])
    ax.axhline([-1],lw=0.65,ls='dashed',color='blue')
    ax.axhline([1],lw=0.65,ls='dashed',color='r')
    ax.axhline([0],lw=0.15,ls='solid',color='k')
    ax.set_title("%s Nino3.4 Index ($\sigma=%.2f \degree C$)" % (expnames_long[ex],stdev))
    figname = "%s%s_Nino34_Index.png" % (figpath,expnames[ex])
    plt.savefig(figname,dpi=150,bbox_inches='tight')
    
    
#%% Plot Hovmullers by looking at the average between 5S and 5N

sst_lonavg = [ds.mean('lat') for ds in sst_anom]

for ex in range(4):
    
    fig,ax  = plt.subplots(1,1,figsize=(8,6),constrained_layout=True)
    plotvar = sst_lonavg[ex].sst
    pcm     = plt.pcolormesh(plotvar.lon,plotvar.time,plotvar,
                             vmin=-2,vmax=2,cmap='cmo.balance')
    
    ax.set_xlabel("Longitude ($\degree E)$")
    ax.set_ylabel("Year")
    
    cb = viz.hcbar(pcm,ax=ax,fraction=0.1)
    cb.set_label("%s Average Tropical Pacific Temperature (5S to 5N)" % expnames_long[ex])
    
    
#%% Check the monthly variance of the ENSO Index

monvar = [ds.groupby('time.month').var('time') for ds in nino34_ds]

mons3 = proc.get_monstr()
for ii in range(4):
    
    fig,ax = viz.init_monplot(1,1,)
    ax.bar(mons3,monvar[ii])
    
    ax.set_ylim([0,2])
    ax.set_title(expnames_long[ii])
    
    ax.set_xlim([-1,12])

#%%
    
#%%




    
    




