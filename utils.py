#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Shared functions used by scripts in ensobase repository

Currently runs on niu (need to load [amv] dependencies...)
add this line to the top of the script to import module:

----
ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut
----

Function Types (will reorganize later):
    (A)     : AWI-CM3 specific preprocessing
    (c)     : calculations/analysis
    (g)     : general/universal function, can move to proc
    (l)     : loader
    (v)     : visualization

Function                Description
--------                -----------
awi_mean_loader     : (l) load mean/monvar/scycle calculations from calc_mean_patterns_TP 
init_tp_map         : (v) initialize tropical Pacific plot 
swap_rename         : (g) check if variable exists and rename if so
standardize_names   : (A) uses swap_rename to replace variable and dimension names in AWI_CM3 output 
varcheck            : (A) checks and converts variables for AWI-CM3




Created on Wed Oct  8 15:26:57 2025

@author: gliu
"""

import numpy as np
import xarray as xr
import cartopy
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys


#%%

amvpath = "/home/niu4/gliu8/scripts/commons"
sys.path.append(amvpath)
from amv import proc,viz

#%%



def awi_mean_loader(expname,vname,calcname,outpath=None):
    if outpath is None:
        outpath = "/home/niu4/gliu8/projects/scrap/TP_crop/summary/"
    ncname = "%s%s_%s_%s.nc" % (outpath,expname,vname,calcname)
    ds = xr.open_dataset(ncname).load()
    return ds

def init_tp_map(nrow=1,ncol=1,figsize=(12.5,4.5),ax=None):
    bbplot = [120, 290, -20, 20]
    fix_lon = np.hstack([np.arange(120,190,10),np.arange(-180,-60,10)])
    proj   = ccrs.PlateCarree(central_longitude=180)
    projd  = ccrs.PlateCarree()
    
    
    if ax is None:
        fig,axs = plt.subplots(nrow,ncol,figsize=figsize,subplot_kw={'projection':proj})
        newfig = True
    else:
        newfig = False
    if nrow != 1 or ncol != 1:
        for ax in axs.flatten():
            ax.set_extent(bbplot)
            ax     = viz.add_coast_grid(ax,bbox=bbplot,fill_color='k',
                                        proj=ccrs.PlateCarree(),fix_lon=fix_lon,ignore_error=True)
        ax = axs
    else:
        ax = axs
        ax.set_extent(bbplot)
        ax     = viz.add_coast_grid(ax,bbox=bbplot,fill_color='k',
                                    proj=ccrs.PlateCarree(),fix_lon=fix_lon,ignore_error=True)

    if newfig:
        return fig,ax
    return ax

def standardize_names(ds):
    
    ds = swap_rename(ds,'time_counter','time')
    ds = swap_rename(ds,"TIME_COUNTER",'time')
    ds = swap_rename(ds,"LON","lon")
    ds = swap_rename(ds,"LAT","lat")
    ds = swap_rename (ds,"LAT232_409","lat")
    
    # Other preprocessing
    # drop LON_bnds, TIME_COUNTER_bnds
    dropvars = ["LON_bnds","TIME_COUNTER_bnds"]
    for dropvar in dropvars:
        if dropvar in ds:
            ds = ds.drop_vars(dropvar)
    return ds

def swap_rename(ds,chkvar,newvar):
    if chkvar in list(ds.coords):
        print("Renaming [%s] to [%s]" % (chkvar,newvar))
        ds = ds.rename({chkvar:newvar})
    return ds

def varcheck(ds,vname,expname):
    if np.any(ds > 273) and vname == "sst": # Convert to Celsius
        print("Converting from Kelvin to Celsius for %s" % expname)
        ds = ds - 273.15
        
    if vname in ['str','ssr','strc','ssrc','ttr','tsr','ttrc','tsrc','sshf','slhf']: # Accumulation over 3h
        # Conversion for STR and SSR considering 3h Accumulation
        if "TCo319" in expname:
            print("Correction for accumulation over 6 hours for %s" % expname)
            accumulation_hr = 6
        else:
            print("Correction for accumulation over 3 hours for %s" % expname )
            accumulation_hr = 3
        conversion  = 1/(3600 * accumulation_hr)  # 3 h accumulation time...? #1/(24*30*3600)
        # https://forum.ecmwf.int/t/surface-radiation-parameters-joule-m-2-to-watt-m-2/1588
        ds          = ds * conversion
        
    if vname in ["cp","lsp"]: # Convert from [meters/accumulation period] to [mm/day]
        if "TCo319" in expname:
            print("Correction for accumulation over 6 hours for %s" % expname)
            accumulation_hr = 6
        else:
            print("Correction for accumulation over 3 hours for %s" % expname )
            accumulation_hr = 3
        conversion = (24/accumulation_hr) * 1000
    return ds

