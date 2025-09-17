#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 10:30:04 2025

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

from scipy.io import loadmat
import matplotlib as mpl


#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

#%% Indicate paths

figpath         = "/home/niu4/gliu8/figures/bydate/2025-09-16/"
proc.makedir(figpath)

datpath         = "/home/niu4/gliu8/projects/scrap/TP_crop/"

expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090"]

"""
Some issues to fix with vname

tx_surf_1m  (file name) corresponds to tx_sur
D20         (file name) corresponds to nz1
Dmax        (file name) corresponds to nz1

in one file
SSR (variable name) --> ssr, and the coordinate names are somewhat messed up

Also, it seems that I spelled separately wrong in the filename, RIP

"""

ninopath        = "/home/niu4/gliu8/projects/scrap/nino34/"
nexps           = len(expnames)

#%% Indicate the variable of choice and location of the map

vname       = "Dmaxgrad"#"ssr"
regpath     = "/home/niu4/gliu8/projects/scrap/TP_crop/regression/"


"""
TCo319_ssp585_Dmaxgrad_ENSO_regression_allmonths.nc
TCo319_ssp585_Dmaxgrad_ENSO_regression_seperatemonths.nc
"""

#%% Declare some functions for visualization

proj   = ccrs.PlateCarree(central_longitude=180)
projd  = ccrs.PlateCarree()
mons3  = proc.get_monstr()
mpl.rcParams['font.family'] = 'Montserrat'

def init_tp_map():
    bbplot = [120, 290, -20, 20]
    fix_lon = np.hstack([np.arange(120,190,10),np.arange(-180,-60,10)])
    proj   = ccrs.PlateCarree(central_longitude=180)
    projd  = ccrs.PlateCarree()
    fig,ax = plt.subplots(1,1,figsize=(12.5,4.5),subplot_kw={'projection':proj})
    ax.set_extent(bbplot)
    ax     = viz.add_coast_grid(ax,bbox=bbplot,fill_color='k',
                                proj=ccrs.PlateCarree(),fix_lon=fix_lon,ignore_error=True)
    return fig,ax


def get_timestr(ds):
    tstart =  str(np.array((ds.time.data[0])))[:4] #ds.time.isel(time=0)
    tend    = str(np.array((ds.time.data[-1])))[:4] #ds.time.isel(time=0)
    return "%s to %s" % (tstart,tend)



#%% Load the regression patterns

allmon = []
sepmon = []

for ex in tqdm.tqdm(range(nexps)):
    
    # All Months Regression
    ncname = "%s%s_%s_ENSO_regression_allmonths.nc" % (regpath,expnames[ex],vname,)
    ds     = xr.open_dataset(ncname).load()
    allmon.append(ds)
    
    # Separate Month Regression
    ncname = "%s%s_%s_ENSO_regression_seperatemonths.nc" % (regpath,expnames[ex],vname,)
    ds     = xr.open_dataset(ncname).load()
    sepmon.append(ds)

#%% (1) Examine the All-Month Pattern by Simulation



if vname == "sst":
    cints   = np.arange(-2,2.1,.1)
    vunit   = "$\degree C$"
    conversion = None
elif vname == "ssr":
    vunit       = "$W m^{-2}$"
    cints       = np.arange(-50,55,5)
    conversion  = 1/(3*3600) # 3 h accumulation time...? #1/(24*30*3600)
    # https://forum.ecmwf.int/t/surface-radiation-parameters-joule-m-2-to-watt-m-2/1588
elif vname == "str":
    vunit       = "$W m^{-2}$"
    cints       = np.arange(-20,22,2)
    conversion  = 1/(3*3600) # 3 h accumulation time...? #1/(24*30*3600)
elif vname == "tx_sur":
    vunit       = "$m s^{-2}$"
    cints       = np.arange(-0.02,0.021,0.002)#None#np.arange(-1,1,2)
    conversion  = None
elif vname == "D20" or vname == "Dmaxgrad":
    vunit      = "m"
    cints     = np.arange(-45,50,5)#np.arange(-24,26,2)
    conversion = None
    
    

for ex in range(nexps):
    dsplot  = allmon[ex]
    plotvar = dsplot.regression_pattern
    sigmask = dsplot.sigmask
    timestr = get_timestr(dsplot)
    
    if conversion is not None:
        plotvar = plotvar * conversion
    
    # Initialize Plot
    
    fig,ax  = init_tp_map()
    #ax.coastlines()
    
    
    
    # Plot Regression Map
    #pcm = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,cmap='cmo.balance',vmin=-2,vmax=2,transform=projd)
    if cints is None:
        pcm = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,cmap='cmo.balance',transform=projd)
    else:
        #pcm     = ax.contourf(plotvar.lon,plotvar.lat,plotvar,cmap='cmo.balance',levels=cints,transform=projd)
        #clbl    = ax.clabel(pcm,)
        
        pcm = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,cmap='cmo.balance',
                            vmin=cints[0],vmax=cints[-1],transform=projd)
        cl = ax.contour(plotvar.lon,plotvar.lat,plotvar,colors="k",levels=cints,transform=projd,linewidths=0.75)
        clbl    = ax.clabel(cl,)
        
        viz.add_fontborder(clbl,w=2.5)
    cb      = fig.colorbar(pcm,ax=ax,fraction=0.01,pad=0.01)
    cb.set_label(r"%s [%s per 1$\sigma$ $Ni\~no3.4$]" % (vname,vunit))
    
    # Add Significance Mark
    if ex <2: # Low resolution, ok to plot dots
        viz.plot_mask(sigmask.lon,sigmask.lat,sigmask.T,color='gray',markersize=.5,
                      ax=ax,geoaxes=True,proj=projd,reverse=False)
    else: # Plot every few dots
        dotint = 20
        viz.plot_mask(sigmask.lon.data[::dotint],sigmask.lat.data[::dotint],
                      sigmask.data[::dotint,::dotint].T,
                      color='gray',markersize=.5,
                      ax=ax,geoaxes=True,proj=projd,reverse=False)
        
    
    # Add title
    title = "%s regression (%s, %s), All Months" % (vname,expnames_long[ex],timestr)
    ax.set_title(title)
    
    
    #plt.show()
    
    figname = "%s%s_%s_ENSO_Regression_Pattern_AllMonths.png" % (figpath,expnames[ex],vname,)
    plt.savefig(figname,dpi=150,bbox_inches='tight')
    
    plt.close()
    #plt.show()


#%% (2) Perform separately for each month



for ex in range(nexps):
    
    for im in range(12):
        
        dsplot  = sepmon[ex]
        plotvar = dsplot.regression_pattern.isel(mon=im)
        sigmask = dsplot.sigmask.isel(mon=im)
        timestr = get_timestr(dsplot)
        
        # Initialize Plot
        
        if conversion is not None:
            plotvar = plotvar * conversion
            
            
        fig,ax  = init_tp_map()
        ax.coastlines()
        
        
        
        # Plot Regression Map
        #pcm = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,cmap='cmo.balance',vmin=-2,vmax=2,transform=projd)
        #pcm     = ax.contourf(plotvar.lon,plotvar.lat,plotvar,cmap='cmo.balance',levels=cints,transform=projd)
        #clbl    = ax.clabel(pcm,)
        
        pcm = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,cmap='cmo.balance',
                            vmin=cints[0],vmax=cints[-1],transform=projd)
        cl = ax.contour(plotvar.lon,plotvar.lat,plotvar,colors="k",levels=cints,transform=projd,linewidths=0.75)
        clbl    = ax.clabel(cl,)
        
        
        viz.add_fontborder(clbl,w=2.5)
        cb      = fig.colorbar(pcm,ax=ax,fraction=0.01,pad=0.01)
        cb.set_label("%s [%s per 1$\sigma$ $Ni\~no3.4$]" % (vname,vunit))
        
        # Add Significance Mark
        if ex <2: # Low resolution, ok to plot dots
            viz.plot_mask(sigmask.lon,sigmask.lat,sigmask.T,color='gray',markersize=.5,
                          ax=ax,geoaxes=True,proj=projd,reverse=False)
        else: # Plot every few dots
            dotint = 20
            viz.plot_mask(sigmask.lon.data[::dotint],sigmask.lat.data[::dotint],
                          sigmask.data[::dotint,::dotint].T,
                          color='gray',markersize=.5,
                          ax=ax,geoaxes=True,proj=projd,reverse=False)
            
        
        # Add title
        title = "%s %s regression (%s, %s)" % (mons3[im],vname,expnames_long[ex],timestr)
        ax.set_title(title)
        
        figname = "%s%s_%s_ENSO_Regression_Pattern_Mon%02i.png" % (figpath,expnames[ex],vname,im+1)
        plt.savefig(figname,dpi=150,bbox_inches='tight')

#for ex in range(nexps):
    


    
    





