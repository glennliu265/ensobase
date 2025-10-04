#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualize ENSO Composites

Created on Wed Oct  1 10:22:30 2025

Written to run on Astraeus

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


#%%
amvpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/" # amv module
scmpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/"

sys.path.append(amvpath)
sys.path.append(scmpath)

from amv import proc,viz
import scm
import amv.loaders as dl
import cvd_utils as cvd

#%% Helper Functions

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
        for ax in axs:
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


#%% User Edits


datpath         = "/Users/gliu/Downloads/02_Research/01_Projects/07_ENSO/01_Data/TP_Crop/composites/"
expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090"]
vnames          = ["sst","ssr","str","tx_sur","D20","Dmaxgrad"]

vunits          = [r"$\degree C$",r"$\frac{W}{m^2}$",r"$\frac{W}{m^2}$",r"$\frac{m}{s^2}$","m","m"]
vnames_long     = ["SST","Surface Shortwave","Surface Longwave","Zonal Wind Stress","Thermocline (20$\degree$ Isotherm)","Thermocline (Max Vertical Gradient)"]


vmaxes          = [2,40,20,0.02,20,20]



ninoname        = [r"$El$ $Ni\tilde{n}o$",r"$La$ $Ni\tilde{n}a$"]
ninocol         = ["cornflowerblue","firebrick"]
ninoshort       = ['nino','nina']

# Conversion for STR and SSR considering 3h Accumulation
conversion  = 1/(3*3600) # 3 h accumulation time...? #1/(24*30*3600)
# https://forum.ecmwf.int/t/surface-radiation-parameters-joule-m-2-to-watt-m-2/1588

figpath         = "/Users/gliu/Downloads/02_Research/01_Projects/07_ENSO/02_Figures/20251003/"
proc.makedir(figpath)


# BBoxes from Jin et al. 2020 Eqn. 6.6 
bbox_cep        = [150,-130+360,-5,5] # central equatorial pacific, for [tau_x], 
bbox_nino3      = [150,-90+350,-5,5]  # Nino 3 Box: For SST, <tau_x>
bbox_epac       = [155,-80+360,-5,5]  # Eastern Pacific (for h_e calculation)
bbox_wpac       = [120,155,-5,5]      # Western Pacific (for h_w calculation)
 


#%% Plotting Things

proj  = ccrs.PlateCarree()
mons3 = proc.get_monstr()



#%% Load the composites computed in identify_enso_events

nexps = len(expnames)
nvars = len(vnames)

composites_byexp = []
for ex in range(nexps):
    
    composites_byvar = []
    for v in tqdm.tqdm(range(nvars)):
        ncname = "%sENSO_Composites_%s_%s.nc" % (datpath,vnames[v],expnames[ex])
        ds = xr.open_dataset(ncname)[vnames[v]].load()
        
        # For SSR and STR, do a conversion
        if vnames[v] in ['str','ssr']:
            
            conversion  = 1/(3*3600) # 3 h accumulation time...? #1/(24*30*3600)
            # https://forum.ecmwf.int/t/surface-radiation-parameters-joule-m-2-to-watt-m-2/1588
            ds = ds * conversion
        
        composites_byvar.append(ds)
    
    #ds_exp = xr.merge(composites_byvar,join='override',compat='override')
    #composites_byexp.append(ds_exp)
    composites_byexp.append(composites_byvar)

#%%


lags =  composites_byexp[ex][v].lags

fsz_title = 18
fsz_axis  = 14
fsz_tick  = 12

for ex in range(nexps): # Loop by Experiment and Event Type
    for ninotype in range(2):
        expdir  = "%s/%s_%s/" % (figpath,expnames[ex],ninoshort[ninotype])
        proc.makedir(expdir)
        
        itt = 0
        for lag in tqdm.tqdm(lags): # Generate Plot for each Lag
            fig,axs = init_tp_map(nrow=5,ncol=1,figsize=(12.5,16))
            
            for v in range(nvars-1): # Subplots for each variable 
                ax      = axs[v]
                ds      = composites_byexp[ex][v].isel(event=ninotype).sel(lags=lag)
                plotvar = ds
                #plotvar = ds
                
                pcm = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,cmap='cmo.balance',
                                    vmin=-vmaxes[v],vmax=vmaxes[v],transform=proj,
                                    )
                
                cb = fig.colorbar(pcm,ax=ax,fraction=0.01,pad=0.01)
                cb.ax.tick_params(labelsize=fsz_tick)
                cb.set_label(r"%s [%s]"  % (vnames[v],vunits[v]),fontsize=fsz_axis)
            
            plt.suptitle(r"AWI-CM3 %s %s Composite, Lag %i" % (expnames_long[ex],ninoname[ninotype],lag),y=0.90,fontsize=fsz_title)
            
            
            figname = "%s%s_%s_Composite_iter%0i_lag%02i" % (expdir,expnames[ex],ninoshort[ninotype],itt,lag)
            plt.savefig(figname,dpi=150,bbox_inches='tight')
            
            itt += 1
            

#%% DEBUG: Chec kwhat is going on with DART 2090 Nino


ex       = 3
ninotype = 1

testvars = []
for v in range(nvars):
    ds       = composites_byexp[ex][v].isel(event=ninotype)
    testvars.append(ds)


#%% Check Asymmetry by adding El Nino and La Nina Composites


composite_diff = []


for ex in range(nexps):
    diffs = []
    for v in range(nvars):
        
        ninocomp = composites_byexp[ex][v].sel(event='nino')
        ninacomp = composites_byexp[ex][v].sel(event='nina')
        diffs.append(ninocomp + ninacomp)
    composite_diff.append(diffs)
        
#%% Plot Composite Differences for each variable (for all lags)
   
        
for ex in range(nexps): # Loop by Experiment and Event Type
    expdir = "%s/%s_%s/" % (figpath,expnames[ex],"composite_diff")

    proc.makedir(expdir)
    
    itt = 0
    for lag in tqdm.tqdm(lags): # Generate Plot for each Lag
        fig,axs = init_tp_map(nrow=5,ncol=1,figsize=(12.5,16))
        
        
        for v in range(nvars-1): # Subplots for each variable 
            ax      = axs[v]
            ds      = composite_diff[ex][v].sel(lags=lag)
            plotvar = ds
            #plotvar = ds
            
            pcm = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,cmap='cmo.balance',
                                vmin=-vmaxes[v]/2,vmax=vmaxes[v]/2,transform=proj,
                                )
            
            cb = fig.colorbar(pcm,ax=ax,fraction=0.01,pad=0.01)
            cb.ax.tick_params(labelsize=fsz_tick)
            cb.set_label(r"%s [%s]"  % (vnames[v],vunits[v]),fontsize=fsz_axis)
        
        plt.suptitle(r"AWI-CM3 %s Composite Sum (Nino + Nina), Lag %i" % (expnames_long[ex],lag),y=0.90,fontsize=fsz_title)
        
        
        figname = "%s%s_Composite_Diff_iter%0i_lag%02i" % (expdir,expnames[ex],itt,lag)
        plt.savefig(figname,dpi=150,bbox_inches='tight')
        
        itt += 1

#%% Also Make Hovmuller plots

ex     = 1
#v      = 3
ninoid = 0
latavg = 10

for ex in range(nexps):
    for ninoid in [0,1]:
        fig,axs   = plt.subplots(5,1,constrained_layout=True,figsize=(12.5,16))
        
        for v in range(nvars-1):
            ax      = axs[v]
        
            plotvar = composites_byexp[ex][v].isel(event=ninoid).sel(lat=slice(-latavg,latavg)).mean('lat')
            
            
            pcm     = ax.pcolormesh(plotvar.lon,plotvar.lags,plotvar,vmin=-vmaxes[v],vmax=vmaxes[v],cmap='cmo.balance')
            
            ax.set_ylabel("Lag (Months)",fontsize=fsz_axis)
            ax.set_xlabel("Longitude",fontsize=fsz_axis)
            ax.tick_params(labelsize=fsz_tick)
            
            cb = fig.colorbar(pcm,ax=ax,fraction=0.01,pad=0.01)
            cb.ax.tick_params(labelsize=fsz_tick)
            cb.set_label(r"%s [%s]"  % (vnames[v],vunits[v]),fontsize=fsz_axis)
            
            #cb.set_label(r"%s %s: %s [%s]"  % (expnames_long[ex],ninoname[ninoid],vnames[v],vunits[v]),fontsize=fsz_axis)
        plt.suptitle("AWI-CM3 %s $Hovm\ddot{o}ller$ (%s)\n -%s to %s Latitude Average" % (expnames_long[ex],ninoname[ninoid],latavg,latavg),
                     fontsize=fsz_title,y=1.05)
        
        figname = "%s%s_Hovmoller_%s.png" % (figpath,expnames[ex],ninoshort[ninoid])
        plt.savefig(figname,dpi=150,bbox_inches='tight')
        