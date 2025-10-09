#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Visualize Mean Patterns (Tropical Pacific)

- Loads output from calc_mean_patterns_TP and visualizes them

Created on Tue Oct  7 18:31:22 2025

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

#%% 

# Simulation Names -----
expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"]

timecrops       = [[1950,2100],None,None,None,None]


figpath         = "/home/niu4/gliu8/figures/bydate/2025-10-14/"
proc.makedir(figpath)

#expnames_all    = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]


#%% Mean State Variable PLots

sstdict = dict(
    vname      = "sst",
    longname   = "sea surface temperature",
    units      = "$\degree C$",
    cints_mean = np.arange(18,35.5,.5),
    cmap_mean  = "cmo.balance",
    )


sshdict = dict(
    vname      = "ssh",
    longname   = "sea surface height",
    units      = "$meters$",
    cints_mean = np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.deep_r",
    )


strdict = dict(
    vname      = "str",
    longname   = "Surface Thermal Radiation",
    units      = r"$\frac{J}{m^2}$",
    cints_mean = np.arange(-100,-22.5,2.5),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "plasma_r",
    )

ssrdict = dict(
    vname      = "ssr",
    longname   = "Surface Solar Radiation",
    units      = r"$\frac{J}{m^2}$",
    cints_mean = np.arange(120,290,10),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.solar",
    )


txsurdict = dict(
    vname      = "tx_sur",
    longname   = "Zonal Wind Stress to Ocean",
    units      = r"$\frac{m}{s^2}$",
    cints_mean = np.arange(-0.14,0.15,0.01),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.curl",
    )


# indicate variable and corresponding dictionary 
vname = "lcc"
indict = txsurdict#None#sstdict


#%%
datpath     = "/home/niu4/gliu8/projects/scrap/TP_crop/"
outpath     = "/home/niu4/gliu8/projects/scrap/TP_crop/summary/"
calcnames    = ["mean","scycle","var","monvar"]

nexps       = len(expnames)
ncalc       = len(calcnames)

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


#%%
    



ds_var = []
for ex in tqdm.tqdm(range(nexps)):
    
    
    metrics = []
    for nn in range(ncalc):
        try:
            ds = awi_mean_loader(expnames[ex],vname,calcnames[nn],outpath=outpath)
        except:
            print("Could not find %s for %s" % (calcnames[nn],expnames_long[ex]))
            metrics.append(None)
            
        
        metrics.append(ds)
    ds_var.append(metrics)

#%% Plot Time Mean Pattern with changes across resolution

projd  = ccrs.PlateCarree()


def varcheck(ds,vname,expname):
    if np.any(ds > 273) and vname == "sst": # Convert to Celsius
        print("Converting from Kelvin to Celsius for %s" % expname)
        ds = ds - 273.15
        
    if vname == "str" or vname == "ssr": # Accumulation over 3h
        # Conversion for STR and SSR considering 3h Accumulation
        if "TCo319" in expname:
            print("Correction for accumulation over 6 days for %s" % expname)
            accumulation_days = 6
        else:
            print("Correction for accumulation over 3 days for %s" % expname )
            accumulation_days = 3
        conversion  = 1/(3600 * accumulation_days)  # 3 h accumulation time...? #1/(24*30*3600)
        # https://forum.ecmwf.int/t/surface-radiation-parameters-joule-m-2-to-watt-m-2/1588
        ds          = ds * conversion
        
    return ds

vunit = "$\degree C$"
#cints = np.arange(255,305)
for ex in range(nexps):
    
    fig,ax  = init_tp_map()
    try:
        plotvar = ds_var[ex][0][vname]
    except:
        print("Skipping %s" % expnames_long[ex])
        continue
    plotvar = varcheck(plotvar,vname,expnames[ex])
    pcm = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=projd)
    
    cb = viz.hcbar(pcm,ax=ax)
    ax.set_title("%s (%s), %s" % (vname,vunit,expnames_long[ex]))
    plt.show()
    
    
#%% Do a draft plot

#cints    = np.arange(18,35.5,.5)
cints     = np.arange(-0.14,0.15,0.01)#np.arange(120,290,10)  #np.arange(-100,110,10)#np.arange(0,0.95,0.05)
fig,axs   = init_tp_map(2,3,figsize=(24,4.5))

plotorder = [0,2,4,1,3]
cmaptest  = 'cmo.curl'#'cmo.solar'#"plasma_r"

for a in range(nexps):
    print(a)
    
    ax = axs.flatten()[a]

    ex = plotorder[a]
    
    # #
    # ex = 0
    # #
    
    if a < 3:
        #ax.set_xticks([])
        gl = ax.gridlines()
        gl.bottom_labels = False#gl.right_labels = False
        #ax.xaxis.set_major_locator(mpl.ticker.NullLocator())
    
    if a not in [0,3]:
        gl = ax.gridlines()
        gl.left_labels=False
        #ax.set_yticks([])
        #ax.yaxis.set_major_locator(mpl.ticker.NullLocator())
    
    plotvar = ds_var[ex][0][vname]
    plotvar = varcheck(plotvar,vname,expnames[ex])
    
    pcm      = ax.contourf(plotvar.lon,plotvar.lat,plotvar,transform=projd,
                      levels=cints,cmap=cmaptest)
    
    cl = ax.clabel(pcm,levels=cints[:-4][::2],fontsize=8)
    viz.add_fontborder(cl,w=1.5)
    
    ax.set_title(expnames_long[ex],)
    
    ax.tick_params(labelsize=8)

axs[1,2].set_visible(False)
cb = fig.colorbar(pcm,ax=axs[1,2],orientation='horizontal',pad=-5,fraction=0.45)
cb.set_label("Mean %s [%s]" % (indict['longname'],indict['units']))
plt.show()

#%% Plot

cints     = indict['cints_mean']
cmap_in   = indict['cmap_mean']
plotorder = [0,2,4,1,3]

#fig,axs   = init_tp_map(2,3,figsize=(58,4))



bbplot = [120, 290, -20, 20]
fix_lon = np.hstack([np.arange(120,190,10),np.arange(-180,-60,10)])
proj   = ccrs.PlateCarree(central_longitude=180)
projd  = ccrs.PlateCarree()
fig,axs = plt.subplots(2,3,figsize=(18,3.5),subplot_kw={'projection':proj},constrained_layout=True)
for a in range(nexps):
    print(a)
    
    ax = axs.flatten()[a]
    ex = plotorder[a]
    
    # Handle Labeling
    blb = viz.init_blabels()
    if a >= 3:
        blb['lower'] = 1
    if a in [0,3]:
        blb['left'] = 1
        
    
    # Set up Plot
    ax.set_extent(bbplot)
    ax     = viz.add_coast_grid(ax,bbox=bbplot,fill_color='k',blabels=blb,
                                proj=ccrs.PlateCarree(),fix_lon=fix_lon,ignore_error=True)
    
    
    plotvar = ds_var[ex][0][vname]
    plotvar = varcheck(plotvar,vname,expnames[ex])
    
    pcm      = ax.contourf(plotvar.lon,plotvar.lat,plotvar,transform=projd,
                      levels=cints,cmap=cmap_in)
    
    cl = ax.clabel(pcm,levels=cints[:-4][::2],fontsize=8)
    viz.add_fontborder(cl,w=1.5)
    ax.set_title(expnames_long[ex],)
    ax.tick_params(labelsize=8)


axs[1,2].set_visible(False)
cb = fig.colorbar(pcm,ax=axs[1,2],orientation='horizontal',pad=-.5,fraction=0.25)
cb.set_label("Mean %s [%s]" % (indict['longname'],indict['units']))

savename = "%sMeanState_Plot_AWI_%s_Combined.png" % (figpath,vname)
plt.savefig(savename,dpi=150,)#bbox_inches='tight')
plt.show()


