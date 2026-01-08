#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Created on Thu Oct  9 11:52:45 2025

Copied Visualize_Mean_Patterns_TP


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

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

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

tysurdict = dict(
    vname      = "ty_sur",
    longname   = "Meridional Wind Stress to Ocean",
    units      = r"$\frac{m}{s^2}$",
    cints_mean = np.arange(-0.08,0.085,0.005),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.curl",
    )


lccdict = dict(
    vname       = "lcc",
    longname    = "low cloud cover",
    units       ="fraction",
    cints_mean  = np.arange(0,1.1,0.1),
    cmap_mean   = 'cmo.ice',
    )

tccdict = dict(
    vname       = "tcc",
    longname    = "total cloud cover",
    units       ="fraction",
    cints_mean  = np.arange(0,1.1,0.1),
    cmap_mean   = 'cmo.ice',
    )

ttrdict = dict(
    vname      = "ttr",
    longname   = "TOA Net Thermal Radiation",
    units      = r"$\frac{J}{m^2}$",
    cints_mean = np.arange(-300,-195,5),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "plasma_r",
    )

# Just copy to make the clear sky version
ttrcdict = ttrdict.copy()
ttrcdict['vname'] = 'ttrc'
ttrcdict['longname'] = "TOA Clear-Sky Thermal Radiation",



tsrdict = dict(
    vname      = "tsr",
    longname   = "TOA Net Shortwave Radiation",
    units      = r"$\frac{J}{m^2}$",
    cints_mean = np.arange(220,370,10),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.solar",
    )

tsrcdict = dict(
    vname      = "tsrc",
    longname   = "TOA Clear-Sky Shortwave Radiation",
    units      = r"$\frac{J}{m^2}$",
    cints_mean = np.arange(350,384,2),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.solar",
    )

cpdict = dict(
    vname      = "cp",
    longname   = "Convective Precipitation",
    units      = r"$\frac{mm}{day}$",
    cints_mean = np.arange(0,11.5,0.5),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.rain",
    )

lspdict = dict(
    vname      = "lsp",
    longname   = "Stratiform Precipitation",
    units      = r"$\frac{mm}{day}$",
    cints_mean = np.arange(0,11.5,0.5),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.rain",
    )

sshfdict = dict(
    vname      = "sshf",
    longname   = "Surface sensible heat flux",
    units      = r"$\frac{J}{m^2}$",
    cints_mean = np.arange(-100,-22.5,2.5),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.thermal",
    )


slhfdict = dict(
    vname      = "slhf",
    longname   = "Surface latent heat flux",
    units      = r"$\frac{J}{m^2}$",
    cints_mean = np.arange(-200,-45,5),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.thermal",
    )

dicts_byvar = [sstdict,
               sshdict,
               strdict,
               ssrdict,
               txsurdict,
               tysurdict,
               lccdict,
               tccdict,
               ttrdict,
               ttrcdict,
               tsrdict,
               cpdict,
               lspdict,
               slhfdict,
               sshfdict
               ]

vnames_all       = [d["vname"] for d in dicts_byvar]
vdicts           = dict(zip(vnames_all,dicts_byvar))



#%% Additional User Edits

datpath     = "/home/niu4/gliu8/projects/scrap/TP_crop/"
outpath     = "/home/niu4/gliu8/projects/scrap/TP_crop/summary/"
calcnames   = ["mean","scycle","var","monvar"]

# Other Things
projd       = ccrs.PlateCarree()
nexps       = len(expnames)
ncalc       = len(calcnames)
plotorder = [0,2,4,1,3]

#%% Look for monthly plot

for vname in tqdm.tqdm(vnames_all):
    print(vname)
    indict = vdicts[vname]
    
    # Load Metrics
    ds_var = []
    for ex in tqdm.tqdm(range(nexps)):
        
        metrics = []
        for nn in range(ncalc):

            try:
                ds = ut.awi_mean_loader(expnames[ex],vname,calcnames[nn],outpath=outpath)
            except:
                if nn == 0:
                    print("Could not find for %s" % (expnames_long[ex]))
                metrics.append(None)
                
            
            metrics.append(ds)
        ds_var.append(metrics)
        
        
    cints     = indict['cints_mean']
    cmap_in   = indict['cmap_mean']
    cintlab   = 4
    bbplot    = [120, 290, -20, 20]
    fix_lon = np.hstack([np.arange(120,190,10),np.arange(-180,-60,10)])
    proj   = ccrs.PlateCarree(central_longitude=180)
    projd  = ccrs.PlateCarree()
    
    
    for im in range(12):
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
            
            try:
                plotvar = ds_var[ex][1][vname].isel(month=im)
            except:
                print("Skipping %s" % expnames_long[ex])
                ax.set_visible(False)
                continue
            plotvar = ut.varcheck(plotvar,vname,expnames[ex])
            
            pcm      = ax.contourf(plotvar.lon,plotvar.lat,plotvar,transform=projd,
                              levels=cints,cmap=cmap_in)
            
            cl = ax.clabel(pcm,levels=cints[:-4][::cintlab],fontsize=8)
            viz.add_fontborder(cl,w=1.5)
            ax.set_title(expnames_long[ex],)
            ax.tick_params(labelsize=8)
        
        
        axs[1,2].set_visible(False)
        cb = fig.colorbar(pcm,ax=axs[1,2],orientation='horizontal',pad=-.5,fraction=0.25)
        cb.set_label("Mean %s [%s]" % (indict['longname'],indict['units']))
        
        savename = "%sMonMeanState_Plot_AWI_%s_Combined_mon%02i.png" % (figpath,vname,im+1)
        plt.savefig(savename,dpi=150,)#bbox_inches='tight')
        plt.close()
    

