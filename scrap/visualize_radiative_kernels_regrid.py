#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Visualize Radiative Kernels computed with [calc_ccfs_regrid]

Created on Mon Nov 10 16:39:39 2025

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
import climlab

#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% Load Land Mask and set other paths

figpath     = "/home/niu4/gliu8/figures/bydate/2025-11-12/"
landmask    = ut.load_land_mask_awi("TCo319",regrid=True)

#%% Load Output

expnames    = ["TCo1279-DART-1950","TCo2559-DART-1950C"]
datpath     = "/home/niu4/gliu8/projects/scrap/regrid_1x1/ccfs_regression_global/"
ccf_vars    = ["sst","eis","Tadv","r700","w700","WS","ucc"] 
rawpath     = "/home/niu4/gliu8/projects/scrap/regrid_1x1/"
standardize = True
add_ucc     = True
flxname     = "cre"

ds_flxs     = []
ds_all      = []
for ex in range(2):
    
    
    outname         = "%s%s_%s_CCFs_Regression_standardize%i_regrid1x1_adducc%i.nc" % (datpath,expnames[ex],flxname,standardize,add_ucc)
    ds = xr.open_dataset(outname).load()
    ds_all.append(ds)
    
    
    flxnc   = "%s%s_%s_regrid1x1.nc" % (rawpath,expnames[ex],flxname)
    ds      = xr.open_dataset(flxnc).load()
    
    dsflx   = ut.preprocess_enso(ut.standardize_names(ds[flxname]))
    dsflx   = ut.varcheck(dsflx,flxname,expnames[ex])
    
    ds_flxs.append(dsflx)
    


#%% Lets Visualize R2

fig,axs      = ut.init_globalmap(1,2,figsize=(12,4.5))

proj = ccrs.PlateCarree()

for ex in range(2):
    ax          = axs[ex]
    plotvar     = ds_all[ex].r2
    pcm         = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,cmap='cmo.deep',vmin=0,vmax=1)
    ax.set_title(expnames[ex])

cb              = viz.hcbar(pcm,ax=axs.flatten(),fraction=0.045,pad=0.01)
savename        = "%s%s_CCFs_Regression_standardize%i_regrid1x1_adduccc%i_R2.png" % (figpath,flxname,standardize,add_ucc)

plt.savefig(savename,dpi=150,bbox_inches='tight')

#%% Zoom in on region and select a point
# Note that this is more of a debugging one...

fig,axs      = ut.init_tp_map(1,2,figsize=(12,4.5)) #plt.subplots(1,2,constrained_layout=True,figsize=(12,4.5),subplot_kw={'projection':proj})#ut.init_globalmap(1,2,figsize=(12,4.5))
proj         = ccrs.PlateCarree()

for ex in range(2):
    ax          = axs[ex]
    plotvar     = ds_all[ex].r2
    pcm         = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,cmap='cmo.deep',vmin=0,vmax=1)
    ax.set_title(expnames[ex])
    ax.set_extent([120,-70,-30,30])
    ax.plot(360-86,-5,marker="x",markersize=25,transform=proj,c='k')
    ax.plot(360-135,0,marker="x",markersize=25,transform=proj,c='red')
    #ax.set_extent([0,360,-30,30])
    #ax.set_extent([-180,180,-30,30])

cb           = viz.hcbar(pcm,ax=axs.flatten(),fraction=0.045,pad=0.05)

#plt.show()
# savename = "%s%s_%s_CCFs_Regression_standardize%i_regrid1x1_adduccc%i_R2.png" % (figpath,expnames[ex],flxname,standardize,add_ucc)
# plt.savefig(savename,dpi=150,bbox_inches='tight')
# plt.show()

#%% Check Prediction/Timeseries at a point

# lonf = 360-86
# latf = -5

for ii in range(2):
    if ii == 0: # EEP
        lonf = 360-135#360-86
        latf = 0#-5
    else: # SEP
        lonf = 360-86
        latf = -5
        
    
    locfn,loctitle = proc.make_locstring(lonf,latf)
    fig,axs = plt.subplots(2,1,constrained_layout=True,figsize=(12.5,6))
    
    
    for ex in range(2):
        
        ax    = axs[ex]
        dspt  = proc.selpt_ds(ds_all[ex],lonf,latf,)
        flxpt = proc.selpt_ds(ds_flxs[ex],lonf,latf,)
        
        # Plot Target
        plotvar = flxpt
        ax.plot(plotvar.time,plotvar,color="k",label="%s (Raw)" % flxname)
        
        # Plot Prediction
        plotvar = dspt.ypred
        ax.plot(plotvar.time,plotvar,color="limegreen",label="MLR Fit")
        
        if ex == 0:
            ax.legend(ncol=2)
        
        # Set Title
        title  = r"%s: $R^2$=%.2f" % (expnames[ex],dspt.r2.data.item())
        title2 = [r"$\beta_{%s}$=%.2f" % (ccf_vars[ii],dspt.coeffs.data[ii]) for ii in range(len(dspt.ccf.data))]
        title2 = ", ".join(title2)
        ax.set_title(title+"\n"+title2)
        
        ax.set_xlim([plotvar.time.isel(time=0),plotvar.time.isel(time=-1)])
        
        ax.axhline([0],c='k',ls='dashed',lw=0.75)
        ax.set_ylabel(r"%s [$\frac{W}{m^2}$]" % flxname)
    plt.suptitle("%s MLR Fit @ %s" % (flxname,loctitle))
        
        
    
    savename = "%s%s_CCFs_Regression_Fit_adducc%i_%s.png" % (figpath,flxname,add_ucc,locfn)
    plt.savefig(savename,dpi=150,bbox_inches='tight')
    #plt.show()
    
    #% Check the Shape of the Residual for selected points
    

    
    fig,axs = plt.subplots(2,1,constrained_layout=True,figsize=(6,6),sharex=True)
    
    for ex in range(2):
        
        ax    = axs[ex]
        dspt  = proc.selpt_ds(ds_all[ex],lonf,latf,)
        flxpt = proc.selpt_ds(ds_flxs[ex],lonf,latf,)
        
        diff  = flxpt - dspt.ypred
        
        ax.hist(diff,bins=15,edgecolor="k")
        ax.set_title(expnames[ex])
        ax.set_xlabel("Residual [W/m2]")
    plt.suptitle("%s MLR Fit @ %s" % (flxname,loctitle))
    
    savename = "%s%s_CCFs_Residual Distribution_adducc%i_%s.png" % (figpath,flxname,add_ucc,locfn)
    plt.savefig(savename,dpi=150,bbox_inches='tight')
    #plt.show()

#%% Visualize the Coefficients


for ex in range(2):
    fig,axs = ut.init_globalmap(nrow=3,ncol=2,figsize=(12,14))
    
    nccfs = len(ds_all[ex].ccf)
    if add_ucc:
        nccfs = nccfs-1
    
    
    for cc in range(nccfs):
        
        ax      = axs.flatten()[cc]
        
        plotvar = ds_all[ex].coeffs.isel(ccf=cc) #/ dtmon
        title  = r"d(%s) / d (%s)" % (flxname,plotvar.ccf.item())
        ax.set_title(title)
        pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,
                            vmin=-5,vmax=5,cmap='cmo.balance',
                            )
    
    plt.suptitle("%s (%s) MLR Coefficients" % (flxname,expnames[ex]),y=0.95)
    cb = viz.hcbar(pcm,ax=axs.flatten(),pad=0.01,fraction=0.015)
    cb.set_label(r"%s Feedback [W $m^{-2}$ $\sigma^{-1}$]" % flxname)
    savename = "%s%s_%s_CCFs_Map_Combine_adducc%i.png" % (figpath,flxname,expnames[ex],add_ucc)
    plt.savefig(savename,dpi=150,bbox_inches='tight')
    #plt.show()

