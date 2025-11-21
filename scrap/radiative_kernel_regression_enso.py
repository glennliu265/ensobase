#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Radiative Kernel Regression ENSO

1. Compute CCF-related change in selected flux.
2. Regress to ENSO
3. Visualize the Pattern

Copied upper section of visualize_radiative_kernels_regrid.py

Created on Wed Nov 12 10:52:28 2025

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

#%% Load Coefficients and Flux

expnames    = ["TCo1279-DART-1950","TCo2559-DART-1950C"]
datpath     = "/home/niu4/gliu8/projects/scrap/regrid_1x1/ccfs_regression_global/"
ccf_vars    = ["sst","eis","Tadv","r700","w700","WS","ucc"] 
rawpath     = "/home/niu4/gliu8/projects/scrap/regrid_1x1/"
standardize = True
add_ucc     = False
flxname     = "cre"
ds_flxs     = []
ds_all      = []
for ex in range(2):
    
    # Load the coefficients
    outname         = "%s%s_%s_CCFs_Regression_standardize%i_regrid1x1_adducc%i.nc" % (datpath,expnames[ex],flxname,standardize,add_ucc)
    ds = xr.open_dataset(outname).load()
    ds_all.append(ds)
    
    # Load the Flux
    flxnc   = "%s%s_%s_regrid1x1.nc" % (rawpath,expnames[ex],flxname)
    ds      = xr.open_dataset(flxnc).load()
    
    # Detrend and Deseason the Flux
    dsflx   = ut.preprocess_enso(ut.standardize_names(ds[flxname]))
    dsflx   = ut.varcheck(dsflx,flxname,expnames[ex])
    ds_flxs.append(dsflx)

#%% Load ENSO indices

# First, Load Nino3.4 Index
standardize = False
ensoid_name = "nino34"
ensoids     = [ut.load_ensoid(expname, ensoid_name,
                          standardize=standardize) for expname in expnames]

#%% Load each predictor variable (copied from calc_ccfs_regrid)


ccf_vars = ["sst","eis","Tadv","r700","w700","WS","ucc"] 
if add_ucc is False:
    ccf_vars = ccf_vars[:-1]
def reduce_time(ds,dsst):
    dsnew,dsst = proc.match_time_month(ds,dsst)
    return dsnew

# Load CCFs
dsbyexp = []
for ex in range(2):
    dsvars = []
    for v in range(len(ccf_vars)):
        
        ncname = "%s%s_%s_regrid1x1.nc" % (rawpath,expnames[ex],ccf_vars[v])
        ds     = xr.open_dataset(ncname)[ccf_vars[v]].load()
        ds     = ut.standardize_names(ds)
        
        if expnames[ex] == 'TCo2559-DART-1950C':
            if ccf_vars[v] == "w700": # Duplicate a month for the missing variables
                w700data = ds.data.squeeze()
                w700data_duplicate_jan1950 = np.concatenate([w700data[[0],...],w700data],axis=0)
                newtime = dsvars[v-1].time
                coords  = dict(time=newtime,lat=ds.lat,lon=ds.lon)
                w700new = xr.DataArray(w700data_duplicate_jan1950,coords=coords,dims=coords,name='w700')
                ds = w700new
        
        print("%s, %s" % (expnames[ex],ccf_vars[v]))
        print(ds.shape)
        print("")
        dsvars.append(ds.squeeze())
    
    # SST is only 8 years, so reduce the time...
    if expnames[ex] == 'TCo2559-DART-1950C': 
        dsvars = [reduce_time(ds,dsvars[0]) for ds in dsvars]
    
    dsbyexp.append(dsvars)

# anomalize and detrend
dsbyexp_anoms = []
for ex in range(2):
    dsvars_anoms = []
    
    for v in tqdm.tqdm(range(len(ccf_vars))):
        
        dsin    = dsbyexp[ex][v]
        dsanoms = ut.preprocess_enso(dsin)
        dsvars_anoms.append(dsanoms)
    
    dsbyexp_anoms.append(dsvars_anoms)
    
#%% Part 1: Reconstruct the Component Fluxes


flxccf_byexp = []
for ex in range(2):
    
    flx_byccf = []
    
    for cc in tqdm.tqdm(range(len(ccf_vars))):
        
        ccfname      = ccf_vars[cc]
        ccf_anom     = dsbyexp_anoms[ex][cc]
        ccf_coeff    = ds_all[ex].coeffs.sel(ccf=ccfname)
        
        flx_ccf_comp = (ccf_coeff * ccf_anom).rename(ccfname)
        
        flx_byccf.append(flx_ccf_comp)
    
    flxccf_byexp.append(flx_byccf)
        

#%% Debug Chunk (Erase Later) -<0>- -- -<0>- -- -<0>- -- -<0>- -- -<0>- 

# Manual Check (to make sure xarray is doing what i expect....)
testarr     = np.stack([ds.data for ds in invar])
test_anom   = ccf_anom.data
test_coeff  = ccf_coeff.data
test_mult   = test_anom * test_coeff[:,:,None]
test_mult == testarr[cc,:,:]
diff = test_mult - testarr[cc,:,:]
# -- -<0>- -- -<0>- -- -<0>- -- -<0>- -- -<0>- -- -<0>- -- -<0>- -- -<0>-

#%% Part 2: Regression to ENSO

ex      = 1
leadlags = np.array([0,])



regpat_byexp = []
for ex in range(2):
    invar   = flxccf_byexp[ex]
    ensoin  = ensoids[ex]
    
    regpats = []
    for cc in range(len(ccf_vars)):
        
        inflx  = invar[cc]
        lagout = ut.calc_leadlag_regression_2d(ensoin,inflx,leadlags)
        regpats.append(lagout)
        
    regpat_byexp.append(regpats)
        
#%% Make a Plot

proj = ccrs.PlateCarree()

for ex in range(2):
    
    fig,axs = ut.init_globalmap(nrow=3,ncol=2,figsize=(12,14))
    
    for cc in range(len(ccf_vars)):
        
        ax      = axs.flatten()[cc]
        ccfname = ccf_vars[cc]
        
        plotvar = regpat_byexp[ex][cc][ccfname].squeeze()
        title  = r"d(%s) / d (%s)" % (flxname,ccfname)
        
        ax.set_title(title)
        pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,
                            vmin=-1,vmax=1,cmap='cmo.balance',
                            )
    
    plt.suptitle("%s (%s) ENSO Regression Coefficients" % (flxname,expnames[ex]),y=0.95)
    cb = viz.hcbar(pcm,ax=axs.flatten(),pad=0.01,fraction=0.015)
    cb.set_label(r"%s Feedback [W $m^{-2}$ $\sigma^{-1}$]" % flxname)
    savename = "%s%s_%s_CCFs_ENSO_Regression_%s_adducc%i.png" % (figpath,flxname,expnames[ex],ensoid_name,add_ucc)
    plt.savefig(savename,dpi=150,bbox_inches='tight')
    
    plt.show()


#%% Checkout Individual CCFs
cc = 2
ex = 0
fig,ax  = ut.init_globalmap(nrow=1,ncol=1,figsize=(12,4.5))
ccfname = ccf_vars[cc]

# For Tadv
regpat_byexp[ex][cc][ccfname].plot(ax=ax,transform=proj,vmin=-5e-6,vmax=5e-6,cmap='cmo.balance')
#regpat_byexp[ex][cc][ccfname].plot(ax=ax,transform=proj,vmin=-1,vmax=1,cmap='cmo.balance')


title  = r"%s (%s) Component, %s" % (flxname,ccfname, expnames[ex])

cb     = viz.hcbar(pcm,ax=axs.flatten(),pad=0.01,fraction=0.015)
cb.set_label(r"%s Feedback [W $m^{-2}$ $\sigma^{-1}$]" % flxname)
ax.set_title(title)

plt.show()



#%% 

lagout.eis.squeeze().plot(vmin=-5,vmax=5,cmap='cmo.balance'),plt.show()

#%% Debug: Check EIS Differences

eis_diff = ds_all[1].coeffs.sel(ccf='eis') - ds_all[0].coeffs.sel(ccf='eis')



eis_diff.plot(),plt.plot(275,0,marker="x",color="k",markersize=24),plt.show()




#%% Make Scatterplot

lonf = 275
latf = 0

keis = ccf_vars.index('eis')

fig,axs = plt.subplots(1,2,figsize=(8,4.5),constrained_layout=True)

for ex in range(2):
    
    ax = axs[ex]
    plotx = proc.selpt_ds(ds_flxs[ex],lonf,latf)
    ploty = proc.selpt_ds(dsbyexp_anoms[ex][keis],lonf,latf)
    #ploty = proc.selpt_ds(flxccf_byexp[ex][keis],lonf,latf)
    plotx,ploty = proc.match_time_month(plotx,ploty)
    
    coeff = proc.selpt_ds(ds_all[ex]['coeffs'].isel(ccf=keis),lonf,latf).item()
    r2    = proc.selpt_ds(ds_all[ex]['r2'],lonf,latf).item()
    lab   = r"$\beta = %.2f, r^2=%.2f$" % (coeff,r2)
    
    ax.scatter(plotx,ploty,label=lab)
    
    
    ax.set_xlabel("CRE (W/m2)")
    ax.set_ylabel("EIS")
    ax.legend()
plt.show()


#%% Make Lineplot


fig,axs = plt.subplots(2,1,figsize=(12,4.5),constrained_layout=True)

for ex in range(2):
    
    ax = axs[ex]
    plotx = proc.selpt_ds(ds_flxs[ex],lonf,latf)
    ploty = proc.selpt_ds(dsbyexp_anoms[ex][keis],lonf,latf)
    #ploty = proc.selpt_ds(flxccf_byexp[ex][keis],lonf,latf)
    plotx,ploty = proc.match_time_month(plotx,ploty)
    
    
    coeff = proc.selpt_ds(ds_all[ex]['coeffs'].isel(ccf=keis),lonf,latf).item()
    r2    = proc.selpt_ds(ds_all[ex]['r2'],lonf,latf).item()
    lab   = r"$\beta = %.2f, r^2=%.2f$" % (coeff,r2)
    
    
    ax.plot(plotx.time,plotx,label="CRE",color="red")
    ax2 = ax.twinx()
    ax2.plot(plotx.time,ploty,label="EIS",color="blue")
    
    #ax.scatter(plotx,ploty,label=lab)
    
    
    #ax.set_xlabel("CRE (W/m2)")
    #ax.set_ylabel("EIS")
    ax.legend()
plt.show()


#%%