#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 11:08:50 2025

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
import importlib

#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% Load Land Mask and set other paths

figpath     = "/home/niu4/gliu8/figures/bydate/2026-02-11/"
landmask    = ut.load_land_mask_awi("ERA5",regrid=True)

#%% Load CCFs

# Kernel Information
expname      = "CERES_EBAF_ERA5_2001_2024"
datpath      = "/home/niu4/gliu8/projects/ccfs/regrid_1x1/"
ccf_vars     = ["sst","eis","Tadv","r700","w700","ws10"]
flxname      = "cre"
standardize  = True
add_ucc      = False
selmons_loop = [[12,1,2],[3,4,5],[6,7,8],[9,10,11]] # [[]]#.None
regrid1x1    = True

# Update
#outpath     = "/home/niu4/gliu8/projects/ccfs/regrid_1x1/%s_components/" % flxname

figpath      = "/home/niu4/gliu8/figures/bydate/2026-02-11/"
proc.makedir(figpath)

#% Make Experiment Directories
radpath      = "/home/niu4/gliu8/projects/ccfs/radiative_components/regrid_1x1/%s/" % expname
outpath      = "/home/niu4/gliu8/projects/ccfs/enso_regression_patterns/by_ccf/regrid_1x1/%s/" % expname
proc.makedir(outpath)

#%% Additional Variables

vnames_new = [flxname + "_" + ccfname for ccfname in ccf_vars]

# =========================
#%% Load all months
# =========================

nccfs       = len(ccf_vars)
dsall       = []
dsallmon    = []

for cc in tqdm.tqdm(range(nccfs)):
    
    vname_new      = vnames_new[cc]
    ccf_name       = ccf_vars[cc]
    
    # Load All Months
    rname_out      = "%s%s_component.nc" % (radpath,vname_new,)
    ds             = xr.open_dataset(rname_out)[ccf_name].load()
    dsall.append(ds)
    
    # Load Seasonal Calculations
    try:
        rname_out      = "%s%s_component_seasonal.nc" % (radpath,vname_new,)
        ds = xr.open_dataset(rname_out)[ccf_name].load()
        dsallmon.append(ds)
    except:
        rname_out      = "%s%s_component_seasonal.nc" % (radpath,vname_new,)
        ds = xr.open_dataset(rname_out).load()
        ds = ds.rename({'__xarray_dataarray_variable__': vname_new})
        print("Could not find variable name in %s" % vname_new)

# =========================
#%% Plot a single timestep (debugging)
# =========================

# itime = 5
# proj = ccrs.PlateCarree()


# fig,axs = ut.init_globalmap(nrow=3,ncol=2,figsize=(12,14))

# for cc in range(nccfs):
    
#     ax      = axs.flatten()[cc]
#     plotvar = dsall[cc].isel(valid_time=itime).squeeze()
    
#     title  = r"%s" % (plotvar.name)
#     ax.set_title(title)
#     pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,
#                         cmap='cmo.balance',
#                         )

#     cb = viz.hcbar(pcm,ax=ax,pad=0.01,fraction=0.015)

# plt.show()


# =====================
#%% Load EOF-based ENSO
# =====================

#ensoeof_path = "/home/niu4/gliu8/projects/ccfs/enso_eof/"
#ncname       = "%sERA5_1979_2024_ENSO_EOF.nc" % (ensoeof_path,)
#dsenso       = xr.open_dataset(ncname).load()

ensopath    = "/home/niu4/gliu8/projects/scrap/nino34/"
if "ERA5" in expname:
    print("Loading Rotated EOF for ERA5...")
    expname_nino = "ERA5_1979_2024"
else:
    expname_nino = expname
ncname      = "%s%s_enso_eof_rotated.nc" % (ensopath,expname_nino)
dsenso      = xr.open_dataset(ncname).load()

cp          = dsenso.cp
ep          = dsenso.ep

# Check to make sure they are the same length as the inputs
cp,_ = proc.match_time_month(cp,dsall[0])
ep,_ = proc.match_time_month(ep,dsall[0])

X           = np.array([ep,cp])
X_std       = np.array([ep/np.std(ep), cp/np.std(cp)])

# =========================
#%% Perform MLR fit to ENSO
# =========================
standardize_nino = False # Set to True to Standardize ENSO Index

if standardize_nino:
    X_in = X_std
else:
    X_in = X

for ss in range(2):
    
    if ss == 0: # Do for all months
        dsin   = dsall
        monstr = "" 
    else: # Do for seasonal estimates
        dsin   = dsallmon
        monstr = "_seasonal"
    
    for cc in range(nccfs):
        rin         = dsin[cc].squeeze()
        rin         = ut.standardize_names(rin)
        rin         = rin.sortby('time') # Temp Fix (will add this to ccf_radiation computation code)
        rin,_       = proc.match_time_month(rin,ep)
        
        ccfname     = ccf_vars[cc]
        rin         = rin.transpose('lat','lon','time')
        nlat,nlon,ntime = rin.shape
        
        beta_cp     = np.zeros((nlat,nlon)) * np.nan
        beta_ep     = beta_cp.copy()
        r2_ccf      = beta_cp.copy()
        ypred       = np.zeros((nlat,nlon,ntime)) * np.nan
        yerr        = ypred.copy()
        for a in tqdm.tqdm(range(nlat)):
            
            for o in range(nlon):
                
                rpt = rin.isel(lat=a,lon=o)#.data[a,o,:]
                
                if np.any(np.isnan(rpt)):
                    continue
                
                mlr_out = ut.mlr_ccfs(X_in,rpt,standardize=False,fill_value=0)
                
                beta_ep[a,o] = mlr_out['coeffs'][0]
                beta_cp[a,o] = mlr_out['coeffs'][1]
                r2_ccf[a,o]  = mlr_out['r2'] 
                ypred[a,o,:] = mlr_out['pred']
                yerr[a,o,:]  = mlr_out['err']
         
        # %
        lat = rin.lat
        lon = rin.lon
        times = rin.time
                
        coords_latlon   = dict(lat=lat,lon=lon)
        coords_pred     = dict(lat=lat,lon=lon,time=times)
        da_r2           = xr.DataArray(r2_ccf,coords=coords_latlon,dims=coords_latlon,name='r2')
        da_ep           = xr.DataArray(beta_ep,coords=coords_latlon,dims=coords_latlon,name='beta_ep')
        da_cp           = xr.DataArray(beta_cp,coords=coords_latlon,dims=coords_latlon,name='beta_cp')
        da_ypred        = xr.DataArray(ypred,coords=coords_pred,dims=coords_pred,name='pred')
        da_yerr         = xr.DataArray(yerr,coords=coords_pred,dims=coords_pred,name='err')
        
        ds_out          = xr.merge([da_r2,da_cp,da_ep,da_ypred,da_yerr])
        edict           = proc.make_encoding_dict(ds_out)
        #outname         = "%s%s_%s_%s_ENSO_regression_stdnino%i.nc" % (outpath,expname,flxname,ccfname,standardize_nino)
        outname         = "%s%s_ccf_%s_ENSO_regression_stdnino%i%s.nc" % (outpath,flxname,ccfname,standardize_nino,monstr)
        ds_out.to_netcdf(outname,encoding=edict)
    
#%% Visualize CCF radiation (move to another script)
# Note: Moved to python notebook

# ds_ensoccf = []

# for cc in range(nccfs):
#     ccfname = ccf_vars[cc]
    
#     outname = "%s%s_%s_%s_ENSO_regression_stdnino%i.nc" % (outpath,expname,flxname,ccfname,standardize_nino)
#     # if standardize_nino:
#     #     outname = proc.addstrtoext(outname,"_standardizeNino",adjust=-1)
#     ds = xr.open_dataset(outname)#[ccfname].load()
#     ds_ensoccf.append(ds)
    
# #%% Plot Each CCF

# proj    = ccrs.PlateCarree()
# fig,axs = ut.init_globalmap(nrow=1,ncol=6,figsize=(24,10))

# for cc in range(nccfs):
    
#     ax      = axs.flatten()[cc]
#     plotvar = ds_ensoccf[cc].beta_cp#dsall[cc].isel(valid_time=itime).squeeze()
    
#     title   = r"%s" % (ccf_vars[cc])
#     ax.set_title(title)
#     pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,
#                         cmap='cmo.balance',
#                         )
    
#     cb = viz.hcbar(pcm,ax=ax,pad=0.01,fraction=0.015)
# #cb.set_label(r"%s Feedback [W $m^{-2}$ $\sigma^{-1}$]" % flxname)
# plt.show()