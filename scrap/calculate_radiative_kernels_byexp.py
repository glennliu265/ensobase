#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Copied [calculate_radiative_kernels.py] and streamed with consistent directory structure
 


Created on Fri Nov 28 11:49:47 2025

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

from sklearn.linear_model import LinearRegression
import sklearn

#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% Other Functions

landmask = ut.load_land_mask_awi("TCo319",regrid=True)

# =================
#%% User Selections
# =================

# Path to Data and Experiments


# CERES-FBCT with ERA5 ----- (calculate_unobscured_fluxes_CERES_FBCT)
expname         = "CERES_FBCT_ERA5"
datpath         = "/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_FBCT_ERA5/anom_detrend1/" #/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/"#"/home/niu4/gliu8/projects/scrap/regrid_1x1/"
outpath         = "/home/niu4/gliu8/projects/ccfs/kernels/regrid_1x1/"#%s/" % expname #/home/niu4/gliu8/projects/ccfs/kernels/regrid_1x1/"
anomalize       = False # Kept for legacy. Input should be anomalized before using `anom_detrend1' shellscripts
tstart          = '2002-07-01'
tend            = '2023-02-01'


# # CERES-EBAF with ERA5
# expname         = "CERES_EBAF_ERA5"
# datpath         = "/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_EBAF_ERA5/anom_detrend1/" #/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/"#"/home/niu4/gliu8/projects/scrap/regrid_1x1/"
# outpath         = "/home/niu4/gliu8/projects/ccfs/kernels/regrid_1x1/"#%s/" % expname #/home/niu4/gliu8/projects/ccfs/kernels/regrid_1x1/"
# anomalize       = False # Kept for legacy. Input should be anomalized before using `anom_detrend1' shellscripts
# tstart          = '2000-03-01'
# tend            = '2024-12-31'

# # # ERA5 1979-2024
# # expname         = "ERA5_1979_2024"
# # datpath         = "/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/%s/anom_detrend1/" % expname#/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/"#"/home/niu4/gliu8/projects/scrap/regrid_1x1/"
# # outpath         = "/home/niu4/gliu8/projects/ccfs/kernels/regrid_1x1/"#%s/" % expname #/home/niu4/gliu8/projects/ccfs/kernels/regrid_1x1/"
# # anomalize       = False # Kept for legacy. Input should be anomalized before using `anom_detrend1' shellscripts
# # tstart          = '1979-01-01'
# # tend            = '2024-12-31'

# # ERAS5 2000-2024
# expname         = "ERA5_2000_2024"
# datpath         = "/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/ERA5_2000_2024/anom_detrend1/" #/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/"#"/home/niu4/gliu8/projects/scrap/regrid_1x1/"
# outpath         = "/home/niu4/gliu8/projects/ccfs/kernels/regrid_1x1/"#%s/" % expname #/home/niu4/gliu8/projects/ccfs/kernels/regrid_1x1/"
# anomalize       = False # Kept for legacy. Input should be anomalized before using `anom_detrend1' shellscripts
# tstart          = '2000-03-01'
# tend            = '2024-12-31'

# # 9km 1950 simulation
# expname         = "TCo1279-DART-1950"
# datpath         = "/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/%s/anom_detrend1/" % expname#/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/"#"/home/niu4/gliu8/projects/scrap/regrid_1x1/"
# outpath         = "/home/niu4/gliu8/projects/ccfs/kernels/regrid_1x1/"#%s/" % expname #/home/niu4/gliu8/projects/ccfs/kernels/regrid_1x1/"
# anomalize       = False # Kept for legacy. Input should be anomalized before using `anom_detrend1' shellscripts
# tstart          = None#'1979-01-01'
# tend            = None#'2024-12-31'

# Variables
flxname         = 'tscretotal'#['allsky','clearsky','cre']  # Loop for fluxes
ccf_vars        = ["sst","eis","Tadv","r700","w700","ws10",]#"ucc"] 

selmons_loop    = [None,[12,1,2],[3,4,5],[6,7,8],[9,10,11]] #[[12,1,2],[3,4,5],[6,7,8],[9,10,11]] # [None,]# # Set to None to do 

# MLR Calculation Options
standardize     = True # Set to True to standardize predictors before MLR
fill_value      = 0    # Replace NaN values with <fill_value>
#add_ucc         = False # Set to True to include upper cloud concentration as a predictor

#%% Now Load each predictor variable to do CCFs calculation
# Note that some were computed in [calc_ccfs_regrid.py]
# Others were preprocessed using remapbil in cdo

"""
Searches for datasets in 
    .../ccfs/input_data/regrid_1x1/<expname>/anom_detrend1/

Outputs computed kernels in
    .../ccfs/kernels/regrid_1x1/<expname>/
"""

# Prepare Output Path
outpath_kernel = outpath + "%s/" % expname
proc.makedir(outpath_kernel)

# Load the CCFs
nccfs = len(ccf_vars)
dsvars_anom = []
for v in range(nccfs):
    
    ncname = "%s%s.nc" % (datpath,ccf_vars[v])
    ds     = xr.open_dataset(ncname)[ccf_vars[v]].load()
    
    # Standardize Names
    ds     = ut.standardize_names(ds)
    
    # <Need to add additional fixes/automatic checks>
    if expname == 'TCo2559-DART-1950C':
        if ccf_vars[v] == "w700": # Duplicate a month for the missing variables
            w700data = ds.data.squeeze()
            w700data_duplicate_jan1950 = np.concatenate([w700data[[0],...],w700data],axis=0)
            newtime = dsvars_anom[v-1].time
            coords  = dict(time=newtime,lat=ds.lat,lon=ds.lon)
            w700new = xr.DataArray(w700data_duplicate_jan1950,coords=coords,dims=coords,name='w700')
            ds = w700new
        
        print("%s, %s" % (expname,ccf_vars[v]))
        print(ds.shape)
        print("")
        dsvars_anom.append(ds.squeeze())
    else:
        dsvars_anom.append(ds.squeeze())
    
    # SST is only 8 years, so reduce the time...
    if expname == 'TCo2559-DART-1950C': 
        dsvars_anom = [ut.reduce_time(ds,dsvars_anom[0]) for ds in dsvars_anom]

# Load the fluxes
dsflx   = xr.open_dataset("%s%s.nc" % (datpath,flxname))[flxname].load()
if "TCo" in expname: # Convert Fluxes
    print("Converting fluxes based on accumulation period for AWI-CM3")
    dsflx   = ut.varcheck(dsflx,flxname,expname)
elif expname == "ERA5_1979_2024":
    print("Considering 1-day flux accumulation period for ERA_1979_2024")
    dtday   = 3600*24
    dsflx   = dsflx / dtday
dsflx   = ut.standardize_names(dsflx)

# Shift CERES data (move function to ut)
import pandas as pd
def shift_time_monthstart(dsin,timename='time'):
    oldtime         = dsin[timename]
    tstart          =  str(oldtime[0].data)[:7] + "-01"
    tend            =  str(oldtime[-1].data)[:7] + "-01"
    newtime         = pd.date_range(start=tstart,end=tend,freq="MS")
    dsin[timename]    = newtime
    print("New time dimension between %s and %s" % (dsin[timename][0].data,dsin[timename][-1].data))
    return dsin
if "CERES" in expname:
    dsflx = ut.shift_time_monthstart(dsflx,timename='time')
    
# Limit to time
if tstart or tend is not None:
    dsflx       = dsflx.sel(time=slice(tstart,tend))
    dsvars_anom = [ds.sel(time=slice(tstart,tend)) for ds in dsvars_anom]
else:
    print("No time range will be applied")
    
# for ii in range(nccfs): (note match time month does not play nice with pd.date_range, need to figure it out)
#     dsvars_anom[ii],_=proc.match_time_month(dsvars_anom[ii],dsflx)

# Check Time Dimension
ntimes_predictors = [len(ds.time) for ds in dsvars_anom]
ntimes_flux       = len(dsflx.time)

# =====================================
#%% Part (2): Compute Radiative Kernels
# =====================================

st = time.time()

# Subset months
for selmons in selmons_loop:
    
    dsexp_sel  = dsvars_anom
    dsexp_flx  = dsflx
    
    if selmons is not None:
        dsexp_flx = proc.selmon_ds(dsexp_flx,selmons)
        dsexp_sel = [proc.selmon_ds(ds,selmons) for ds in dsexp_sel]
    else:
        print("Calculating for all months!")
    
    # Pre-allocate
    dsexp_flx       = dsexp_flx.transpose('lat','lon','time')
    lon             = dsexp_flx.lon.data
    lat             = dsexp_flx.lat.data
    nlat,nlon,ntime = dsexp_flx.shape
    nccfs           = len(dsexp_sel)
    coeffs          = np.zeros((nlat,nlon,nccfs)) * np.nan # [ Lat x Lon x CCFs ]
    ypred           = np.zeros(dsexp_flx.shape) * np.nan   # [ Lat x Lon x Time ]
    r2              = np.zeros((nlat,nlon)) * np.nan       # [ Lat x Lon ]
    yerr            = ypred.copy()
    
    # Do a silly loop (took 5 min 17 sec)
    for o in tqdm.tqdm(range(nlon)):
        lonf = lon[o]
        
        for a in range(nlat):
            latf = lat[a]
            
            chkland = proc.selpt_ds(landmask,lonf,latf).data
            if np.isnan(chkland):
                continue
            
            # Check for NaN in predictor
            dspts  = [proc.selpt_ds(ds,lonf,latf) for ds in dsexp_sel]
            chknan = [np.any(np.isnan(ds.data)) for ds in dspts]
            if np.any(chknan):
                iinan = np.where(chknan)[0][0]
                #print("NaN detected for variables %s, lon (%.2f), lat (%.2f)... skipping." % (chknan,lonf,latf))
                continue
            # Check for NaN in target
            flxpt  = proc.selpt_ds(dsexp_flx,lonf,latf)
            if np.any(np.isnan(flxpt.data)):
                #print("NaN detected for Flux, lon (%.2f), lat (%.2f)... skipping." % (lonf,latf))
                continue
            
            # Do calculations
            mlr_out = ut.mlr_ccfs(dspts,flxpt,standardize=standardize,verbose=False)
            
            # Calculate error
            yerr_pt = flxpt - mlr_out['pred']
            
            r2[a,o]       = mlr_out['r2']
            ypred[a,o,:]  = mlr_out['pred']
            coeffs[a,o,:] = mlr_out['coeffs']
            yerr[a,o,:]   = yerr_pt
            
            # if add_ucc:
            #     ccfnames = ccf_vars
            
            # else:
            #     if "ucc" in ccf_vars:
            #         ccfnames = ccf_vars[:-1]
            #     ccfnames = ccf_vars
            
    ccfnames       = ccf_vars
    
    coords_r2       = dict(lat=lat,lon=lon)
    coords_coeffs   = dict(lat=lat,lon=lon,ccf=ccfnames)
    coords_pred     = dict(lat=lat,lon=lon,time=dsexp_flx.time)
    
    da_r2           = xr.DataArray(r2,coords=coords_r2,dims=coords_r2,name='r2')
    da_coeffs       = xr.DataArray(coeffs,coords=coords_coeffs,dims=coords_coeffs,name='coeffs')
    da_pred         = xr.DataArray(ypred,coords=coords_pred,dims=coords_pred,name='ypred')
    da_yerr         = xr.DataArray(yerr,coords=coords_pred,dims=coords_pred,name='yerr')
    ds_out          = xr.merge([da_r2,da_coeffs,da_pred,da_yerr])
    edict           = proc.make_encoding_dict(ds_out)
    
    outname         = "%s%s_kernels_standardize%i.nc" % (outpath_kernel,flxname,standardize)
    if selmons is not None:
        selmonstr    = proc.mon2str(np.array(selmons)-1)
        outname      = proc.addstrtoext(outname,"_"+selmonstr,adjust=-1)
    
    ds_out.to_netcdf(outname,encoding=edict)
    print("Completed CCF kernel calculation for %s (%s) in %.2fs" % (flxname,expname,time.time()-st))

    