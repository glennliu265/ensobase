#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate the PC-based NAO Index by performing an EOF analysis separately for
anomalies of each month.

Created on Fri Feb 20 14:38:38 2026

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
import matplotlib as mpl

import xeofs as xe

#%% Helper Functions

def lon360to180_xr(ds,lonname='lon'):
    # Based on https://stackoverflow.com/questions/53345442/about-changing-longitude-array-from-0-360-to-180-to-180-with-python-xarray
    ds.coords[lonname] = (ds.coords[lonname] + 180) % 360 - 180
    ds = ds.sortby(ds[lonname])
    return ds

def swap_rename(ds,chkvar,newvar):
    if chkvar in list(ds.coords):
        print("Renaming [%s] to [%s]" % (chkvar,newvar))
        ds = ds.rename({chkvar:newvar})
    return ds

def standardize_names(ds):
    
    ds = swap_rename(ds,'valid_time','time')
    ds = swap_rename(ds,'time_counter','time')
    ds = swap_rename(ds,"TIME_COUNTER",'time')
    ds = swap_rename(ds,"LON","lon")
    ds = swap_rename(ds,"LAT","lat")
    ds = swap_rename(ds,"longitude","lon")
    ds = swap_rename(ds,"latitude","lat")
    ds = swap_rename (ds,"LAT232_409","lat")
    
    # Other preprocessing
    # drop LON_bnds, TIME_COUNTER_bnds
    dropvars = ["LON_bnds","TIME_COUNTER_bnds"]
    for dropvar in dropvars:
        if dropvar in ds:
            ds = ds.drop_vars(dropvar)
    return ds

def xrdeseason(ds,check_mon=True):
    """ Remove seasonal cycle, given an Dataarray with dimension 'time'"""
    if check_mon:
        try: 
            if ds.time[0].values.item().month != 1:
                print("Warning, first month is not Jan...")
        except:
            print("Warning, not checking for feb start")
    
    return ds.groupby('time.month') - ds.groupby('time.month').mean('time')

def xrdetrend(ds,timename='time',verbose=True):
    
    detrendlin = lambda x: sp.signal.detrend(x)
    result = xr.apply_ufunc(
        detrendlin,
        ds,
        input_core_dims=[["time"]],
        output_core_dims=[["time"]],
        vectorize=True,
        dataset_fill_value=np.nan
        )
    result['time'] = ds['time']
    return result

def selmon(ds,mon):
    return ds.sel(time=ds.time.dt.month.isin([mon]))

#%% Input Information

bbox     = [-90,40,20,80] # Assumes with degrees west #[-90+360, 40, 20, 80]

# Information for SLP Data (Input)
vname    = "msl" # Name of the variable
datpath  = "/home/niu4/gliu8/projects/scrap/processed_global/" # Input Path
expname  = "TCo319-DART-ssp585d-gibbs-charn" # Name of the experiment (for naming)
expnames = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-2060","TCo1279-DART-1950","TCo1279-DART-2090_msl"] # "TCo319-DART-ctl1950d-gibbs-charn","TCo95-hi1950d","TCo95-ssp585d",

for expname in expnames:
    
    # Additional Input Information (for expname loop)
    ncname   = "%s_%s.nc" % (expname,vname)# Name of Netcdf
    
    # Info for output
    outpath     = "/home/niu4/gliu8/projects/scrap/nao_indices/"
    ncname_out  = "%s%s_NAO_Indices.nc" % (outpath,expname,)
    ncname_out_allmon  = "%s%s_NAO_Indices_AllMonths.nc" % (outpath,expname,)
    
    
    
    #%% Load the data
    
    # Open a View
    ds      = xr.open_dataset(datpath+ncname)[vname]
    
    # Select the NAO Region
    if np.any(ds.lon) > 180:
        print("Flipping Longitude...")
        ds      = lon360to180_xr(ds) # Correct Longitude (assuming it has values over 360)
    dsreg   = ds.sel(lon=slice(bbox[0],bbox[1]),lat=slice(bbox[2],bbox[3]))
    
    # Load the variable (took ~400 seconds)
    st      = time.time()
    dsreg   = dsreg.load()
    print("Loaded in %.2fs" % (time.time()-st))
    
    dsreg   = standardize_names(dsreg)
    
    #%% Perform Preprocessing (Detrend, Deseason)
    
    # Remove seasonal cycle and detrend
    st      = time.time()
    dsa     = xrdeseason(dsreg) # Remove mean seasonal cycle
    dsa_dt  = xrdetrend(dsa) # Remove simple linear trend
    print("Deseason/Detrend in %.2fs" % (time.time()-st))
    
    #%% Do NAO Calculations
    
    N_mode    = 3
    months    = np.arange(1,13,1)
    
    eofmon    = []
    pcmon     = []
    varexpmon = []
    for mon in tqdm.tqdm(months):
        dsmon = selmon(dsa_dt,mon)
        
        # Use xEOFs to compute necessary information
        model           = xe.single.EOF(use_coslat=True,n_modes=N_mode)
        st              = time.time()
        model.fit(dsmon,dim='time')
        
        eofall          = model.components()
        pcall           = model.scores()
        varexpall       = model.explained_variance_ratio()
        print("Computed EOF in %.2fs" % (time.time()-st))
    
        # Need to set this from {} --> 'none' for to_netcdf to work later
        eofall.attrs['solver_kwargs']='none'
        pcall.attrs['solver_kwargs']='none'
        varexpall.attrs['solver_kwargs']='none'
        
        # Flip Signs where necessary
        spgbox     = [-60,20,45,80]
        eapbox     = [-60,20,45,60] # Shift Box west for EAP
        bbox_check = [spgbox,eapbox,]    
        print("Flipping boxes based on [bbox_check]")
        nmode_check = len(bbox_check)
        for N in range(nmode_check):
            chkbox = bbox_check[N]
            
            sumflx = eofall.isel(mode=N).sel(lon=slice(chkbox[0],chkbox[1]),lat=slice(chkbox[2],chkbox[3])).mean().data.item()
            
            if sumflx > 0:
                print("Flipping sign for SLP, mode %i" % (N+1))
                eofall[N,:,:] *= -1
                pcall[N,:] *= -1
        
        pcall['month'] = mon
        eofall['month'] = mon
        varexpall['month'] = mon
        
        # Append Output
        eofmon.append(eofall)
        pcmon.append(pcall)
        varexpmon.append(varexpall)
    
    # Concatenate and Save Output
    eofmon    = xr.concat(eofmon,dim='month')
    pcmon     = xr.concat(pcmon,dim='month')
    varexpmon = xr.concat(varexpmon,dim='month')
    dsout     = xr.merge([eofmon.rename('eof'),
                      pcmon.rename('pc'),
                      varexpmon.rename('varexp')])
    
    dsout.to_netcdf(ncname_out)
        
                
    # Repeat Again for All Months Together ========================================
    
    # Use xEOFs to compute necessary information
    model           = xe.single.EOF(use_coslat=True,n_modes=N_mode)
    st              = time.time()
    model.fit(dsa_dt,dim='time')
    eofall          = model.components()
    pcall           = model.scores()
    varexpall       = model.explained_variance_ratio()
    print("Computed EOF in %.2fs" % (time.time()-st))
    
    # Flip Signs where necessary
    spgbox     = [-60,20,45,80]
    eapbox     = [-60,20,45,60] # Shift Box west for EAP
    bbox_check = [spgbox,eapbox,]    
    print("Flipping boxes based on [bbox_check]")
    nmode_check = len(bbox_check)
    for N in range(nmode_check):
        chkbox = bbox_check[N]
        
        sumflx = eofall.isel(mode=N).sel(lon=slice(chkbox[0],chkbox[1]),lat=slice(chkbox[2],chkbox[3])).mean().data.item()
        
        if sumflx > 0:
            print("Flipping sign for SLP, mode %i" % (N+1))
            eofall[N,:,:] *= -1
            pcall[N,:] *= -1
    
    # Need to set this from {} --> 'none' for to_netcdf to work later
    eofall.attrs['solver_kwargs']='none'
    pcall.attrs['solver_kwargs']='none'
    varexpall.attrs['solver_kwargs']='none'
    
    dsout     = xr.merge([eofall.rename('eof'),
                      pcall.rename('pc'),
                      varexpall.rename('varexp')])
    dsout.to_netcdf(ncname_out_allmon)
    
    #%%
    # # Apply area weight
    # dsa_dt = dsa_dt.transpose('time','lat','lon')
    # wgt    = np.sqrt(np.cos(np.radians(dsa_dt.lat.values))) # [Lat]
    # dswgt  = dsa_dt.data * wgt[None,:,None]






        



