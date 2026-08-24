#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 19:13:04 2026

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

from tqdm.notebook import tqdm

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut


#%% Open view of OISST (8GB for global)

st     = time.time()
ncpath = "/home/niu4/gliu8/share/OISST/mergetest/anom_detrend2_19820101-20251231/"
ncname = "oisst_day_1_sst_anom.nc"
dsview = xr.open_dataset(ncpath+ncname).load()
print("Loaded Data in %.2fs" % (time.time()-st))
print(dsview)

#%% Set Function

# User Edits


def get_doy(timeseries,doy,days):
    idchoose = np.where(np.isin(doy,days))[0]
    datasel  = timeseries[idchoose]
    return datasel

def lagcorr_day(timeseries,doy,doy_start,nlag,verbose=False):
    if np.any(np.isnan(timeseries)):
        return np.zeros(len(np.arange(0,nlags))) * np.nan
    
    # Func Start Here
    doy_end   = doy_start + nlags
    
    # Part 1, Get Days of year
    if doy_end > 365:
        year_cross = True
        doy_end    = doy_end %  365
        if verbose:
            print("End Date crosses Dec 31")
            print("End Day is now %i" % doy_end)
        
    
        days_sel_y0 = np.arange(doy_start,366,1)  # From day to Dec 31
        days_sel_y1 = np.arange(1,doy_end+1,1)    # From Jan 1 to end (next year)
        
        tsin_y0 = get_doy(timeseries,doy,days_sel_y0) # Get Data Indexing by Day of Year
        tsin_y1 = get_doy(timeseries,doy,days_sel_y1) # 
        if verbose:
            print("Selecting Days %i to %i (Y0)" % (doy_start,365))
            print("Selecting Days %i to %i (Y1)" % (1,doy_end))
            print("Total Selected Days: %i" % (len(days_sel_y0) + len(days_sel_y1)))
        
    else:
        year_cross =False
        days_sel = np.arange(doy_start,doy_end+1,1)   # Indicate which days are being selected
        
        tsin = get_doy(timeseries,doy,days_sel)       # Get Data Indexing by Day of Year
        if verbose:
            print("Selecting Days %i to %i" % (doy_start,doy_end)) 
            print("Total Selected Days: %i" % (len(days_sel))) 
            print("\n")
    
    # Part 2, Concatenate and perform calculations
    if year_cross:
        # First Timesries
        ndays0 = len(days_sel_y0)        # Number of days selected
        nyrs0  = int(len(tsin_y0)/ndays0) # Number of years
        tsin0  = tsin_y0.reshape(nyrs0,ndays0) # Reshape to yr x day
    
        # Second Timeseries
        ndays1 = len(days_sel_y1)
        nyrs1  = int(len(tsin_y1)/ndays1)
        tsin1  = tsin_y1.reshape(nyrs1,ndays1)
    
        if nyrs1 != nyrs0:
            print("Warning! Number of years is not equal! (%i, %i)" % (nyrs0,nyrs1))
    
        # Stack along Day of Year, lagging the year
        tsin_rs = np.concatenate([tsin0[:(nyrs0-1),:],
                               tsin1[1:,:],
                              ],axis=1)
        ndays = ndays1 + ndays0
    else:
    
        ndays = len(days_sel)
        nyrs  = int(len(tsin)/ndays)
        tsin_rs  = tsin.reshape(nyrs,ndays)
    
    corrout  = []
    lags = np.arange(0,nlags)
    for ll in range(nlags):
        lag    = lags[ll]
        tsbase = tsin_rs[:,:(ndays-lag)]
        tslag  = tsin_rs[:,lag:]
        corrout.append(np.corrcoef(tsbase.flatten(),tslag.flatten())[0,1])
    corrout = np.array(corrout)
    return corrout

def calc_all_doy_lags(timeseries,doy,nlags=300):
    corr_bydoy = []
    for dd in np.arange(1,366):
        corrout = lagcorr_day(timeseries,doy,dd,nlags,verbose=False)
        corr_bydoy.append(corrout)
    corr_bydoy = np.array(corr_bydoy) # [doy,lag]
    return corr_bydoy


# Use Above as Numpy -- It is much faster (~2 min instead of 20 min...)
def seldoy_np(timeseries,doys,doysel):
    idsel = np.where(np.isin(doys,doysel))[0]
    return timeseries[idsel]


# NOTE THIS IS THE WORKING VERSION RN for NO WINDOW
def lagcorrdaily_nowindow_np(timeseries,years,doys,nlags,verbose=False):
    if np.any(np.isnan(timeseries)):
        lags          = np.arange(1,nlags+1,1)
        nlagloop      = len(lags)
        lagcorr_byday = np.zeros((365,nlagloop)) * np.nan
        
        return lagcorr_byday
    
    # Get Year Start and Year End
    year_start,year_end = years[0],years[-1]
    
    # Set Up Lags and Start Days
    lags          = np.arange(1,nlags+1,1)
    nlagloop      = len(lags)
    doy_all       = np.arange(1,366,1)
    
    # Loop
    lagcorr_byday = np.zeros((365,nlagloop))
    for dd in tqdm(range(365)):
        doy_start     = doy_all[dd]

        
        year_shift    = 0
        shift_counter = 0
        for ll in range(nlagloop):
            lag       = lags[ll]
            
            # Determine Start and Lag Day
            lag       = lags[ll]
            doy_base  = doy_start
            doy_lag   = doy_start + lag
            if doy_lag > 365: # Move to Next Year
                doy_lag = doy_lag % 365
                # Add to Year Shift for First Crossing, or Nearing 1-year since last
                if (shift_counter == 0) or (shift_counter == 364): 
                    year_shift += 1   
                    shift_counter = 0 # Reset to 0
                    if verbose:
                        print("Year shift detected at Lag %i" % lag)
            # if verbose:
            #     print("Lag Day: %i" % doy_lag)
            
            # Select Days
            dsbase    = seldoy_np(timeseries,doys,doy_base)
            dslag     = seldoy_np(timeseries,doys,doy_lag)
            baseyears = seldoy_np(years,doys,doy_base)
            lagyears  = seldoy_np(years,doys,doy_lag)
            
            # Select Years (if Year Crossing)
            if year_shift > 0: # Apply Shift in Years for Lag
                base_years = [year_start,year_end-year_shift]# np.arange(year_start,year_end-year_shift+1)
                lag_years  = [year_start+year_shift,year_end]#np.arange(year_start+year_shift,year_end+1)
    
                dsbase     = dsbase[ (baseyears>=base_years[0]) & (baseyears<=base_years[1]) ]
                dslag      = dslag[ (lagyears>=lag_years[0])    & (lagyears<=lag_years[1]) ]

                shift_counter += 1
            try:
                corrout = np.corrcoef(dsbase,dslag)[0,1]
            except:
                print("Failure, sizes are:")
                print(dsbase.shape)
                print(dslag.shape)
                corrout = np.nan
            lagcorr_byday[dd,ll] = corrout.copy() 
    return lagcorr_byday


#%%

st            = time.time()
doy_start     = 1   # Day of Year to Start
nlags         = 300 # # of days to include
timeseries    = dsview#.data # Input Timeseries
doy           = dsview.time.dt.dayofyear
years         = dsview.time.dt.year
xrname        = '__xarray_dataarray_variable__'
timeseries    = timeseries[xrname].squeeze()
use_xrfunc    = True
winsize       = 0


# This is the Loop Version
if not use_xrfunc:
    doy               = doy.squeeze()
    timeseries        = timeseries.transpose('time','lat','lon',)
    ntime,nlat,nlon   = timeseries.shape
    years             = years.data
    doys              = doy.data
    
    nday              = 365
    nlags             = nlags+1
    decorr_timescales = np.zeros((nday,nlags,nlat,nlon)) * np.nan # 57 GB for 300 Lags...
    
    if winsize == 0:
        print("Window Size is Zero")
        for a in tqdm(range(nlat)):
            for o in range(nlon):
                tspt     = timeseries.isel(lat=a,lon=o).data
                if np.any(np.isnan(tspt)):
                    continue
                decorrpt = lagcorrdaily_nowindow_np(tspt,years,doys,nlags,verbose=False)
                decorr_timescales[:,:,a,o] = decorrpt.copy()
        
    else: # Old Version
        print("Calc using old func")
        for a in tqdm(range(nlat)):
            for o in range(nlon):
                tspt     = timeseries.isel(lat=a,lon=o).data.squeeze()
                if np.any(np.isnan(tspt)):
                    continue
                decorrpt = calc_all_doy_lags(tspt,doy,nlags=nlags)
                decorr_timescales[:,:,a,o] = decorrpt.copy()
            
    dictout = dict(
        day_of_year = np.arange(1,366,1),
        lag = np.arange(1,nlags+1,1),
        lat = timeseries.lat.data,
        lon = timeseries.lon.data,
        )
    
    ds_doy = xr.DataArray(decorr_timescales,dims=dictout,coords=dictout,name='decorrelation_timescales')
    
    ds_doy.to_netcdf("doy_test_calc_pointwise_winsize%i.nc" % winsize)
    print("Completed Calculation in %.2fs" % (time.time()-st))  
    
else:
    print("Use xrfunc")
    
    if winsize == 0:
        
        print("\tUsing Winsize = 0")
        
        ds_doy= xr.apply_ufunc(
            lagcorrdaily_nowindow_np,
            timeseries,
            doy,
            input_core_dims=[['time'],['time']],
            output_core_dims=[['doy','lags',]],
            vectorize=True,
            )
        
        
    else:
        print("\tUsing Winsize = 300")
        
        # This is the xrfunc Version (took 5 days to run for OISST :(...)
        ds_doy= xr.apply_ufunc(
            calc_all_doy_lags,
            timeseries,
            years,
            doy,
            nlags,
            input_core_dims=[['time'],['time'],['time'],[]],
            output_core_dims=[['doy','lags',]],
            vectorize=True,
            )
    
    ds_doy.to_netcdf("doy_test_calc_pointwise_winsize%i.nc" % winsize)
    
print("Completed Calculation in %.2fs" % (time.time()-st))
    
    

