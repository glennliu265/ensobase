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
print("Loaded Data in %.2fs" % (time.time()))
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
        

    

st = time.time()
doy_start     = 1   # Day of Year to Start
nlags         = 300 # # of days to include
timeseries    = dsview#.data # Input Timeseries
doy           = dsview.time.dt.dayofyear

ds_doy= xr.apply_ufunc(
    calc_all_doy_lags,
    timeseries,
    doy,
    input_core_dims=[['time'],['time']],
    output_core_dims=[['doy','lags',]],
    vectorize=True,
    )


ds_doy.to_netcdf("doy_test_calc.nc")

print("Completed Calculation in %.2fs" % (time.time()-st))



