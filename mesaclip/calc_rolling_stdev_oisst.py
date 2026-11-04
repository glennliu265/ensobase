#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Rolling Standard Deviatino by Day of Year

Copied from [calc_rolling_threshold_oisst.py]
Created on Tue Nov 03 14:43 2026

@author: gliu

"""


# import sys
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
from datetime import datetime

from tqdm import tqdm

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

st_all = time.time()

#%% Functions

def rolling_stdev_threshold(timeseries,doy,ibefore,iafter,
                               sigma=1,return_count=False,noleap=False):
    
    if np.any(np.isnan(timeseries)):
        nanout = np.nan * np.ones((365))
        if return_count:
            return nanout,0
        return nanout
    
    
    # Get Date Ranges
    if noleap:
        ndays = 365
        days   = np.arange(1,366,1)
    else:
        ndays = 366
        days   = np.arange(1,367,1)
    dranges = []
    for dd in days:
        dmin = dd - ibefore
        dmax = dd + iafter
        if dmin <= 0:  # Crossing Jan 1
            dmin   = dmin + ndays
            drange = np.hstack([np.arange(dmin,ndays+1),np.arange(1,dmax+1)])
        elif dmax > ndays:  # Crossing Dec 31
            dmax = dmax - ndays
            drange = np.hstack([np.arange(dmin,ndays+1),np.arange(1,dmax+1)])
        else:
            drange = np.arange(dmin,dmax+1)
        dranges.append(drange)
    dranges = np.array(dranges) # [Day, Range]
    
    # Looping for each date, calculate the percentiles
    nsamples_byday = []
    stdev_byday = []
    for dd in range(len(days)):
        daysin   = dranges[dd,:]
        idchoose = np.where(np.isin(doy,daysin))[0]
        datasel  = timeseries[idchoose]
        qtl      = np.nanstd(datasel) * sigma
        stdev_byday.append(qtl)
        if return_count:
            nsamples_byday.append(len(datasel))
    stdev_byday = np.array(stdev_byday) # [Day, Quantile]
    if return_count:
        nsamples_byday = np.array(nsamples_byday)
        return stdev_byday,nsamples_byday
    return stdev_byday


#%% Open view of OISST (8GB for global)

st         = time.time()
procname   = "anom_detrend1_19820101-20251231_climatology1985to2014"  # anom_detrend2_19820101-20251231
ncpath     = "/home/niu4/gliu8/share/OISST/mergetest/%s/" % procname
ncname     = "oisst_day_1_sst_anom.nc"
dsview     = xr.open_dataset(ncpath+ncname).load()
xrname     = '__xarray_dataarray_variable__'
dsview     = dsview.convert_calendar('noleap') # Remove Leap Year
print("Loaded Data in %.2fs" % (time.time()-st))
print(dsview)

sigma      = 1 # number of standard deviations

outdir     = "/home/niu4/gliu8/projects/mesaclip/thresholds/%s/" % procname
proc.makedir(outdir)

# Threshold Selections
winsize    = 15
quantiles  = [0.10,0.90]
thresname  = "winsize%0i_std%i" % (winsize,sigma)

# NOTE NEED TO MANUALLY CHANGE OUTNAME
outname    = "%soisst_day_1_rolling_threshold_%s.nc" % (outdir,thresname) 

#%% Set Up the Function

doy=dsview.time.dt.dayofyear
st = time.time()

get_thres = lambda XX,YY: rolling_stdev_threshold(XX,YY,winsize,winsize,
                                                  sigma=sigma,
                                                  return_count=False,
                                                  noleap=True)

ds_thresholds   = xr.apply_ufunc(
    get_thres,
    dsview[xrname],
    doy,
    input_core_dims=[['time'],['time']],
    output_core_dims=[['doy']],
    vectorize=True,
    )
print("Computed threshold in %.2fs" % (time.time()-st))



ds_thresholds['doy']   = np.arange(1,366,1)
ds_thresholds['stdev'] = sigma
ds_thresholds          = ds_thresholds.squeeze()

#%% Save Output
ds_thresholds.to_netcdf(outname)
print("Calculated Threshold in %.2fs" % (time.time()-st_all))

