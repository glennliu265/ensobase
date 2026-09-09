#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Rolling Percentile Thresholds for MHW/MCW Identification for OISST

Copied loading sections from [calc_doy_xrcorr.py]
Funtions + workflow copied from [Test_Threshold_Sensitivity.ipynb]


Created on Tue Sep  8 15:06:55 2026

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


def rolling_quantile_threshold(timeseries,doy,ibefore,iafter,
                               quantiles=[0.1,0.9],return_count=False,noleap=False):
    
    if np.any(np.isnan(timeseries)):
        nanout = np.nan * np.ones((365,2))
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
    quantile_byday = []
    for dd in range(len(days)):
        daysin   = dranges[dd,:]
        idchoose = np.where(np.isin(doy,daysin))[0]
        datasel  = timeseries[idchoose]
        qtl      = np.quantile(datasel,quantiles)
        quantile_byday.append(qtl)
        if return_count:
            nsamples_byday.append(len(datasel))
    quantile_byday = np.array(quantile_byday) # [Day, Quantile]
    if return_count:
        nsamples_byday = np.array(nsamples_byday)
        return quantile_byday,nsamples_byday
    return quantile_byday


#%% Open view of OISST (8GB for global)

st     = time.time()
ncpath = "/home/niu4/gliu8/share/OISST/mergetest/anom_detrend2_19820101-20251231/"
ncname = "oisst_day_1_sst_anom.nc"
dsview = xr.open_dataset(ncpath+ncname).load()
xrname = '__xarray_dataarray_variable__'
dsview = dsview.convert_calendar('noleap') # Remove Leap Year
print("Loaded Data in %.2fs" % (time.time()-st))
print(dsview)


outdir     = "/home/niu4/gliu8/projects/mesaclip/thresholds/anom_detrend2_19820101-20251231/"
proc.makedir(outdir)

# Threshold Selections
winsize    = 15
quantiles  = [0.10,0.90]
thresname  = "winsize%0i_pct%03i-%03i" % (winsize,quantiles[0]*100,quantiles[1]*100)

# NOTE NEED TO MANUALLY CHANGE OUTNAME
outname    = "%soisst_day_1_rolling_threshold_%s.nc" % (outdir,thresname) 





#%% Set Up the Function

doy=dsview.time.dt.dayofyear
st = time.time()

get_thres = lambda XX,YY: rolling_quantile_threshold(XX,YY,winsize,winsize,
                                                  quantiles=quantiles,
                                                  return_count=False,
                                                  noleap=True)

ds_thresholds   = xr.apply_ufunc(
    get_thres,
    dsview[xrname],
    doy,
    input_core_dims=[['time'],['time']],
    output_core_dims=[['doy','quantile']],
    vectorize=True,
    )
print("Computed threshold in %.2fs" % (time.time()-st))



ds_thresholds['doy'] = np.arange(1,366,1)
ds_thresholds['quantile'] = quantiles
ds_thresholds = ds_thresholds.squeeze()

#%% Save Output
ds_thresholds.to_netcdf(outname)
print("Calculated Threshold in %.2fs" % (time.time()-st_all))

#%%


hey =  rolling_quantile_threshold(xx[xrname],doy,winsize,winsize,
                                                  return_count=False,
                                                  noleap=True)
