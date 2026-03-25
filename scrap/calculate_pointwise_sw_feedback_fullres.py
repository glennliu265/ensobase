#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Copied from regridded script, compute the pointwise sw feedback

Pipeline: <manual processing> --> [crop_SEP_sst.sh] --> [anom_detrend1_awiloop.sh]

Main Differences
    - Working at higher resolution for each simulation with pre-cropped data
    - Data has already been anomalized and linearly detrended via cdo [anom_detrend1_awiloop.sh]
    - 


Created on Tue Feb 24 10:26:30 2026

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

from sklearn.linear_model import LinearRegression
import sklearn

#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% Indicate Paths

# Get the components
expnames = ["TCo319_ctl1950d", "TCo319_ssp585", "TCo1279-DART-1950", "TCo1279-DART-2090", "TCo2559-DART-1950C"]
dpath    = "/home/niu4/gliu8/projects/scrap/SEP_crop/anom_detrend1/"
vnames   = ["sst","tscre"]

figpath  = "/home/niu4/gliu8/figures/bydate/2026-03-03/"
proc.makedir(figpath)

bbox_sep = [240,290,-50,10]
bbox     = bbox_sep


selmon    = [12,1,2]
if selmon is None:
    selmonstr = ""
    leadlags = np.arange(-12,13,1)
else:
    selmonstr = "_DJF"
    leadlags = np.arange(-6,7,1)



#%% Load Variables
nvar    = len(vnames)
nexp    = len(expnames)
#for ex in tqdm.tqdm(range(nexp)):
for ex in tqdm.tqdm(range(4,nexp+1)):

    
    expname = expnames[ex]
    
    ds_exp = []
    for vv in range(nvar):
        
        ncname = "%s%s_%s.nc" % (dpath,expnames[ex],vnames[vv])
        ds     = xr.open_dataset(ncname)[vnames[vv]].load()
        
        if 'time' in ds.dims:
            if len(ds['time']) < 2:  # Drop time if it is only 1 dimension long...
                print("Dropping existing time variable with length 1")
                ds = ds.squeeze().drop_vars('time')
        
        ds     = ut.standardize_names(ds)
        ds     = ut.remove_duplicate_times(ds)
        ds     = ut.varcheck(ds,vnames[vv],expnames[ex])
        
        ds     = proc.sel_region_xr(ds,bbox_sep)
        
        ds_exp.append(ds)
        
    #%% Preprocessing, then Calculate Lead Lag regression
    
    # Change NaN points to zero to make sure the ufunc doesnt break...
    
    def check_nan(ds):
        if np.any(np.isnan(ds)):
            print("NaN points detected, replacing with zeros")
            ds = xr.where(ds == 0.,np.nan,ds)
        return ds
    ds_exp=[check_nan(ds) for ds in ds_exp]
    
    # Check Time Matching
    sst,tscre   = ds_exp
    sst,tscre   = proc.match_time_month(sst,tscre)
    if np.any(sst['time'].data != tscre['time'].data):
        print("Warning, time is not matching exactly. Using time from first variable...")
        tscre['time'] = sst['time']
    
    # Select Months if applicable
    if selmon is not None:
        sst   = proc.selmon_ds(sst,[12,1,2])
        tscre = proc.selmon_ds(tscre,[12,1,2])
        
    
    st          = time.time()
    llreg       = ut.calc_leadlagreg_pointwise(tscre,sst,leadlags)
    proc.printtime(st,"Calculated")
    
    #%% Save the output
    llreg       = llreg.rename('beta') # I forgot to rename it so it might be '__xarray_dataarray_variable__'
    st          = time.time()
    outpath = "/home/niu4/gliu8/projects/scrap/SEP_crop/metrics"
    outname = "%s/Lag_Regression_%s_vs_%s_%s_LagMax%02i%s.nc" % (outpath,vnames[0],vnames[1],expname,leadlags[-1],selmonstr)
    llreg.to_netcdf(outname)
    proc.printtime(st,"Saved")
