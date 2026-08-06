#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Anomalize and Detrend OISST Output

Modeled after `merge_anom_detrend_mesaclip.py` and `merge_anomalized_regridded_MESACLIP.ipynb`

Created on Tue Aug  4 11:26:07 2026

@author: gliu
"""

import time
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import glob
import os

#%% Helper Functions

def remove_duplicate_times(ds,verbose=True,timename='time'):
    # From : https://stackoverflow.com/questions/51058379/drop-duplicate-times-in-xarray
    _, index = np.unique(ds[timename], return_index=True)
    print("Found %i duplicate times. Taking first entry." % (len(ds[timename]) - len(index)))
    return ds.isel({timename:index})

def detrend_dim(da, dim="time", deg=1):
    # Function by Rohit Ghosh: https://github.com/rg568/EERIE_scripts/blob/main/FESOM/ENSO_Z500_DJF_teleconnection_IFS-FESOM.ipynb
    coeffs = da.polyfit(dim=dim, deg=deg)
    trend = xr.polyval(da[dim], coeffs.polyfit_coefficients)
    return da - trend

def dailyclim(ds):
    return ds.groupby("time.dayofyear").mean("time")

def deseason_daily(ds,clim=False):
    dsclim       = dailyclim(ds)#ds.groupby("time.dayofyear").mean("time")
    dsanom       = ds.groupby("time.dayofyear") - dsclim
    dsanom       = dsanom.drop_vars('dayofyear')
    if clim:
        return dsanom,dsclim
    return dsanom

def makedir(expdir):
    """
    Check if "expdir" exists, and creates a directory if it doesn't

    Parameters
    ----------
    expdir : TYPE
        DESCRIPTION.

    """
    checkdir = os.path.isdir(expdir)
    if not checkdir:
        print(expdir + " Not Found! \n\tCreating Directory...")
        os.makedirs(expdir)
    else:
        print(expdir+" was found!")
    return None
    
#%% UserEdits

tstart  = "1982-01-01"
tend    = "2025-12-31" #"2025-12-31"
outpath = "/home/niu4/gliu8/share/OISST/mergetest/"
expname = "oisst"
vname   = "sst"
freq    = "day_1"
ncname  = "oisst_v2.1_daily_198109_20260630_regrid1x1.nc"
deg = 2


rawpath = "/home/niu4/gliu8/share/OISST/mergetest/" #% (scenario)

# Make Output Folder
tstart       = tstart.replace('-','')#npdatetime_to_str(dsanom.time[0]).replace('-','')
tend         = tend.replace('-','')  #npdatetime_to_str(dsanom.time[-1]).replace('-','')
outpath_proc = "%s/anom_detrend%i_%s-%s/" % (outpath,deg,tstart,tend,)
makedir(outpath_proc)

#%% Open File (15.85 Sec)

# Open View, Restrictto Time Slice, and Load
st    = time.time()
dsraw = xr.open_dataset(rawpath+ncname)[vname]
dsraw = dsraw.sel(time=slice(tstart,tend)).load()
print("Loaded file in %.2fs" % (time.time()-st))

#%% Remove Mean Seasonal Cycle (~9sec)
st    = time.time()
dsanom = deseason_daily(dsraw)
dsraw.close()
del dsraw
print("Deseasoned in %.2fs" % (time.time()-st))

#%% Detrend (26.31s)

st     = time.time()
dsanom = detrend_dim(dsanom,dim='time',deg=deg)
print("Detrended in %.2fs" % (time.time()-st))

#%% Save Output

st           = time.time()
outname      = "%s%s_%s_%s_anom.nc" % (outpath_proc,expname,freq,vname)
dsanom.to_netcdf(outname)
print("\tSaved to %s" % outname)
