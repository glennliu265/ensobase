#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Anomalize and Detrend OISST Output

Modeled after `merge_anom_detrend_mesaclip.py` and `merge_anomalized_regridded_MESACLIP.ipynb`

Updates
 - [2026.10.01] - Add Anom Detrend 1 (with selected climatology period)

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

def deseason_daily(ds,clim=False,period=None,verbose=True):
    if period is not None:
        if verbose:
            print("Taking climatology over the provided period: %s" % period)
        dsperiod = ds.sel(time=slice(*period))
        dsclim   = dailyclim(dsperiod)
    else:
        dsclim   = dailyclim(ds)#ds.groupby("time.dayofyear").mean("time")
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
rawpath = "/home/niu4/gliu8/share/OISST/mergetest/" #% (scenario)

expname = "oisst"
vname   = "sst"
freq    = "day_1" #"day_1"


if freq == "month_1":
    ncname  = "oisst_v2.1_monthly_198109_20260630_regrid1x1.nc"
else:
    ncname  = "oisst_v2.1_daily_198109_20260630_regrid1x1.nc"


deg     = 1
climperiod = ['1985-01-01','2014-12-31'] # Set Period to Calculate Climatology



# Make Output Folder
tstart       = tstart.replace('-','')#npdatetime_to_str(dsanom.time[0]).replace('-','')
tend         = tend.replace('-','')  #npdatetime_to_str(dsanom.time[-1]).replace('-','')
procstr = "anom_detrend%i_%s-%s" % (deg,tstart,tend,)
if climperiod is not None:
    climname     = "_climatology%sto%s" % (climperiod[0][:4],climperiod[1][:4])
    procstr      = procstr+climname
if freq == "month_1":
    outpath = "/home/niu4/gliu8/share/OISST/mergetest/monthly/"
else:
    outpath = "/home/niu4/gliu8/share/OISST/mergetest/"
outpath_proc = "%s/%s/" % (outpath,procstr)
makedir(outpath_proc)


#%% Open File (15.85 Sec)

# Open View, Restrictto Time Slice, and Load
st    = time.time()
dsraw = xr.open_dataset(rawpath+ncname)[vname]
dsraw = dsraw.sel(time=slice(tstart,tend)).load()
print("Loaded file in %.2fs" % (time.time()-st))

#%% Remove Mean Seasonal Cycle (~9sec)
st    = time.time()
dsanom = deseason_daily(dsraw,period=climperiod)
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
