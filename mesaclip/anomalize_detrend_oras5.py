#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Copied from Anomalize_Detrend_OISST.py

Repeat similar procedure for ORAS5 Vertical Depth Output
Processes output pre-processed by `regrid_level_oras5.sh` and `consolidate_oras5.sh`

Created on Tue Aug 11 11:17:03 2026

@author: gliu

"""
import time
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import glob
import os

import sys

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut


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

tstart       = "1982-01-01"
tend         = "2025-12-31" #"2025-12-31"
rawpath      = "/home/niu4/gliu8/share/ORAS5/regridded/" #% (scenario)
climperiod = ['1985-01-01','2014-12-31'] # Set Period to Calculate Climatology
outpath      = rawpath
expname      = "oras5"
vname        = "votemper"
vname_new    = "TEMP"
freq         = "month_1" #"day_1"
deg          = 1
levels       = [50,100,200,500,1000]

# Make Output Folder
tstart       = tstart.replace('-','')#npdatetime_to_str(dsanom.time[0]).replace('-','')
tend         = tend.replace('-','')  #npdatetime_to_str(dsanom.time[-1]).replace('-','')
procstr = "anom_detrend%i_%s-%s" % (deg,tstart,tend,)
if climperiod is not None:
    climname     = "_climatology%sto%s" % (climperiod[0][:4],climperiod[1][:4])
    procstr      = procstr+climname
if freq == "month_1":
    outpath = rawpath #rawpath + "monthly/" /home/niu4/gliu8/share/OISST/mergetest/monthly/"
else:
    print("Only monthly resolution (month_1=TEMP) is supported by ORAS5, TEMP 3-D data")
    outpath = rawpath #"/home/niu4/gliu8/share/OISST/mergetest/"
outpath_proc = "%s/%s/" % (outpath,procstr)
makedir(outpath_proc)
ystart       = tstart[:4]
yend         = tend[:4]

makedir(outpath_proc)

#%% Now Loop for each level
nlevels = len(levels)
for ll in range(nlevels):
    
    st        = time.time()
    level     = levels[ll]
    vname_out = "%s%i" % (vname_new,level)
    print("Now Processing %s" % vname_out)
    
    # Open and Load the File
    st        = time.time()
    ncname    = "%s%s_level%s_%s-%s_regridded.nc" % (rawpath,vname,level,ystart,yend)
    dsraw     = xr.open_dataset(ncname)
    renamedict = {vname : vname_out, 'time_counter':'time'} # Rename Dimensions and subset to time
    dsraw      = dsraw.rename(renamedict)
    dsraw     = dsraw.sel(time=slice(tstart,tend)).load()[vname_out]
    print("File loaded in %.2fs" % (time.time()-st))
    
    
    
    # Remove Mean Seasonal Cycle
    st     = time.time()
    dsanom = proc.xrdeseason(dsraw)
    del dsraw
    print("Deseasoned in %.2fs" % (time.time()-st))
    
    # Detrend
    st     = time.time()
    dsanom = detrend_dim(dsanom,dim='time',deg=deg)
    print("Detrended in %.2fs" % (time.time()-st))
    
    # Save Output
    st           = time.time()
    outname      = "%s%s_%s_%s_anom.nc" % (outpath_proc,expname,freq,vname_out)
    print(outname)
    dsanom.to_netcdf(outname)
    print("\tSaved to %s" % outname)




