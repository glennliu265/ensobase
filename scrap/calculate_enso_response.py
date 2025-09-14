#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

1) Load Nino3.4 Index (calc_nino34)
2) Load Target Variable
3) Preprocess Target Variable (Detrend and Deseason)

4) Regress Index onto anomalized variable to obtain instantaneous ENSO response (all months)
5) Regress Index separately month by month to obtain monthwise-response

Created on Fri Sep 12 15:06:47 2025

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


#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

#%% Indicate paths

figpath         = "/home/niu4/gliu8/figures/bydate/2025-09-16/"
proc.makedir(figpath)

datpath         = "/home/niu4/gliu8/projects/scrap/TP_crop/"

expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090"]

vname           = "D20" #"sst"#"str"
"""
Some issues to fix with vname

tx_surf_1m  (file name) corresponds to tx_sur
D20         (file name) corresponds to nz1
Dmax        (file name) corresponds to nz1

in one file
SSR (variable name) --> ssr, and the coordinate names are somewhat messed up

Also, it seems that I spelled separately wrong in the filename, RIP

"""


ninopath        = "/home/niu4/gliu8/projects/scrap/nino34/"

nexps           = len(expnames)

#%% Helper Functions


def swap_rename(ds,chkvar,newvar):
    if chkvar in list(ds.coords):
        print("Renaming [%s] to [%s]" % (chkvar,newvar))
        ds = ds.rename({chkvar:newvar})
    return ds

def standardize_names(ds):
    
    ds = swap_rename(ds,'time_counter','time')
    ds = swap_rename(ds,"TIME_COUNTER",'time')
    ds = swap_rename(ds,"LON","lon")
    ds = swap_rename(ds,"LAT","lat")
    ds = swap_rename (ds,"LAT232_409","lat")
    return ds
    
# def rename_time_dart(ds):
#     if 'time_counter' in list(ds.coords):
#         print("Renaming [time_counter] to [time]")
#         ds = ds.rename({'time_counter':'time'})
#     if 'TIME_COUNTER' in list(ds.coords):
#         print("Renaming [TIME_COUNTER] to [time]")
#         ds = ds.rename({'TIME_COUNTER':'time'})
#     return ds

def preprocess_enso(ds):
    # Remove Mean Seasonal Cycle and the Quadratic Trend
    dsds   = proc.xrdeseason(ds)
    
    # Compute Spectra
    def detrend_quadratic(ds):
        x = np.arange(len(ds))
        y = ds.data
        ydetrended,model=proc.detrend_poly(x,y,2)
        return ydetrended
        
    st = time.time()
    dsanom = xr.apply_ufunc(
        detrend_quadratic,  # Pass the function
        dsds,  # The inputs in order that is expected
        # Which dimensions to operate over for each argument...
        input_core_dims=[['time'],],
        output_core_dims=[['time'],],  # Output Dimension
        vectorize=True,  # True to loop over non-core dims
        )
    print("Detrended in %.2fs" % (time.time()-st))
    #dsanom = proc.xrdetrend_1d(ds,2) 
    return dsanom


def match_time_month(var_in,ts_in):
    
    # Note works for datetime64[ns] format in xr.DataArray
    
    if len(var_in.time) != len(ts_in.time): # Check if they match
        
        # Warning: Only checking Year and Date
        vstart = str(np.array((var_in.time.data[0])))[:7]
        tstart = str(np.array((ts_in.time.data[0])))[:7]
        
        if vstart != tstart:
            print("Start time (v1=%s,v2=%s) does not match..." % (vstart,tstart))
            if vstart > tstart:
                print("Cropping to start from %s" % vstart)
                ts_in = ts_in.sel(time=slice(vstart+"-01",None))
            elif vstart < tstart:
                print("Cropping to start from %s" % tstart)
                var_in = var_in.sel(time=slice(tstart+"-01",None))
        
        vend = str(np.array((var_in.time.data[-1])))[:7]
        tend = str(np.array((ts_in.time.data[-1])))[:7]
        
        
        if vend != tend:
            
            print("End times (v1=%s,v2=%s) does not match..." % (vend,tend))
            
            if vend > tend:
                print("\nCropping to start from %s" % tend)
                var_in = var_in.sel(time=slice())
            elif vend < tend:
                print("\nCropping to start from %s" % vend)
                ts_in = ts_in.sel(time=slice(None,vend+"-31"))
                
        print(len(var_in.time) == len(ts_in.time))  
    return var_in,ts_in
        

#%% Part (1) Load the ENSO Indices

ds_enso = []
for ex in range(nexps):
    
    ninonc      = "%s%s_nino34.nc" % (ninopath,expnames[ex])
    ds          = xr.open_dataset(ninonc).load()
    
    
    ds_enso.append(ds)

#%% Part (2) Load the variable of choice

ds_var = []
for ex in tqdm.tqdm(range(nexps)):
    
    
    if vname == "tx_sur":
        vname_file = "tx_surf_1m"
        ncname = "%s%s_%s.nc" % (datpath,expnames[ex],vname_file)
    else:
        ncname = "%s%s_%s.nc" % (datpath,expnames[ex],vname)
    ds = xr.open_dataset(ncname)
    
    if vname.upper() in list(ds.keys()):
        print("Renaming %s to %s" % (vname.upper(),vname))
        ds = ds.rename({vname.upper():vname})
    
    ds_var.append(ds)

ds_var = [standardize_names(ds) for ds in ds_var]

#%% Part (3) Deseason and Detrend

# Takes
# 1.32 sec
# 13.06 sec (exp2)
# 110.36s
# 94.95s
if vname in ["D20","Dmax"]:
    ds_anoms = [preprocess_enso(ds['nz1']) for ds in ds_var]
else:
    ds_anoms = [preprocess_enso(ds[vname]) for ds in ds_var]



#%% Part (4) Compute the regression pattern, for all months and separately. then save file.

# All Months
for ex in range(4):
    ts_in   = ds_enso[ex].sst
    var_in  = ds_anoms[ex]
    
    var_in,ts_in = match_time_month(var_in,ts_in)
    
    dsout   = proc.detrend_by_regression(var_in, ts_in, regress_monthly=False)
    
    edict   = proc.make_encoding_dict(dsout)
    ncout   = "%s%s_%s_ENSO_regression_allmonths.nc" % (datpath,expnames[ex],vname)
    
    dsout.to_netcdf(ncout,encoding=edict)
    
# Monthly
for ex in range(4):
    ts_in   = ds_enso[ex].sst
    var_in  = ds_anoms[ex]
    
    var_in,ts_in = match_time_month(var_in,ts_in)
    
    dsout   = proc.detrend_by_regression(var_in, ts_in, regress_monthly=True)
    
    edict   = proc.make_encoding_dict(dsout)
    ncout   = "%s%s_%s_ENSO_regression_seperatemonths.nc" % (datpath,expnames[ex],vname)
    
    dsout.to_netcdf(ncout,encoding=edict)

