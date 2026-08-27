#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Daily Autocorrelation using XRCorr

Created on Thu Aug 27 09:46:58 2026

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

from tqdm import tqdm

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut


st_all = time.time()
#%% Open view of OISST (8GB for global)

st     = time.time()
ncpath = "/home/niu4/gliu8/share/OISST/mergetest/anom_detrend2_19820101-20251231/"
ncname = "oisst_day_1_sst_anom.nc"
dsview = xr.open_dataset(ncpath+ncname).load()
xrname = '__xarray_dataarray_variable__'
dsview = dsview.convert_calendar('noleap') # Remove Leap Year
print("Loaded Data in %.2fs" % (time.time()-st))
print(dsview)

outdir  = "/home/niu4/gliu8/projects/mesaclip/memory/oisst_byday/"

winsize = 5

#%%


#doyall  = ds.time.dt.dayofyear
nlags   = 365
lags    = np.arange(nlags+1)
days    = np.arange(1,366)
ds      = dsview.squeeze()[xrname]

years   = ds.time.dt.year
ystart  = years[0].data.item()
yend    = years[-1].data.item()

if winsize == 0:
    
    #lagcorr_alldoys = []
    for doybase in tqdm(np.arange(1,366)):
        
        base_var    = ds.sel(time=ds.time.dt.dayofyear.isin(doybase))
        
        lagcorr_day = []
        for ll in range(nlags):
            
            # Select Lag Variable
            doylag = doybase + ll
            if doylag > 365:
                doylag = doylag % 365
                year_cross = True
            else:
                year_cross = False
            lag_var  = ds.sel(time=ds.time.dt.dayofyear.isin(doylag))
            
            if year_cross:
                base_slice  = ["%04i-01-01" % (ystart), "%04i-12-31" % (yend-1)] 
                lag_slice   = ["%04i-01-01" % (ystart+1), "%04i-12-31" % (yend)] 
                
                base_var_in = base_var.sel(time=slice(*base_slice))
                lag_var_in  = lag_var.sel(time=slice(*lag_slice))
            
            else:
                base_var_in = base_var
                lag_var_in  = lag_var
            
            lag_var_in['time'] = base_var_in['time']
            
            corrout = xr.corr(base_var_in,lag_var_in,dim='time')
            lagcorr_day.append(corrout)
            
        lagcorr_day = xr.concat(lagcorr_day,dim='lag')
        
        outname = "%sdaily_acf_lagmax%i_nowindow_doy%03i.nc" % (outdir,nlags,doybase)
        lagcorr_day.to_netcdf(outname)
else:
    print("Running script with window size")
    
    for doybase in tqdm(np.arange(1,366)):
        
        
        # Construct Base Window
        baseinfo    = proc.construct_window_doy(doybase,winsize,verbose=False)
        
        lagcorr_day = []
        for ll in tqdm(range(nlags)):
            
            # Select Lag Variable
            doylag = doybase + ll
            if doylag > 365:
                doylag = doylag % 365
                year_cross = True
            else:
                year_cross = False
            
            # Construct Lag Widnow
            laginfo = proc.construct_window_doy(doylag,winsize,verbose=False)
            
            # Perform Indexing
            basevar,lagvar = proc.index_year_base_lag_xr(ds,baseinfo,laginfo,ll,verbose=False)
            
            # Calcualte xrcorr
            lagvar['time'] = basevar['time']
            corrout = xr.corr(basevar,lagvar,dim='time')
            lagcorr_day.append(corrout)
        
        lagcorr_day = xr.concat(lagcorr_day,dim='lag')
        outname = "%sdaily_acf_lagmax%i_winsize%02i_doy%03i.nc" % (outdir,nlags,winsize,doybase)
        lagcorr_day.to_netcdf(outname)
        
        # End DOY Loop

print("Completed Script in %.2fs" % (time.time()-st_all))
    

# #%% Quick Test if I just use xrcorr

# doy_monstart = proc.get_doy_monstart()

# st          = time.time()
# ll      = 1

# doybase = ds.sel(time=ds.time.dt.dayofyear.isin(1))

# print("Caluclated Lag in %.2fs" % (time.time()-st))

