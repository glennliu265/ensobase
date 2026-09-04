#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate daily lagged autocorrelation function by ensemble member

Copied output for calc_doy_xrcorr

Created on Fri Sep  4 09:31:44 2026

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


#%% User Edits

# Experiment Information
ncpath        = "/home/niu4/gliu8/share/CESM1/MESACLIP/processed/anom_detrend2_19820101-20251231/"
xrname        = '__xarray_dataarray_variable__'
expname       = "mesaclip_lores" #"mesaclip_hires"
vname         = "SST"
ens_restrict  = np.arange(1,11,1)#None

fill_val      = None

# Output Information
outpath       = "/home/niu4/gliu8/projects/mesaclip/memory/anom_detrend2_19820101-20251231/"

# DOY Calculation Choices
winsize       = 0
nlags         = 365
lags          = np.arange(nlags+1)
doy_selection      = None
doy_selection_name = None
if doy_selection is None:
    days         = np.arange(1,366)
    doystr       = "all"
else:
    days         = doy_selection
    doystr       = doy_selection_name
ndoy = len(days)

        
#%%


# Get Ensemble Numbers to Loop From
if ens_restrict is not None:
    ensnum     = ens_restrict
    print("Restricting to Members: %s" % ens_restrict)
else:
    print("Selecting all Members")
    if "hires" in expname: # Loop for 10
        ensnum = np.arange(1,11,1)
    else: # Loop for 40
        ensnum = np.arange(1,41,1)
        
nens = len(ensnum)
print("Looping for %i members!" % nens)

#%%

for e in range(nens):
    st     = time.time()
    ens    = ensnum[e]
    
    # Load NetCDF
    # ex: mesaclip_lores_day_1_SST_anom_ens001.nc
    ncname = "%s%s_day_1_%s_anom_ens%03i.nc" % (ncpath,expname,vname,ens)
    dsview = xr.open_dataset(ncname).load()
    ds     = dsview.squeeze()[xrname]
    print("Loaded Data in %.2fs" % (time.time()-st))
    
    
    # Set Output Name (add string later)
    outdir       = "%s/%s/ens%03i" % (outpath,expname,ens)
    proc.makedir(outdir)
    outstr_start = "%s/%s_daily_acf_doy_%s_winsize%02i_lagmax%03i" % (outdir,expname,
                                                                              doystr,winsize,nlags,)
    
    
    # Part Below here should be identical to calc_doy_xrcorr
    # Get Necessary Variables
    years   = ds.time.dt.year
    ystart  = years[0].data.item()
    yend    = years[-1].data.item()
    
    if winsize == 0:
        print("Not using any window")
        #lagcorr_alldoys = []
        for dd in tqdm(range(ndoy)):
            doybase    = days[dd]
            
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
                
                # Check if # of years is the same ------------------------
                # If not, either skip (for missing year in doy_base) or fill with fill_val (for missing_year in doy_lag)
                # Note that for mesaclip_lr, 2006-01-02 is missing?
                byr = base_var.time.dt.year
                lyr = lag_var.time.dt.year
                nb = len(byr)
                nl = len(lyr)
                if nb != nl:
                    print("Warning! Uneven year lengths (base: %i, lag: %i)" % (nb,nl))
                    if nb < nl:
                        iyrs             = np.where(~np.isin(lyr,byr))[0]
                        print("Filling %i Missing Times in [base_var] with mean (or [fill_val)..." % (len(iyrs)))
                        fill_val         = base_var.mean('time') 
                        if fill_val is not None:
                            fill_val         = xr.ones_like(fill_val) * fill_val
                        
                        missing_concat   = []
                        for iyr in iyrs:
                            
                            missingyr               = lyr[iyr]
                            byr_day_average         = fill_val.copy()
                            # Get Missing Date (taking mon and day from first element of byr)
                            time_missing            = "%04i-%02i-%02i" % (missingyr,byr.time[0].data.item().month,byr.time[0].data.item().day)
                            # Use xr.date_range to create cftime object
                            tmiss                   = xr.date_range(start=time_missing,periods=1,calendar='noleap',use_cftime=True)[0]
                            # Assign to mean Array
                            byr_day_average['time'] = tmiss
                            missing_concat.append(byr_day_average.copy())
                        base_var_new = xr.concat([base_var,] + missing_concat,dim='time').sortby('time')
                        base_var     = base_var_new.copy()
                        del base_var_new
                    else:
                        iyrs             = np.where(~np.isin(byr,lyr))[0]
                        print("Filling %i Missing Times in [lag_var] with mean (or [fill_val)..." % (len(iyrs)))
                        fill_val         = lag_var.mean('time') 
                        if fill_val is not None:
                            fill_val         = xr.ones_like(fill_val) * fill_val
                        
                        missing_concat   = []
                        for iyr in iyrs:
                            
                            missingyr               = byr[iyr]
                            lyr_day_average         = fill_val.copy()
                            # Get Missing Date (taking mon and day from first element of lyr)
                            time_missing            = "%04i-%02i-%02i" % (missingyr,lyr.time[0].data.item().month,lyr.time[0].data.item().day)
                            # Use xr.date_range to create cftime object
                            tmiss                   = xr.date_range(start=time_missing,periods=1,calendar='noleap',use_cftime=True)[0]
                            # Assign to mean Array
                            lyr_day_average['time'] = tmiss
                            missing_concat.append(lyr_day_average.copy())
                        lag_var_new = xr.concat([lag_var,] + missing_concat,dim='time').sortby('time')
                        lag_var     = lag_var_new.copy()
                        del lag_var_new
                    # End Condition for mismatched base/lag year lengths
                # ------------------------
                
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
            
            # Set Output Name
            outname = "%s_doy%03i.nc" % (outstr_start,doybase)
            lagcorr_day.to_netcdf(outname)
            
    else:
        
        print("Running script with window size")
        
        for dd in tqdm(range(ndoy)):
            
            doybase = days[dd]
            
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
            
            # Set Output Name
            outname = "%s_doy%03i.nc" % (outstr_start,doybase)
            lagcorr_day.to_netcdf(outname)
            
            # End DOY Loop
            
        # End Window Conditional
    
    print("Completed Ens %03i in %.2fs" % (ens,time.time()-st))
    # End Ens Loop

print("Completed Script in %.2fs" % (time.time()-st_all))

