#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

For CESM1 Simulations, Calculate Skewness using SciPy

Created on Tue Jul 28 16:42:07 2026

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

from tqdm.notebook import tqdm

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%%


def selyr(ds,ystart,yend):
    return ds.sel(time=slice(ystart+"-01-01",yend+"-12-31"))
    
def dailyclim(ds):
    return ds.groupby("time.dayofyear").mean("time")

def deseason_daily(ds,clim=False):
    dsclim       = dailyclim(ds)#ds.groupby("time.dayofyear").mean("time")
    dsanom       = ds.groupby("time.dayofyear") - dsclim
    if clim:
        return dsanom,dsclim
    return dsanom

def cubic_detrend(ds,return_all=False):
    # Get Coefficients (50.2 sec)
    cubicfit = ds.polyfit('time',deg=3)
    
    # Recover Model (13.5 sec)
    cubicfunc = xr.polyval(ds.time,cubicfit.polyfit_coefficients) 

    # Calculate Difference
    ds_detrend = ds-cubicfunc
    if return_all:
        return ds_detrend,cubicfit,cubicfunc
    return ds_detrend

def detrend_dim(da, dim="time", deg=1):
    # Function by Rohit Ghosh: https://github.com/rg568/EERIE_scripts/blob/main/FESOM/ENSO_Z500_DJF_teleconnection_IFS-FESOM.ipynb
    coeffs = da.polyfit(dim=dim, deg=deg)
    trend = xr.polyval(da[dim], coeffs.polyfit_coefficients)
    return da - trend
    
nanskew        = lambda x,axis: sp.stats.skew(x,axis=axis,nan_policy='omit')

#%%


# User Edits
scenario      ="BHIST"
expname       = "lores"
freq          = "day_1"
outpath       = "/home/niu4/gliu8/projects/mesaclip/skewness"
detrend_order = 2
vname         = "TS"
ystart        = '1980'
yend          = '2005'
seasonal      = True

# Set Number of Ensembles based on experiment
if expname == "lores":
    ensnum         = np.arange(1,23,1) # Update this Later!!
elif expname == 'hires':
    ensnum         = np.arange(1,11,1)

# Set path to Ensemble Member
rawpath_to_ens = "/home/niu4/gliu8/share/CESM1/MESACLIP/RDA/%s/%s/%s/regrid_1x1" % (expname,scenario,freq)
process_name   = "CESM1_MESACLIP_%s_%s_%s_Crop_%s-%s_detrend%i" % (expname,scenario,freq,ystart,yend,detrend_order)
if seasonal:
    process_name = process_name + "_seasonal"
proc.makedir("%s/%s/" % (outpath,process_name))

# Loop by Ensemble Mmeber
for ens in ensnum:
    
    # TS_18500101-18541231_regrid1x1.nc
    rawpath  = "%s/ens%03i" % (rawpath_to_ens,ens)
    ncsearch = "%s/%s_*_regrid1x1.nc" % (rawpath,vname)
    nclist   = glob.glob(ncsearch)
    nclist.sort()
    nfiles   = len(nclist)
    print("Found %i files for Ens %0i" % (nfiles,ens))
    
    # Open files and concatenate by time
    ds = xr.open_mfdataset(nclist,concat_dim='time',combine='nested')[vname]
    
    # Drop to selected time slice
    ds = selyr(ds,ystart,yend)
    print(len(ds.time))
    
    # Load the Dataset (All Time)
    st = time.time()
    ds = ds.load()
    print("Dataset Loaded in %.2fs" % (time.time()-st))
    
    # Preprocessing ===========================================================
    # LoRes Time (~53.97 sec 1980-2005)
    
    # Remove Mean Seasonal Cycle
    st  = time.time()
    dsa = deseason_daily(ds)
    
    # Remove Fitted Trend
    dsa = detrend_dim(dsa,deg=detrend_order)
    print("Anomalized in %.2fs" % (time.time()-st))
    
    # Compute Skewness using SciPy ============================================
    if seasonal:
        # LoRes time (~257.55s 1980-2005)
        st       = time.time()
        dsa_skew = dsa.groupby('time.season').reduce(func=nanskew,dim='time')
        print("Computed Skewness in %.2fs" % (time.time()-st))
        
    else:
        # LoRes time (~257.55s 1980-2005)
        st       = time.time()
        dsa_skew = dsa.reduce(func=nanskew,dim='time')
        print("Computed Skewness in %.2fs" % (time.time()-st))
    
    # Save Output =============================================================
    
    # Set Output Name
    st  = time.time()
    outname="%s/%s/%s_ens%03i.nc" %(outpath,process_name,vname,ens)
    # if seasonal:
    #     outname = proc.addstrtoext(outname,"_seasonal")
    dsa_skew.to_netcdf(outname)
    print("Saved in %.2fs" % (time.time()-st))
    print(outname)
    
    
    
    
    
    
    
    





