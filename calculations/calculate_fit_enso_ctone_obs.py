#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Fit ENSO Combination Mode to Flux Anomalies

Copied wholeperiod script

Created on Mon Sep 14 11:27:15 2026

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
import scipy as sp

from tqdm import tqdm

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

tbxpath = "/home/niu4/gliu8/scripts/commons/tbx"
sys.path.append(tbxpath)
import tbx as tbx

#%% Functions

# Do some preprocessing
def preproc(ds):
    dsanom = proc.xrdeseason(ds)
    dsanom = proc.xrdetrend_dim(dsanom,dim='time',deg=2)
    return dsanom

def preproc_dataset(ds,varlist):
    dsanoms = []
    for vname in varlist:
        dsa = preproc(ds[vname])
        dsanoms.append(dsa.rename(vname))
    return xr.merge(dsanoms)

def apply_to_dataset(ds,varlist,func):
    # Apply function that takes data array as argument to all specific variables of dataset
    dsnew = []
    for vname in varlist:
        dsa = func(ds[vname])
        dsnew.append(dsa.rename(vname))
    return xr.merge(dsnew)

def fit_enso_ctone_phi(dsin,ninoin,fillval=0,save_model=True):
    # Given anomaly timeseries [dsin] and nino3.4 index [ninoin], all np.arrays,
    # Compute the regression coefficient for beta*nino34(t)*cos(wt+phi) for a 
    # set of 12 phis (1 for each month), and 1/w = seasonal cycle
    
    y = dsin                 # Monthly Anomaly Timeseries to Fit
    x = np.arange(len(dsin)) # Indices representing time dimension

    # Replace NaNs with fillval (default is zero)
    ninoin = np.where(np.isnan(ninoin),fillval,ninoin)
    y      = np.where(np.isnan(y),fillval,y)
    
    # Set Phis
    # Currently fits 12 phis, at intervals of pi/6
    phis   = np.arange(12) * np.pi/6
    
    # Set Freq (Omega), set to seasonal cycle
    omega   = (2*np.pi)/12
    
    betas  = [] # Regression Slope [12]
    r2s    = [] # Fit to Original Timeseries [12]
    ypreds = [] # Modeled Timeseries [12 x ntime]
    for phi in phis:
        
        # Make the Function
        def funcfit(t,beta):
            return beta * ninoin * np.cos( omega*t - phi)
        
        # Use Scipy Optimize to obtain beta
        params, covariance = sp.optimize.curve_fit(funcfit, x, y)

        # Calculate the fit
        ymodel = funcfit(x,params[0])
        r2     = np.corrcoef(y,ymodel)[0,1]**2
        
        # Save Variables
        betas.append(params[0].item())
        r2s.append(r2.item())
        if save_model:
            ypreds.append(ymodel)
    if save_model:
        return np.array(betas),np.array(r2s),np.array(ypreds)
    return np.array(betas),np.array(r2s) # Otherwise, Don't Output Model

def pointwise_fit_enso_ctone_phi(anom_in,ninoin,save_model=False):
    st = time.time()
    if save_model:
        funcin  = lambda a,b: fit_enso_ctone_phi(a,b,save_model=True)
        outcore = [['month'],['month'],['month','time']]
    else:
        funcin  = lambda a,b: fit_enso_ctone_phi(a,b,save_model=False)
        outcore = [['month'],['month'],]
    
    # Apply to Each Point
    cfitout = xr.apply_ufunc(
        funcin,
        anom_in,
        ninoin,
        input_core_dims=[['time'],['time']],
        output_core_dims=outcore,
        vectorize=True,
        )
    
    # Make into DataSet
    betas   = cfitout[0].rename('beta')
    r2s     = cfitout[1].rename('r2')
    dsout   = [betas,r2s]
    if save_model:
        ymodels = cfitout[2].rename('ymodel')
        dsout   = dsout + [ymodels,]
    dsout = xr.merge(dsout)
    dsout['month'] = np.arange(1,13,1) # Assign Proper Months
    print("Completed Fit in %.2fs" % (time.time()-st))
    return dsout

def concat_exp(ds_by_variable):
    # Copied from compare_sst_cre_ssp585.ipynb
    testcat = xr.concat(ds_by_variable,dim='time')
    ntime_cat = len(testcat.time)
    dummy_time = xr.date_range('0000-01-01',periods=ntime_cat,freq="MS",calendar='noleap',use_cftime=True)
    testcat['time'] = dummy_time
    return testcat


def pointwise_lp(ds_raw,cutoffmon,order=6):
    # Copied from simple_mode_ctone.ipynb
    # xrfunc Application of movmean
    st        = time.time()
    apply_lp  = lambda ds: proc.lp_butter(ds.data,cutoffmon,order)
    ds_smooth= xr.apply_ufunc(
        apply_lp,
        ds_raw,
        input_core_dims=[['time']],
        output_core_dims=[['time']],
        vectorize=True,
        )
    ds_smooth['time'] = ds_raw['time']
    ds_smooth         = ds_smooth.transpose(*ds_raw.dims) # Make sure Dimensions Match Original...
    print("Smoothed in %.2fs" % (time.time()-st))
    return ds_smooth


# =============================================================================
#%% Part 1. Variable Loading and Preprocessing
# =============================================================================
stall             = time.time()

#%% 1.1 Set Variable and Experiment Loops (User Edits here)

lpf               = 6 # None to not LPF, otherwise cutoff month of LPF
vnames            = ["cre","tscre","ttcre","sst"]
compute_ids       = np.arange(0,3) #  Compute for all but sst

# Note: Important to Keep in this order, where 1-2 are full funs, 3-5 are 2055-2100 only...
expname           = "CERES_EBAF"

# Count Experients
nvars             = len(vnames)

# Output Path
outpath = "/home/niu4/gliu8/projects/ccfs/enso_ctone_fits/"

#%% 1.2 Load Variables

st        = time.time()
nvars     = len(vnames)
varsobs = []
for vv in range(nvars):
    vname = vnames[vv]
    if vname == "sst":
        expin = "ERA5"
    else:
        expin = "CERES_EBAF"
    ds    = ut.loadregrid(expin,vname,reformat=True)
    varsobs.append(ds)

varsobs = xr.merge(varsobs)
print("Loaded variables in %.2fs" % (time.time()-st))


#%% 1.3 Anomalize Each Separately

st           = time.time()
varsobs_anom = preproc_dataset(varsobs,vnames)
print("Preprocessed variables in %.2fs" % (time.time()-st))

#%% Optinally Apply Low PAss Filter
if lpf is not None:
    st           = time.time()
    varsobs_anom = pointwise_lp(varsobs_anom,lpf)
    print("Low-Pass Filtered in %.2fs" % (time.time()-st))

#%% 1.4 Calculate Nino3.4

bbox_nino34    = [-170+360,-120+360,-5,5]     # Nino3.4 Box
nino34         = proc.aavg(varsobs_anom['sst'],bbox_nino34)

# =============================================================================
#%% Part 2. Calculation of Regression Coefficients
# =============================================================================

# Loop by Variable
for ii in tqdm(range(len(compute_ids))):
    
    stl  = time.time()
    vv   = compute_ids[ii] # Choose Index of Varible to loop
    vname= vnames[vv]
    print("Starting calculations for %s..." % vname)
    
    # Get Variables
    anom_in = varsobs_anom[vname]
    ninoin  = nino34
    
    # Perform Fit
    fitout  = pointwise_fit_enso_ctone_phi(anom_in,ninoin,save_model=True)
    
    # Save Variable
    outname = "%sENSO_CTONE_PHI_Fits_ConcatExp_Obs_WholePeriod_%s.nc" % (outpath,vname)
    if lpf is not None:
        outname = proc.addstrtoext(outname,'lpf%0i' % lpf)
    fitout.to_netcdf(outname)

    
    print("\tCompleted calculations for %s in %.2fs." % (vname,time.time()-stl))

print("Script ran to completion in %.2fs." % (time.time()-stall))

