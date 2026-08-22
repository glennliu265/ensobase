#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate the Seasonal Fits to Nino3.4 * cos(wt*phi) for a selected variable
for Early/Late SSP.585 Periods

Copied from `compare_sst_cre_ssp585_awi.ipynb`

Created on Fri Aug 21 15:17:43 2026

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

# =============================================================================
#%% Part 1. Variable Loading...
# =============================================================================
stall             = time.time()

#%% 1.1 Set Variable and Experiment Loops (User Edits here)

timerange_early   = ['2015-01-01','2050-12-31']
timerange_late    = ['2065-01-01','2100-12-31']

vnames            = ["cre","tscre","ttcre","sst"]
compute_ids       = np.arange(0,3) #  Compute for all but sst

# Note: Important to Keep in this order, where 1-2 are full funs, 3-5 are 2055-2100 only...
expnames  = ["TCo319-DART-ssp585d-gibbs-charn",
             "TCo319_ssp585",
             "TCo319_ssp585_ens01",
             "TCo319_ssp585_ens02",
             "TCo319_ssp585_ens03",]

# Count Experients
nexps             = len(expnames)
nvars             = len(vnames)

# Designate_Names
earlynames        = ["%s_early" % expnames[ex] for ex in np.arange(2)]
latenames         = ["%s_late" % expnames[ex] for ex in np.arange(5)]

earlyperiodname   = "%s-%s" % (timerange_early[0][:4],timerange_early[1][:4])
lateperiodname    = "%s-%s" % (timerange_late[0][:4],timerange_late[1][:4])
print(earlyperiodname)
print(lateperiodname)

# Output Path
outpath = "/home/niu4/gliu8/projects/ccfs/enso_ctone_fits/"

#%% 1.2 Load Variables and Merge into Dataset

st        = time.time()
varsbyexp = []
for ex in tqdm(range(nexps)):
    expname = expnames[ex]

    vbv = []
    for vv in range(nvars):
        vname = vnames[vv]
        ds    = ut.loadregrid(expname,vname,reformat=True)
        vbv.append(ds)
    varsbyexp.append(vbv)
varsawi  = [xr.merge(ds) for ds in varsbyexp]
print("Loaded variables in %.2fs" % (time.time()-st))

#%% 1.3 Separate into Early and Late Periods

# Early Period (Only 2 Runs)
vars_awi_early = [varsawi[ex].sel(time=slice(*timerange_early)) for ex in np.arange(0,2)]

# Late Period (All 4 Runs)
vars_awi_late  = [varsawi[ex].sel(time=slice(*timerange_late)) for ex in np.arange(5)]

# =============================================================================
#%% Part 2. Preprocessing
# =============================================================================

#%% 1.1 Deseason and Detrend

st = time.time()
varsanom_awi_early = [preproc_dataset(ds,vnames) for ds in vars_awi_early]
varsanom_awi_late  = [preproc_dataset(ds,vnames) for ds in vars_awi_late]
print("Preprocessed variables in %.2fs" % (time.time()-st))

#%% 1.2 Calculate Niño3.4

# Briefly Compute Nino3.4 for each base
bbox_nino34  = [-170+360,-120+360,-5,5]     # Nino3.4 Box
nino34_early = [proc.aavg(ds['sst'],bbox_nino34) for ds in varsanom_awi_early]
nino34_late  = [proc.aavg(ds['sst'],bbox_nino34) for ds in varsanom_awi_late]

# =============================================================================
#%% Part 3. Calculation of Regression Coefficients
# =============================================================================

# Loop by Variable
for ii in range(len(compute_ids)):
    stl  = time.time()
    vv   = compute_ids[ii] # Choose Index of Varible to loop
    vname= vnames[vv]
    print("Starting calculations for %s..." % vname)
    
    
    # First, Compute for Late Period (5 experiments)
    for ex in tqdm(range(5)):
        
        expname = latenames[ex]
        
        anom_in = varsanom_awi_late[ex][vname]
        ninoin  = nino34_late[ex]
        
        fitout  = pointwise_fit_enso_ctone_phi(anom_in,ninoin,save_model=True)
        
        outname = "%sENSO_CTONE_PHI_Fits_%s_%s_%s.nc" % (outpath,expname,lateperiodname,vname)
        fitout.to_netcdf(outname)
    
    # Next, Loop for Early Period (2 experiments)
    for ex in tqdm(range(2)):
        
        expname = earlynames[ex]
        
        anom_in = varsanom_awi_early[ex][vname]
        ninoin  = nino34_early[ex]
        
        fitout  = pointwise_fit_enso_ctone_phi(anom_in,ninoin,save_model=True)
        
        outname = "%sENSO_CTONE_PHI_Fits_%s_%s_%s.nc" % (outpath,expname,earlyperiodname,vname)
        fitout.to_netcdf(outname)
    
    print("\tCompleted calculations for %s in %.2fs." % (vname,time.time()-stl))

print("Script ran to completion in %.2fs." % (time.time()-stall))
     
        
        
        
        
    
    




#%%

#%%