#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

In SSP5.85 Simulation, Estimate Radiative Kernels using a sliding window


(1) Load CCFs and TOA Fluxes (Raw)
(2) Subset into Windows
(3) Looping for each window...
    (a) Deseason and Detrend
    (b) Estimate Radiative Kernel (fit)
(4) Package output and save...


Scripts Copied from
    - sliding_spectra_ccf.py
    
    

Created on Tue May 19 13:57:06 2026

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
from tqdm import tqdm

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut


#%% User Edits

# Expeirment Selections
expname      = "TCo319_ssp585"
ccf_vars     = ["sst","eis","Tadv","r700","w700","ws10"]
tstart       = '2015-01-01'
tend         = '2100-12-31'

rawpath      = "/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/TCo319_ssp585/raw/"


# Analysis Options
bbox_sep     = [-90+360,-75+360,-40,-15] # Southeast Tropical Pacific Box from Kang et al. 2026
bbsel        = bbox_sep
bbname       = "SEP"
load_global  = True


nyr_sliding  = 30 # Select Sliding Window Length
flxname      = "tscre"

# Output Options
outpath = "/home/niu4/gliu8/projects/ccfs/metrics/regrid_1x1/scrap/sliding_kernels/"
outname = "%sSlidingKernels_%s_%s_AllMonth_%02iyrwindow.nc" % (outpath,expname,flxname,nyr_sliding,)


# ============================================
#%% Part 1. Load (raw) CCFs for an experiment
# ============================================







ds_ccf       = []
nccfs        = len(ccf_vars)
for cc in range(nccfs):
    vname   = ccf_vars[cc]
    ncname  = "%s%s.nc" % (rawpath,vname)
    ds      = xr.open_dataset(ncname)[vname]#.load()
    if load_global:
        dsreg = ds.load()
        bbname = "Global"
    else:
        dsreg   = proc.sel_region_xr(ds,bbsel).load()
    dsreg   = ut.standardize_names(dsreg)
    dsreg   = dsreg.sel(time=slice(tstart,tend))
    print("%s" % vname)
    print("\tStart time is %s" % dsreg.time.isel(time=0))
    print("\tEnd time is %s" % dsreg.time.isel(time=-1))
    
    ds_ccf.append(dsreg)

# ============================================
#%% Part 2. Load (raw) fluxes for an experiment
# ============================================

flxnames = ["cre","tscre","ttcre"]
ds_flx       = []
for ff in range(3):
    vname   = flxnames[ff]
    ncname  = "%s%s.nc" % (rawpath,vname)
    ds      = xr.open_dataset(ncname)[vname]
    
    if load_global:
        dsreg = ds.load()
    else:
        dsreg      = proc.sel_region_xr(ds,bbsel).load()
    dsreg      = ut.standardize_names(dsreg)
    dsreg      = ut.varcheck(dsreg,vname,expname)
    
    dsreg   = dsreg.sel(time=slice(tstart,tend))
    
    dsflx   = dsreg.copy()
    ds_flx.append(dsreg)

# Select TOA Flux
ffsel   = flxnames.index(flxname)
flxin   = ds_flx[ffsel]

# ============================================
#%% Part 3. Estimate Radiative Kernels
# ============================================

# Subset into periods and preprocess


# Generate Periods for CCFs, then preprocess
st           = time.time()
subsets_ccfs = [ut.generate_periods(ds,nyr_sliding)[0] for ds in ds_ccf]
ccf_anoms    = [ut.preprocess_byperiod(ds) for ds in subsets_ccfs]
proc.printtime(st,"Anomalized CCF by period")

# Generate Periods for flux
st           = time.time()
subsets_flx  = ut.generate_periods(flxin,nyr_sliding)[0]
flx_anom     = ut.preprocess_byperiod(subsets_flx)
proc.printtime(st,"Anomalized Flux by period")

#% Calculate the radiative kernels
nperiods = len(flx_anom)
mlrout_byperiod = []
for nw in tqdm(range(nperiods)):
    
    # Get CCFs and reshape into [time x predictor x lat x lon]
    ccf_period = [cc[nw] for cc in ccf_anoms] #ccf_anoms[]
    flx_period = flx_anom[nw]
    
    # Reshape CCFs and merge along predictor dimension
    ccf_names  = [cc.name for cc in ccf_period]
    
    # Reassign Time to Make Sure they all match (honestly can do this before resampling..)
    match_flx_time   = lambda ds : proc.match_time_month(ds,flx_period)[0]
    ccf_period       = [match_flx_time(ds) for ds in ccf_period]
    # Replace time If Needed
    timeref    = flx_period.time
    for cc in range(nccfs):
        
        timeccf = ccf_period[cc].time
        chk     = timeref.data == timeccf.data
        if np.any(~chk):
            print("Warning, replacing time for CCF %s" % ccf_vars[cc])
            ccf_period[cc]['time'] = timeref
            
    # Now Merge Things
    ccf_period_merge = xr.concat(ccf_period,dim='predictors')
    ccf_period_merge = ccf_period_merge.transpose('time','predictors','lat','lon')
    
    # Perform pointwise MLR
    mlrout = proc.pointwise_mlr(ccf_period_merge,flx_period,
                                fill_value=0,standardize=True,
                                predictor_names=ccf_names)
    
    mlrout_byperiod.append(mlrout)

#%% Before Appending, Need to fix the time dimension
# Assign shared "Time Index" and store actual time ranges in another Window

time_index      = np.arange(len(mlrout_byperiod[nw].time))
time_byperiod   = []
for nw in tqdm(range(nperiods)):
    
    ptime = mlrout_byperiod[nw].time.copy()
    time_byperiod.append(ptime)
    mlrout_byperiod[nw]['time'] = time_index

time_byperiod    = np.array(time_byperiod)
coords_time      = dict(period=np.arange(nperiods),time=time_index)
ds_time_byperiod = xr.DataArray(time_byperiod,name='timestamp',
                                coords=coords_time,dims=coords_time)

    
# Now Concatenate By Period
ds_out = xr.concat(mlrout_byperiod,dim='period')
ds_out = xr.merge([ds_out,ds_time_byperiod])

#%% Save Output


st = time.time()
edict = proc.make_encoding_dict(ds_out)
ds_out.to_netcdf(outname,encoding=edict)
proc.printtime(st,"Saved output")



#%% Below is Scrap. I can Delete
# #%% Perform period-wise estimate of radiative Kernels

# # Subset CCFs
# aavgs_ccfs   = [proc.area_avg_cosweight(ds) for ds in ds_ccf]
# nyr_ceres    = 24
# subsets_ccfs = [ut.generate_periods(ds,nyr_ceres)[0] for ds in aavgs_ccfs]
# st = time.time()
# ccf_anoms = [ut.preprocess_byperiod(ds) for ds in subsets_ccfs]

# proc.printtime(st,"Preprocessed CCFs")

# # Subset Fluxes
# aavgs_flxs   = [proc.area_avg_cosweight(ds) for ds in ds_flx]
# subsets_flxs = [ut.generate_periods(ds,nyr_ceres)[0] for ds in aavgs_flxs]
# st = time.time()
# flx_anoms    = [ut.preprocess_byperiod(ds) for ds in subsets_flxs]
# proc.printtime(st,"Preprocessed Fluxes")

# # Calculate the radiative kernels
# nw     = 0
# ff     = 1

# # Need to change CCFs from [ccf][period][time] to [period][ccf x time]
# ccf_anoms_arr = np.array(ccf_anoms)

# def reshape_ccfs_predictor(ccf_anoms,nw):
#     # ccf_anoms is [ccf][period][time]
#     # nw is the index of the period
#     # takes time dim of first ccf (note that xr.concat yielded non-matching times...)
#     # Standardizes along the predictor dimension
#     ccf_names        = [cc[nw].name for cc in ccf_anoms]
#     ccf_window       = np.array([cc[nw].data for cc in ccf_anoms])
#     coords           = dict(predictor=ccf_names,time=ccf_anoms[0][nw].time)
#     ccf_window       = xr.DataArray(ccf_window,dims=coords,coords=coords)#xr.concat(ccf_window,dim='predictor')
#     ccf_window       = ccf_window / ccf_window.std('time')
#     return ccf_window

# nperiods  = len(ccf_anoms[0])
# coeffs    = []
# r2s       = []
# residuals = []

# for nw in tqdm(range(nperiods)):
    
#     ccfanoms_in   = reshape_ccfs_predictor(ccf_anoms,nw).transpose('time','predictor')
#     mlrout        = proc.mlr_point(ccfanoms_in,flx_anoms[ff][nw])
    
#     #[cc.time.isel(-1) for cc in ccf_anoms]
#     coeffs.append(mlrout[0])
#     r2s.append(mlrout[1])
#     residuals.append(mlrout[2])

# # Get Tcenters for plotting
# _,tranges = ut.generate_periods(aavgs_ccfs[0],nyr_ceres)
# tcenters  = [ut.get_center_time(tt) for tt in tranges]

# # Calxulate RMSE
# rmse_byperiod = []
# for nw in range(nperiods):
#     rmse = np.nanmean((flx_anoms[ff][nw] - residuals[nw])**2)**(0.5)
#     rmse_byperiod.append(rmse.item())
# rmse_byperiod = np.array(rmse_byperiod)
    

