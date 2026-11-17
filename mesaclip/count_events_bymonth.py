#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Script to Count MHW Events by Month
Also Identify Month of Maximum

Created on Mon Nov 16 14:06:58 2026

@author: gliu
"""

# import sys
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
from datetime import datetime

from tqdm import tqdm

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

import metpy.calc as mpcalc

#%% Functions

# Inputs
def count_months_by_event(startid,endid,times,verbose=False):
    # Limit to Event Only and Convert to Int
    ntotal  = len(startid)
    nevents = len(startid[~np.isnan(startid)])
    sid     = startid[:nevents].astype(int)
    eid     = endid[:nevents].astype(int)
    
    # Quick Check to make sure sid and eid are the same
    nevents_end = len(endid[~np.isnan(endid)])
    if (nevents_end != nevents) and verbose:
        print("WARNING! # of ending indices (%i) does not match starting indices (%i). Errors may result..." % (nevents_end,nevents))
    
    
    # For Each Event, Determine the count
    month_by_event = np.zeros((12,ntotal)) * np.nan
    for ee in range(nevents):
        
        # Get Months of Event Between Start and End Date
        eventids    = np.arange(sid[ee],eid[ee]+1)
        eventdates  = times[eventids]
        eventmonths = np.array([tt.month for tt in eventdates]) #eventdates.month
        
        # Get Counts and write to array
        emons,ecounts = np.unique(eventmonths,return_counts=True)
        for nm in range(len(emons)): #emon in emons:
            emon  = emons[nm]
            count = ecounts[nm]
            month_by_event[emon-1,ee] = count
    return month_by_event


def calculate_save_maxmon(ncname,dstime_daily):
    
    # Load Output 
    st            = time.time()
    dsin          = xr.open_dataset(ncname).load()
    print("Data Loaded in %.2fs" % (time.time()-st))


    #%% Count Months for each Event

    st            = time.time()
    countens = xr.apply_ufunc(
            count_months_by_event,# Pass the function
            dsin.id_first.squeeze(),
            dsin.id_last.squeeze(),
            dstime_daily,
            input_core_dims=[['eventid'],['eventid'],['time']],
            output_core_dims=[['month','eventid'],],  # Output Dimension
            vectorize=True,  # True to loop over non-core dims
        )
    print("Counted Monthly Frequency for each event in %.2fs" % (time.time()-st))

    #% Calculate and Save Max Count

    st             = time.time()
    countens_nonan = xr.where(np.isnan(countens),0,countens) # Replace NaN with 0 to prevent errors
    maxcountmon    = countens_nonan.reduce(np.nanargmax,dim='month') + 1
    maxcountmon    = xr.where(np.isnan(dsin.id_first.squeeze()),np.nan,maxcountmon)
    print("Determined Max Month for each event in %.2fs" % (time.time()-st))

    #% Save the Output


    st             = time.time()
    outname_counts = proc.addstrtoext(ncname,"_CountByMonth")
    outname_maxmon = proc.addstrtoext(ncname,"_CountMaxMonth")
    countens.to_netcdf(outname_counts)
    maxcountmon.to_netcdf(outname_maxmon)
    print("Saved Output in %.2fs" % (time.time()-st))

    return countens,maxcountmon

def calculate_save_maxmon_dummy(ncname,dstime_daily):
    # For Checking Names, Delete Later
    # Load Output 
    st            = time.time()
    # dsin          = xr.open_dataset(ncname).load()
    # print("Data Loaded in %.2fs" % (time.time()-st))


    # #% Count Months for each Event

    # st            = time.time()
    # countens = xr.apply_ufunc(
    #         count_months_by_event,# Pass the function
    #         dsin.id_first.squeeze(),
    #         dsin.id_last.squeeze(),
    #         dstime_daily,
    #         input_core_dims=[['eventid'],['eventid'],['time']],
    #         output_core_dims=[['month','eventid'],],  # Output Dimension
    #         vectorize=True,  # True to loop over non-core dims
    #     )
    # print("Counted Monthly Frequency for each event in %.2fs" % (time.time()-st))

    # #% Calculate and Save Max Count

    # st             = time.time()
    # countens_nonan = xr.where(np.isnan(countens),0,countens) # Replace NaN with 0 to prevent errors
    # maxcountmon    = countens_nonan.reduce(np.nanargmax,dim='month') + 1
    # maxcountmon    = xr.where(np.isnan(dsin.id_first.squeeze()),np.nan,maxcountmon)
    # print("Determined Max Month for each event in %.2fs" % (time.time()-st))

    # #% Save the Output


    # st             = time.time()
    outname_counts = proc.addstrtoext(ncname,"_CountByMonth")
    outname_maxmon = proc.addstrtoext(ncname,"_CountMaxMonth")
    print(outname_counts)
    print(outname_maxmon)
    print("\n")
    # countens.to_netcdf(outname_counts)
    # maxcountmon.to_netcdf(outname_maxmon)
    # print("Saved Output in %.2fs" % (time.time()-st))

    return None,None
    


#%% User Edits

datpath                        = "/home/niu4/gliu8/projects/mesaclip/events/anom_detrend2_19820101-20251231"
sign                           = "positive" # Just Load 1 Sign for Now
cumulative_intensity_criterion = np.arange(4)


# Other Necessary Variables/Paths

# Date Time (MAKE SURE IS NO LEAP--must match with original indices from id_events_oisst....etc)
dstime_daily   = xr.date_range('1982-01-01','2025-12-31',freq="D",use_cftime=True,calendar="noleap")
dstime_monthly = xr.date_range('1982-01-01','2025-12-31',freq="MS",use_cftime=True,calendar="noleap")

figpath        = "/home/niu4/gliu8/figures/bydate/2026-11-17-CGC-Baseline_update/"
proc.makedir(figpath)

#procname      = "winsize0_monthres"
winsize        = 15
min_duration   = 5

#%% Load data

# # Start by Testing with the hobday Defi
# procpath      = "Metrics_monthlybaseline1_rolling15_tol02_10to90Pct_cumulativeintensity0_BaselineOnly"
# ncname        = "%s/%s/oisst_day_1_sst_anom_tol02_%s_events.nc" % (datpath,procpath,sign)

# _,_ = calculate_save_maxmon(ncname,dstime_daily)

#%% Now Also Process Cumulative Indices

# Just Load Positive Events

baseline_only=True
for cc in tqdm(range(4)):
    stall = time.time()
    if baseline_only and cc > 0:
        procpath = "Metrics_monthlybaseline1_rolling15_efolding_winsize15_10to90Pct_cumulativeintensity%i_BaselineOnly" % (cc)
    else:
        procpath = "Metrics_monthlybaseline1_rolling15_efolding_winsize15_10to90Pct_cumulativeintensity%i" % (cc)
    ncname = "%s/%s/oisst_day_1_sst_anom_efolding_winsize%02i_%s_events.nc" % (datpath,procpath,winsize,sign)
    
    _,_, = calculate_save_maxmon(ncname,dstime_daily)
    
    print("Processed Set %i in %.2fs" % (cc,time.time()-stall))
    
    

