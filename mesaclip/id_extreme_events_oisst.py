#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Copied from id_extreme_events_oisst_lores.py
    - Use output from anomalize_detrend_oisst.py

Created on Mon Aug  3 21:51:19 2026

@author: gliu

"""

import time
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import glob
import os

import sys

#%% Indicate Paths


mhwpath = "/home/niu4/gliu8/scripts/mhw_tools/"
sys.path.append(mhwpath)
import utils as mut

#%% Helper Functions

def ds_dropvars(ds,keepvars):
    '''Drop variables in ds whose name is not in the list [keepvars]'''
    # Drop unwanted dimension
    dsvars = list(ds.variables)
    remvar = [i for i in dsvars if i not in keepvars]
    ds = ds.drop_vars(remvar)
    return ds

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

def addstrtoext(name,addstr,adjust=0):
    """
    Add [addstr] to the end of a string with an extension [name.ext]
    Result should be "name+addstr+.ext"
    -4: 3 letter extension. -3: 2 letter extension
    """
    # Searches for 2 letter extension
    if (name[-2:] == "nc") or (name[-3] == "."): # Adjust for 2-letter extension
        print("2-letter extension detected")
        adjust = -1
    return name[:-(4+adjust)] + addstr + name[-(4+adjust):]

#%% User Edits

tstart              = "1982-01-01"
tend                = "2025-12-31" #"2025-12-31"
expname             = "oisst"
vname               = "sst"
freq                = "day_1" # "day_1"#"month_1"
procname            = "anom_detrend1_19820101-20251231_climatology1985to2014"
integratedintensity = False # Set to True to apply Integrated Intensity threshold

# Amplitude Thresholds
# See `calc_rolling_threshold_oisst.py
thresnc        = "/home/niu4/gliu8/projects/mesaclip/thresholds/%s/oisst_day_1_rolling_threshold_winsize15_pct010-090.nc" % procname# None
thresname      = "rolling15"

if freq=="month_1":
    thresnc   = None # No Rolling Threshold for Monthly Data
    thresname = "monthlythres"

# Duration Thresholds
combine_tol    = 2    # Set Fixed Tolerance (doesn't matter if efolding_tol is True)
efolding_tol   = False #False #False
winsize        = 15
if winsize == 0:
    efolding_nc    = "/home/niu4/gliu8/projects/mesaclip/memory/oisst_byday/daily_efolding_timescale_lagmax365_nowindow.nc"
else:
    efolding_nc    = "/home/niu4/gliu8/projects/mesaclip/memory/oisst_byday/daily_efolding_timescale_lagmax365_winsize15.nc"
efolding_vname = "efolding_timescale"

# see `visualize_efolding_timescales

# Paths
tstart       = tstart.replace('-','')#npdatetime_to_str(dsanom.time[0]).replace('-','')
tend         = tend.replace('-','')  #npdatetime_to_str(dsanom.time[-1]).replace('-','')
rawpath      = "/home/niu4/gliu8/share/OISST/mergetest/" #% (scenario)
outpath      = "/home/niu4/gliu8/share/OISST/mergetest/"

# Make Output Directory
if freq == "month_1":
    outpath = outpath + "monthly/"
outpath_proc = "%s/%s/" % (outpath,procname)
makedir(outpath_proc)

# Calculation Options
monthly    = True # Group Quantiles using Monthly Baseline
verbose    = False
if efolding_tol:
    print("Loading and using e-folding timescale tolerance")
    tolname="efolding_winsize%i" % winsize
    dstol = xr.open_dataset(efolding_nc)[efolding_vname].load()
    tol   = None
else:
    print("Using fixed tolerance: %i (%s)" % (combine_tol,freq))
    tol = combine_tol
    tolname="tol%02i" % combine_tol
    
if thresnc is not None:
    print("Loading custom threshold: %s" % thresname)
    dsthres = xr.open_dataset(thresnc).load()
    xrname  ='__xarray_dataarray_variable__'
    dsthres = dsthres[xrname].squeeze()
else:
    print("Regular Monthly Threshold will be used")

# Make Output Directory
outpath_event  = "/home/niu4/gliu8/projects/mesaclip/events/%s/" % (procname)
outdir_metrics = "%sMetrics_monthlybaseline%i_%s_%s_10to90Pct_integratedintensity%i/" % (outpath_event,monthly,thresname,tolname,integratedintensity)
makedir(outdir_metrics)

# Set Output Options for xrfunc, based on id_extremes_arr output
outdims_xrfunc = [["eventid"],
                  ["eventid"],
                  ["eventid"],
                  ["eventid"],
                  ['eventid'],
                  ['eventid'],
                  ['eventid'],
                  [],
                  ]

metric_names       = ['id_max',
                  'event_max',
                  'duration',
                  'event_mean',
                  'cumulative_intensity',
                  'id_first',
                  'id_last',
                  'nevents',
                  ]

#%% Start Script

start_all    = time.time()
loadname     = "%s%s_%s_%s_anom.nc"    % (outpath_proc,expname,freq,vname)
outname      = "%s%s_%s_%s_anom_%s.nc" % (outdir_metrics,expname,freq,vname,tolname)

# Load the Variable
st           = time.time()
dsload       = xr.open_dataset(loadname)
dsload       = dsload.load()['__xarray_dataarray_variable__']
print("Loaded in %.2fs" % (time.time()-st))

# Convert to No Leap
dsload       = dsload.convert_calendar('noleap')

# Calculate Rolling Threshold (~40 sec), 134.29s on Niu
# ~73 seconds for pre-loaded threshold
st                = time.time()
if thresnc is None:
    print("Calculating thresholds...")
    thresholds_global = mut.get_rolling_threshold(dsload,quantiles=[0.10,0.90],monthly=True)
else:
    print("Tiling existing thresholds")
    # Rename doy to dayofyear for groupby operation
    renamedict = dict(doy='dayofyear')
    dsthres    = dsthres.rename(renamedict)
    thresholds_global = xr.ones_like(dsload.squeeze()).groupby('time.dayofyear') * dsthres
print("\tThreshold Calculated in %.2fs" % (time.time()-st))


# Make Function to include combine tolerance
if not efolding_tol:
    func_in       = lambda ds,thres,sign : mut.id_extremes_arr(ds,thres,sign,tol=tol,integrated_threshold=integratedintensity)
else:
    func_in       = lambda ds,thres,sign,tolsel : mut.id_extremes_arr(ds,thres,sign,tol=tolsel,integrated_threshold=integratedintensity)

# Moved Above
#outdims_xrfunc = [["eventid"],["eventid"],["eventid"],["eventid"] ,['eventid']           ,[]]
#metric_names       = ['id_max'   ,'event_max','duration' ,'event_mean','cumulative_intensity','nevents']

# First, calculate for positive ===========================================
st         = time.time()
thresin    = thresholds_global.isel(quantile=1)
positive   = True
if efolding_tol:
    # Apply E-folding Threshold
    events_pos = xr.apply_ufunc(
        func_in,
        dsload,
        thresin,
        positive,
        dstol,
        input_core_dims=[["time"],["time"],[],['doy']],
        output_core_dims=outdims_xrfunc,
        vectorize=True,
    )
else:
    events_pos = xr.apply_ufunc(
        func_in,
        dsload,
        thresin,
        positive,
        input_core_dims=[["time"],["time"],[]],
        output_core_dims=outdims_xrfunc,
        vectorize=True,
    )
print("\t(+) Events Found in %.2fs" % (time.time()-st))

# Postprocess Output
dsout    = xr.merge([events_pos[ii].rename(metric_names[ii]) for ii in range(len(events_pos))])

# Reduce NaN
nmax = np.nanmax(dsout.nevents)
# Check to see that they are all NaN after the last event
chk_all_nan = [np.all(np.isnan(ds.isel(eventid=nmax).data)) for ds in events_pos[:-1]]
if np.all(chk_all_nan):
    print("\tReducing to Last Event ID (i=%i)" % (nmax))
    dsout = dsout.isel(eventid=slice(0,nmax))

# Save Positive
st          = time.time()
outname_pos = addstrtoext(outname,"_positive_events")
dsout.to_netcdf(outname_pos)
print("\tSaved (+) Events in %.2f" % (time.time()-st))
    
del events_pos,dsout

# Now Calculate for Negative ===========================================
st         = time.time()
thresin    = thresholds_global.isel(quantile=0)
positive   = False
if efolding_tol:
    events_neg = xr.apply_ufunc(
        func_in,
        dsload,
        thresin,
        positive,
        dstol,
        input_core_dims=[["time"],["time"],[],['doy']],
        output_core_dims=outdims_xrfunc,
        vectorize=True,
    )
else:
    events_neg = xr.apply_ufunc(
        func_in,
        dsload,
        thresin,
        positive,
        input_core_dims=[["time"],["time"],[]],
        output_core_dims=outdims_xrfunc,
        vectorize=True,
    )
print("\t(-) Events Found in %.2fs" % (time.time()-st))

# Postprocess Output
dsout    = xr.merge([events_neg[ii].rename(metric_names[ii]) for ii in range(len(events_neg))])
# Reduce NaN
nmax     = np.nanmax(dsout.nevents)
# Check to see that they are all NaN after the last event
chk_all_nan = [np.all(np.isnan(ds.isel(eventid=nmax).data)) for ds in events_neg[:-1]]
if np.all(chk_all_nan):
    print("\tReducing to Last Event ID (i=%i)" % (nmax))
    dsout = dsout.isel(eventid=slice(0,nmax))
    
# Save Netat9ve
st          = time.time()
outname_neg = addstrtoext(outname,"_negative_events")
dsout.to_netcdf(outname_neg)
print("\tSaved (-) Events in %.2f" % (time.time()-st))

print("Completed n %.2fs" % (time.time()-start_all))
del events_neg,dsout,dsload,thresholds_global
# = addstrtoext(outname,"_")

    