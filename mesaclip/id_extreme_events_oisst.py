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

#%% Helper Functions
    
def combine_consecutive_events(timeseries,event_indices,tol=1,verbose=True):
    # Given a tolerance level, merge consecutive events
    # Separate into discrete events
    # Copied from combine_events from enso event id
    nevents        = len(event_indices)
    if verbose:
        print("Original Starting Count: %i Events" % nevents)
        print("\tCombining Events <= [%i] timesteps apart" % tol)
    
    # Looping through events
    event_combine = []
    for ii in range(nevents+1): 
        
        if ii == (nevents): # Don't Perform Check for the last event
            # if verbose:
            #     print("Merging last event: %s" % event_merge)
            event_combine.append(event_merge) # This is defined later
            continue
        
        # Get the Index of the event
        ievent = event_indices[ii].item()
        
        if ii == 0: # Start Array with event to calculate Distances
            prev_id     = ievent
            event_merge = [ievent,]
            continue
        
        if (ievent - prev_id) <= tol: # Consecutive Event
            event_merge.append(ievent)
            # if verbose:
            #     print("%i is consecutive to previous events (%s)" % (ievent,event_merge))
        else: # Otherwise, just add event and merge
            event_combine.append(event_merge)
            event_merge = [ievent,] # Make a new one
            # if verbose:
            #     print("Making new event sequence at %i" % (ievent))
        prev_id = ievent
    nevents_combined = len(event_combine)
    if verbose:
        print("\tCulled to %i events!" % nevents_combined)
    return event_combine

def ds_dropvars(ds,keepvars):
    '''Drop variables in ds whose name is not in the list [keepvars]'''
    # Drop unwanted dimension
    dsvars = list(ds.variables)
    remvar = [i for i in dsvars if i not in keepvars]
    ds = ds.drop_vars(remvar)
    return ds

def retrieve_event_metrics_arr(event_combine,timeseries,tname='max'):
    """
    Processes merged events and calculates basic statistics (max, min, mean, stdev, duration, cumulative sum)
    
    Inputs
        event_combine : list of lists, where each element is an event and each inner list contains the indices corresponding to the event
        timeseries    : xr.DataArray , timeseries containing values

    Returns
        xr.DataSet Containing Time, Indices, and Summary Stats for each event, numbered by eventid.
    
    See `develop_MHW_id_code.ipynb` for debgging script

    """
    # Given a list of lists containing indices of combined events 
    nevents_combined = len(event_combine)
    
    # Part (1): Perform A Loop through Events =================================
    # Event Timing
    duration      = np.zeros((nevents_combined)) * np.nan
    
    # Indexing
    id_center     = duration.copy()
    id_first      = duration.copy()
    id_last       = duration.copy()
    id_max        = duration.copy()
    id_min        = duration.copy()
    
    # Event Stats
    event_mean    = duration.copy()
    event_std     = duration.copy()
    event_cumu    = duration.copy()
    for ie in range(nevents_combined):
        
        eventid_loop = event_combine[ie]
        
        # Get Variables  ------------------------------------------------
        intensities    = timeseries[eventid_loop]
        nconsecutive   = len(eventid_loop)
        
        # Record Some Metrics -------------------------------------------
        # Determine Some Indices for Metrics
        # Index within event chunk
        idmax            = np.argmax(np.abs(intensities))
        idmin            = np.argmin(np.abs(intensities))
        idfirst          = 0 #eventid_loop[0]
        idlast           = -1 #eventid_loop[-1]
        idcenter         = nconsecutive // 2 # Middle is just divided by 2
        if not nconsecutive & 0x1: # If Evenn, shift earlier
            idcenter     = idcenter - 1
        
        # Index from full timeseries
        id_center[ie]    = eventid_loop[idcenter]
        id_first[ie]     = eventid_loop[idfirst]
        id_last[ie]      = eventid_loop[idlast]
        id_max[ie]       = eventid_loop[idmax]
        id_min[ie]       = eventid_loop[idmin]
        
        # Timing (Note this saves to unreadable number...)
        duration[ie]     = nconsecutive            # Duration (Note assumes regular spacing...)
        
        # Statistics
        event_mean[ie]   = np.nanstd(intensities) # Mean
        event_std[ie]    = np.nanstd(intensities) # Standard Deviation
        event_cumu[ie]   = np.nansum(intensities) # Cumulative Values

    if tname == "max":
        id_out = id_max
    elif tname == "min":
        id_out = id_min
    elif tname == "center":
        id_out = id_center
    elif tname == "start":
        id_out = id_first
    elif tname == "end":
        id_out = id_last
    values_out = [timeseries[dd.astype(int)] for dd in id_out]
    #group2           = [duration,event_mean,event_std,event_cumu]
    
    return id_out,values_out,duration,event_mean,event_std,event_cumu

def get_rolling_threshold(timeseries,quantiles=[0.10,0.90],monthly=True):
    "Compute climatologically-varying percentile threshold and tile to original timeseries"
    
    if monthly: # Compute Quantiles Grouping by Month
        thres_bymon = timeseries.groupby('time.month').quantile(quantiles,dim='time')
    else:       # Compute Quantiles Grouping by Day of Year
        thres_bymon = timeseries.groupby('time.dayofyear').quantile(quantiles,dim='time')
    thresholds = []
    nq         = len(quantiles)
    for qq in range(nq):
        thres_in = thres_bymon.isel(quantile=qq)
        if monthly:
            thres = xr.ones_like(timeseries).groupby('time.month') * thres_in
        else:
            thres = xr.ones_like(timeseries).groupby('time.dayofyear') * thres_in
        thres = thres.drop_vars('quantile')
        thresholds.append(thres)
    
    thresholds = xr.concat(thresholds,dim='quantile')
    thresholds['quantile'] = quantiles
    return thresholds

    # # Combine into List and drop day of year
    # times_and_values = event_times + event_values
    # for dd in range(len(times_and_values)):
    #     if 'dayofyear' in list(times_and_values[dd].coords.keys()):
    #         times_and_values[dd] = times_and_values[dd].drop_vars('dayofyear')
    
def pad_nan(indata,nmax):
    ndata = len(indata)
    return np.pad(indata,(0,nmax-ndata),'constant',constant_values=np.nan)

def id_extremes_arr(timeseries,thres,positive,eventid_max=None,tol=1,verbose=False,tname='max'):
    if eventid_max is None:
        eventid_max = int(len(timeseries) * 0.25)
    # If NaN, just Continue
    if np.any(np.isnan(timeseries)):
        #print("Skipping NaN")
        dummy=np.zeros(eventid_max) * np.nan
        #output = *[dummy,]*5
        return dummy,dummy,dummy,dummy,0
    
    # Use Thresholds to find Events 
    if positive == True:
        below      = False
        event_indices = np.where(timeseries > thres)[0]
    
    else:
        below      = True
        # if verbose:
        #     print("Looking for events below threshold")
        event_indices = np.where(timeseries < thres)[0]
    
    nevents         = len(event_indices)
    # if verbose:
    #     print("Identified %i Events" % nevents)
    
    event_combine   = combine_consecutive_events(timeseries,event_indices,tol=tol,verbose=verbose)
    metrics_out     = retrieve_event_metrics_arr(event_combine,timeseries,tname=tname)
    # #id_out,values_out,duration,event_mean,event_std,event_cumu = metrics_out
    
    # For Output Variables, Pad with NaN #Enter the Amount
    nevents         = len(metrics_out[0])
    #npad            = eventid_max-nevents
    metrics_out     = [pad_nan(arr,eventid_max) for arr in metrics_out]
    id_out,values_out,duration,event_mean,event_std,event_cumu = metrics_out



    return id_out,values_out,duration,event_mean,nevents


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

tstart  = "1982-01-01"
tend    = "2025-12-31" #"2025-12-31"
outpath = "/home/niu4/gliu8/share/OISST/mergetest/"
expname = "oisst"
vname   = "sst"
freq    = "day_1"

deg          = 2

rawpath      = "/home/niu4/gliu8/share/OISST/mergetest/" #% (scenario)

tstart       = tstart.replace('-','')#npdatetime_to_str(dsanom.time[0]).replace('-','')
tend         = tend.replace('-','')  #npdatetime_to_str(dsanom.time[-1]).replace('-','')
outpath_proc = "%s/anom_detrend%i_%s-%s/" % (outpath,deg,tstart,tend,)

# Calculation Options
monthly    = True
tol        = 2
verbose    = False

# Make Output Directory
outdir_metrics = "%sMetrics_monthlybaseline%i_tol%02i_10to90Pct/" % (outpath,monthly,tol)
makedir(outdir_metrics)


if expname == "lores":
    enslist = np.arange(1,41,1)
else:
    enslist = np.arange(1,11,1)
nens = len(enslist)


#%% Ensemble Loop

start_all    = time.time()
outname      = "%s%s_%s_%s_anom.nc" % (outpath_proc,expname,freq,vname)

# Load the Variable
st           = time.time()
dsload       = xr.open_dataset(outname)
dsload       = dsload.load()['__xarray_dataarray_variable__']
print("Loaded in %.2fs" % (time.time()-st))


# Calculate Rolling Threshold (~40 sec), 134.29s on Niu
st           = time.time()
thresholds_global = get_rolling_threshold(dsload,quantiles=[0.10,0.90],monthly=True)
print("\tThreshold Calculated in %.2fs" % (time.time()-st))


# First, calculate for positive ===========================================
st         = time.time()
thresin    = thresholds_global.isel(quantile=1)
positive   = True
func_in     = lambda ds,thres,sign : id_extremes_arr(ds,thres,sign,tol=tol)
events_pos = xr.apply_ufunc(
    func_in,
    dsload,
    thresin,
    positive,
    input_core_dims=[["time"],["time"],[]],
    output_core_dims=[["eventid"],["eventid"],["eventid"],["eventid"],[]],
    vectorize=True,
)
print("\t(+) Events Found in %.2fs" % (time.time()-st))

# Postprocess Output
outnames = ['id_max','event_max','duration','event_mean','nevents']
dsout    = xr.merge([events_pos[ii].rename(outnames[ii]) for ii in range(len(events_pos))])

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
events_neg = xr.apply_ufunc(
    func_in,
    dsload,
    thresin,
    positive,
    input_core_dims=[["time"],["time"],[]],
    output_core_dims=[["eventid"],["eventid"],["eventid"],["eventid"],[]],
    vectorize=True,
)
print("\t(-) Events Found in %.2fs" % (time.time()-st))

# Postprocess Output
dsout    = xr.merge([events_neg[ii].rename(outnames[ii]) for ii in range(len(events_neg))])
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
    

    

    
    
    
    