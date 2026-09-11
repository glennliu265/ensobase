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
    # Make some adjustments based on the tolerance
    # Note: Asssumes timeseries has an even number of years or is evenly divisible...
    if type(tol) != int: # Tolerance is NOT just 1 value
        if len(tol) == 12: # Monthly E-Folding Tolerance
            # Tolerance Applied by Month
            nyrs     = int(len(timeseries)/12)
            tolcheck = np.tile(np.arange(1,13,1),nyrs)
            if len(timeseries)%12 != 0:
                print("Warning, Timeseries is not evenly divisible by 12, indexing errors may result.")
                print(tolcheck.shape)
                print(timeseries.shape)
        elif len(tol) == 365: # Daily E-folding Tolerance
            # Tolerance Applied by Day of Year
            nyrs     = int(len(timeseries)/365)
            tolcheck = np.tile(np.arange(1,366,1),nyrs)
            if len(timeseries)%365 != 0:
                print("Warning, Timeseries is not evenly divisible by 365, indexing errors may result.")
                print(tolcheck.shape)
                print(timeseries.shape)
        else:
            print("Invalid size for Tol... must be either a constant or array of 12 (monthly) or 365 (doy)")
    else:
        tolcheck = None
    
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

        if tolcheck is None: # Tolerance is just a single value
            tol_in = tol
        else: # Tolerance is an Array (day of year, or month)
            event_category = tolcheck[ievent] # Check day or month of event
            tol_in         = tol[event_category-1]    # -1 for Python Indexing
            
        
        if (ievent - prev_id) <= tol_in: # Consecutive Event
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

    

    return id_out,values_out,duration,event_mean,event_cumu,nevents


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

tstart         = "1982-01-01"
tend           = "2025-12-31" #"2025-12-31"
expname        = "oisst"
vname          = "sst"
freq           = "day_1"#"month_1"
deg            = 2 # Detrend Degree

# Amplitude Thresholds
# See `calc_rolling_threshold_oisst.py`
thresnc        = "/home/niu4/gliu8/projects/mesaclip/thresholds/anom_detrend2_19820101-20251231/oisst_day_1_rolling_threshold_winsize15_pct010-090.nc" # None
thresname      = "rolling15"

# Duration Thresholds
combine_tol    = 2 # Set Fixed Tolerance (doesn't matter if efolding_tol is True)
efolding_tol   = False #False
winsize        = 15
if winsize == 0:
    efolding_nc    = "/home/niu4/gliu8/projects/mesaclip/memory/oisst_byday/daily_efolding_timescale_lagmax365_nowindow.nc"
else:
    efolding_nc    = "/home/niu4/gliu8/projects/mesaclip/memory/oisst_byday/daily_efolding_timescale_lagmax365_winsize15.nc"
efolding_vname = "efolding_timescale"

# see `visualize_efolding_timescales

# Paths
rawpath      = "/home/niu4/gliu8/share/OISST/mergetest/" #% (scenario)
outpath      = "/home/niu4/gliu8/share/OISST/mergetest/"

# Make Output Directory
tstart       = tstart.replace('-','')#npdatetime_to_str(dsanom.time[0]).replace('-','')
tend         = tend.replace('-','')  #npdatetime_to_str(dsanom.time[-1]).replace('-','')
if freq == "month_1":
    outpath_proc = "%s/monthly/anom_detrend%i_%s-%s/" % (outpath,deg,tstart,tend,)
else:
    outpath_proc = "%s/anom_detrend%i_%s-%s/" % (outpath,deg,tstart,tend,)




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
outpath_event = "/home/niu4/gliu8/projects/mesaclip/events/anom_detrend%i_%s-%s/" % (deg,tstart,tend,)
outdir_metrics = "%sMetrics_monthlybaseline%i_%s_%s_10to90Pct/" % (outpath_event,monthly,thresname,tolname)
makedir(outdir_metrics)

# Set Number of Ensembles (only Relevant for MESACLIP)
if expname == "lores":
    enslist = np.arange(1,41,1)
else:
    enslist = np.arange(1,11,1)
nens = len(enslist)

#%% Ensemble Loop

start_all    = time.time()
loadname     = "%s%s_%s_%s_anom.nc" % (outpath_proc,expname,freq,vname)

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
    thresholds_global = get_rolling_threshold(dsload,quantiles=[0.10,0.90],monthly=True)
else:
    print("Tiling existing thresholds")
    # Rename doy to dayofyear for groupby operation
    renamedict = dict(doy='dayofyear')
    dsthres    = dsthres.rename(renamedict)
    thresholds_global = xr.ones_like(dsload.squeeze()).groupby('time.dayofyear') * dsthres
print("\tThreshold Calculated in %.2fs" % (time.time()-st))


# Make Function to include combine tolerance
if not efolding_tol:
    func_in       = lambda ds,thres,sign : id_extremes_arr(ds,thres,sign,tol=tol)
else:
    func_in       = lambda ds,thres,sign,tolsel : id_extremes_arr(ds,thres,sign,tol=tolsel)


outdims_xrfunc = [["eventid"],["eventid"],["eventid"],["eventid"],['eventid'],[],]
outnames       = ['id_max','event_max','duration','event_mean','cumulative_intensity','nevents']

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

    