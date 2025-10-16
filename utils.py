#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Shared functions used by scripts in ensobase repository

Currently runs on niu (need to load [amv] dependencies...)
add this line to the top of the script to import module:

----
ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut
----

Function Types (will reorganize later):
    (A)     : AWI-CM3 specific preprocessing
    (c)     : calculations/analysis
    (g)     : general/universal function, can move to proc
    (l)     : loader
    (v)     : visualization

Function                Description
--------                -----------
awi_mean_loader         : (l) load mean/monvar/scycle calculations from calc_mean_patterns_TP 
calc_lag_regression_1d  : (g) Compute lead lag regression for 1d timeseries
combine_events          : (g) Given identified events, combine similar events and get other traits (duration, etc)
init_tp_map             : (v) initialize tropical Pacific plot 
load_ensoid             : (l) load enso indices calculated by calc_nino34.py
mcsample                : (g) Monte Carlo Sampler to repeat function
preprocess_enso         : (c) detrend (quadratic) and deseasonalize for ENSO calculations
remove_duplicate_times  : (g) Remove duplicate times from a DataArray
swap_rename             : (g) check if variable exists and rename if so
stack_events            : (g) For 1-d timeseries, stack events along specified leads/lags
stack_events_2d         : (g) Stack Events, but applied to 2D case with time x lat x lon....
standardize_names       : (A) uses swap_rename to replace variable and dimension names in AWI_CM3 output 
varcheck                : (A) checks and converts variables for AWI-CM3

Created on Wed Oct  8 15:26:57 2025

@author: gliu
"""

import numpy as np
import xarray as xr
import cartopy
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import time

#%%

amvpath = "/home/niu4/gliu8/scripts/commons"
sys.path.append(amvpath)
from amv import proc,viz

#%%

def awi_mean_loader(expname,vname,calcname,outpath=None):
    if outpath is None:
        outpath = "/home/niu4/gliu8/projects/scrap/TP_crop/summary/"
    ncname = "%s%s_%s_%s.nc" % (outpath,expname,vname,calcname)
    ds = xr.open_dataset(ncname).load()
    return ds

def init_tp_map(nrow=1,ncol=1,figsize=(12.5,4.5),ax=None):
    bbplot = [120, 290, -20, 20]
    fix_lon = np.hstack([np.arange(120,190,10),np.arange(-180,-60,10)])
    proj   = ccrs.PlateCarree(central_longitude=180)
    projd  = ccrs.PlateCarree()
    
    if ax is None:
        fig,axs = plt.subplots(nrow,ncol,figsize=figsize,subplot_kw={'projection':proj})
        newfig = True
    else:
        newfig = False
    if nrow != 1 or ncol != 1:
        for ax in axs.flatten():
            ax.set_extent(bbplot)
            ax     = viz.add_coast_grid(ax,bbox=bbplot,fill_color='k',
                                        proj=ccrs.PlateCarree(),fix_lon=fix_lon,ignore_error=True)
        ax = axs
    else:
        ax = axs
        ax.set_extent(bbplot)
        ax     = viz.add_coast_grid(ax,bbox=bbplot,fill_color='k',
                                    proj=ccrs.PlateCarree(),fix_lon=fix_lon,ignore_error=True)

    if newfig:
        return fig,ax
    return ax
    
def calc_lag_regression_1d(var1,var2,lags): # CAn make 2d by mirroring calc_lag_covar_annn
    # Calculate the regression where
    # (+) lags indicate var1 lags  var2 (var 2 leads)
    # (-) lags indicate var1 leads var2 (var 1 leads)
    
    ntime = len(var1)
    betalag = []
    poslags = lags[lags >= 0]
    for l,lag in enumerate(poslags):
        varlag   = var1[lag:]
        varbase  = var2[:(ntime-lag)]
        
        # Calculate correlation
        beta = np.polyfit(varbase,varlag,1)[0]   
        betalag.append(beta.item())
    
    neglags = lags[lags < 0]
    neglags_sort = np.sort(np.abs(neglags)) # Sort from least to greatest #.sort
    betalead = []
    
    for l,lag in enumerate(neglags_sort):
        varlag   = var2[lag:] # Now Varlag is the base...
        varbase  = var1[:(ntime-lag)]
        # Calculate correlation
        beta = np.polyfit(varlag,varbase,1)[0]   
        betalead.append(beta.item())
    
    # Append Together
    return np.concat([np.flip(np.array(betalead)),np.array(betalag)])


def combine_events(var_in,id_in,tol=1,verbose=True):
    
    # Separate into discrete events
    nevents        = id_in.data.sum().item()
    eventids       = np.where(id_in)[0]
    event_combine = []
    for ii in range(nevents+1): #event_id:
        if ii == (nevents):
            if verbose:
                print("Merging last event: %s" % event_merge)
            event_combine.append(event_merge)
            continue
        
        ievent = eventids[ii].item()#.data
        
        if ii == 0:
            prev_id     = ievent
            event_merge = [ievent,]
            continue
        
        if (ievent - prev_id) <= tol: # Consecutive Event
            event_merge.append(ievent)
            if verbose:
                print("%i is consecutive to previous events (%s)" % (ievent,event_merge))
        else: # Otherwise, just add event and merge
            event_combine.append(event_merge)
            event_merge = [ievent,] # Make a new one
            if verbose:
                print("Making new event sequence at %i" % (ievent))
        prev_id = ievent
    if verbose:
        print("Identified %i events!" % len(event_combine))

    # Identify Event Centers
    ncevent         = len(event_combine)
    event_time      = []
    center_ids      = []
    event_max       = []
    for ie in range(ncevent):
        idsin = event_combine[ie]
        if len(idsin) == 1: # Only 1 step
            event_time.append(var_in.time.isel(time=idsin[0]))
            event_max.append(var_in.isel(time=idsin[0]).data.item())
            center_ids.append(idsin[0])
        else:
            amplitudes = var_in.isel(time=idsin) #.argmax()
            idmax      = np.argmax(np.abs(amplitudes.data)).item()
            event_max.append(np.nanmax(np.abs(amplitudes.data)).item())
            event_time.append(var_in.time.isel(time=idsin[idmax]))
            center_ids.append(idsin[idmax])

    def get_mon(ds):
        return [da.time.dt.month.data.item() for da in ds]
    
    def get_duration(ds):
        return [len(ts) for ts in ds]
    
    durations = get_duration(event_combine)
    months    = get_mon(event_time)
    
    outdict = dict(event_time=event_time,
                   center_ids=center_ids,
                   event_combine=event_combine,
                   event_max=event_max,
                   durations=durations,
                   eventmonths=months)
    return outdict
    
def load_ensoid(expname,ninoid_name='nino34',datpath=None,standardize=True):
    # Load Enso indices calculated with calc_nino34.py
    if datpath is None:
        datpath = "/home/niu4/gliu8/projects/scrap/nino34/"
    ncname = "%s%s_%s.nc" % (datpath,expname,ninoid_name)
    ds = xr.open_dataset(ncname).load()
    if standardize:
        return ds.sst
    return ds.sst * ds['std'].data.item()

def mcsampler(ts_full,sample_len,mciter,preserve_month=True,scramble_year=False,target_timeseries=None):
    # Given a monthly timeseries [time] and sample length (int), take [mciter] # of random samples.
    # if preserve_month = True, preserve the 12-month sequence as a chunk
    # if scramble_year = True, randomize years that you are selecting from (do not preserve year order)
    # if target_timeseries is not None: also select random samples from list of timeseries (must be same length as ts_full)
    
    # Function Start
    ntime_full        = len(ts_full)
    
    # 1 -- month agnostic (subsample sample length, who cares when)
    if not preserve_month:
        
        print("Month with not be preserved.")
        istarts    = np.arange(ntime_full-sample_len)
        
        sample_ids = []
        samples    = []
        for mc in range(mciter):
            # ts_full[istarts[-1]:(istarts[-1]+sample_len)] Test last possible 
            iistart = np.random.choice(istarts)
            idsel   = np.arange(iistart,iistart+sample_len) 
            msample = ts_full[idsel]
            
            
            sample_ids.append(idsel)
            samples.append(msample) # [iter][sample]
            
        samples = np.array(samples) # [iter x sample]
        # Returns 
            
    elif preserve_month:
        # 2 -- month aware (must select starting points of January + maintain the chunk, preserving the month + year to year autocorrelation)
        if not scramble_year:
            
            # Only start on the year  (to preserve month sequence)
            istarts    = np.arange(0,ntime_full-sample_len,12)
            
            # -------------------- Same as Above
            sample_ids = []
            samples    = []
            for mc in range(mciter):
                # ts_full[istarts[-1]:(istarts[-1]+sample_len)] Test last possible 
                iistart = np.random.choice(istarts)
                idsel   = np.arange(iistart,iistart+sample_len) 
                msample = ts_full[idsel]
                
                sample_ids.append(idsel)
                samples.append(msample) # [var][iter][sample]
            samples = np.array(samples) # [var x iter x sample]
            # -------------------- 
            
        # 3 -- month aware, year scramble (randomly select the year of each month, but preserve each month)
        elif scramble_year: # Scrample Year and Month
            
            # Reshape to the year and month
            nyr_full        = int(ntime_full/12)
            ts_yrmon        = ts_full.reshape(nyr_full,12)
            ids_ori         = np.arange(ntime_full)
            ids_ori_yrmon   = ids_ori.reshape(ts_yrmon.shape)
            
            nyr_sample      = int(sample_len/12)
            sample_ids      = []
            samples         = []
            for mc in range(mciter): # For each loop
                
                # Get start years
                startyears = np.random.choice(np.arange(nyr_full),nyr_sample)
                # Select random years equal to the sample length and combine
                idsel      = ids_ori_yrmon[startyears,:].flatten() 
                # ------
                msample    = ts_full[idsel]
                sample_ids.append(idsel)
                samples.append(msample) # [var][iter][sample]
            samples = np.array(samples) # [var x iter x sample]
            # -----
    
    outdict = dict(sample_ids = sample_ids, samples=samples)    
    if target_timeseries is not None:
        
        sampled_timeseries = []
        for ts in target_timeseries:
            if len(ts) != len(ts_full):
                print("Warning... timeseries do not have the same length")
            randsamp = [ts[sample_ids[mc]] for mc in range(mciter)]
            randsamp = np.array(randsamp)
            sampled_timeseries.append(randsamp) # [var][iter x time]
    outdict['other_sampled_timeseries'] = sampled_timeseries
    
    return outdict

def preprocess_enso(ds):
    # Remove Mean Seasonal Cycle and the Quadratic Trend
    dsds   = proc.xrdeseason(ds)
    
    # Make Function
    def detrend_quadratic(ds):
        x = np.arange(len(ds))
        y = ds.data
        if np.any(np.isnan(y)):
            return np.ones(y.shape)*np.nan
        ydetrended,model=proc.detrend_poly(x,y,2)
        return ydetrended
        
    st = time.time()
    dsanom = xr.apply_ufunc(
        detrend_quadratic,  # Pass the function
        dsds,  # The inputs in order that is expected
        # Which dimensions to operate over for each argument...
        input_core_dims=[['time'],],
        output_core_dims=[['time'],],  # Output Dimension
        vectorize=True,  # True to loop over non-core dims
        )
    print("Detrended in %.2fs" % (time.time()-st))
    return dsanom

def remove_duplicate_times(ds,verbose=True):
    # From : https://stackoverflow.com/questions/51058379/drop-duplicate-times-in-xarray
    _, index = np.unique(ds['time'], return_index=True)
    print("Found %i duplicate times. Taking first entry." % (len(ds.time) - len(index)))
    return ds.isel(time=index)
 
def standardize_names(ds):
    
    ds = swap_rename(ds,'time_counter','time')
    ds = swap_rename(ds,"TIME_COUNTER",'time')
    ds = swap_rename(ds,"LON","lon")
    ds = swap_rename(ds,"LAT","lat")
    ds = swap_rename (ds,"LAT232_409","lat")
    
    # Other preprocessing
    # drop LON_bnds, TIME_COUNTER_bnds
    dropvars = ["LON_bnds","TIME_COUNTER_bnds"]
    for dropvar in dropvars:
        if dropvar in ds:
            ds = ds.drop_vars(dropvar)
    return ds


def stack_events(target_var,eventids,ibefore,iafter):
    # Stack events between -ibefore and +iafter months
    
    nevents        = len(eventids)
    plotlags       = np.hstack([np.flip((np.arange(0,ibefore+1) * -1)),np.arange(1,iafter+1,1)])
    stacked_events = np.zeros((nevents,len(plotlags))) * np.nan
    ntime          = len(target_var)
    for ie in range(nevents):
        
        ievent    = eventids[ie]
        istart    = ievent-ibefore
        iend      = ievent+iafter
        
        if (istart >=0) and (iend < ntime):
            stacked_events[ie,:] = target_var[istart:(iend+1)]
            
        elif iend >= ntime:
            filler = np.zeros( (iend-ntime+1)) * np.nan
            subset = np.hstack([target_var[istart:],filler])
            stacked_events[ie,:] = subset
        elif istart < 0: # Note havent tested this
            filler  = np.zeros(np.abs(istart)) * np.nan
            subset  = np.hstack([filler,target_var[:(iend+1)],])
            stacked_events[ie,:] = subset
    return stacked_events

def stack_events_2d(invar,eventids,ibefore,iafter,times_da=None):
    
    invar = invar.transpose('time','lat','lon')
    ntime,nlat,nlon = invar.shape
    nevents         = len(eventids)
    leadlags        = np.arange(-ibefore,iafter+1)
    nlags           = len(leadlags)
    temp_var        = np.zeros((nevents,nlags,nlat,nlon)) * np.nan
    
    if times_da is not None:
        event_times = np.zeros((nevents,nlags)) * np.nan
    
    for ie in range(nevents):
        
        ievent    = eventids[ie]
        istart    = ievent-ibefore
        iend      = ievent+iafter
        
        # Corrections for end-case indices
        if istart < 0:
            print("istart is at %s" % istart)
            insert_start = np.abs(istart)#(lags - np.abs(istart)).item()
            istart       = 0
        else:
            insert_start = 0
        
        if iend > ntime:
            print("iend is at %s" % (iend))
            insert_end = nlags + (ntime-ievent) #(lags*2+1) - (iend - ntime)
            iend       = ntime
        else:
            insert_end = nlags#lags*2+1
            
        
        indices_in = np.arange(insert_start,insert_end)
        temp_var[ie,indices_in,:,:] = invar[istart:(iend+1),:,:]
        if times_da is not None:
            event_times[ie,indices_in]  = times_da[istart:(iend+1)]

    
    coords             = dict(eventid=eventids,lag=leadlags,lat=invar.lat,lon=invar.lon)
    invar_subset       = xr.DataArray(temp_var,coords=coords,dims=coords)
    if times_da is not None:
        coords_time        = dict(eventid=eventids,lag=leadlags)
        event_times_subset = xr.DataArray(event_times.astype('datetime64[ns]'),coords=coords_time,dims=coords_time)
        return invar_subset,event_times_subset   
    return invar_subset


def swap_rename(ds,chkvar,newvar):
    if chkvar in list(ds.coords):
        print("Renaming [%s] to [%s]" % (chkvar,newvar))
        ds = ds.rename({chkvar:newvar})
    return ds

def varcheck(ds,vname,expname):
    if np.any(ds > 273) and vname == "sst": # Convert to Celsius
        print("Converting from Kelvin to Celsius for %s" % expname)
        ds = ds - 273.15
        
    if vname in ['str','ssr','strc','ssrc','ttr','tsr','ttrc','tsrc','sshf','slhf']: # Accumulation over 3h
        # Conversion for STR and SSR considering 3h Accumulation
        if "TCo319" in expname:
            print("Correction for accumulation over 6 hours for %s" % expname)
            accumulation_hr = 6
        else:
            print("Correction for accumulation over 3 hours for %s" % expname )
            accumulation_hr = 3
        conversion  = 1/(3600 * accumulation_hr)  # 3 h accumulation time...? #1/(24*30*3600)
        # https://forum.ecmwf.int/t/surface-radiation-parameters-joule-m-2-to-watt-m-2/1588
        ds          = ds * conversion
        
    if vname in ["cp","lsp"]: # Convert from [meters/accumulation period] to [mm/day]
        if "TCo319" in expname:
            print("Correction for accumulation over 6 hours for %s" % expname)
            accumulation_hr = 6
        else:
            print("Correction for accumulation over 3 hours for %s" % expname )
            accumulation_hr = 3
        conversion = (24/accumulation_hr) * 1000
        ds         = ds * conversion
    return ds



