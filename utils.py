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
awi_mean_loader             : (l) load mean/monvar/scycle calculations from calc_mean_patterns_TP 
band_avg_spectra            : (g) Take Band-average around ENSO and Combination Mode Frequencies
calc_enso_index             : (g) Calculate ENSO Index
calc_grad_centered          : (g) Calculate centered-difference for spatial gradients
calc_lag_regression_1d      : (g) Compute lead lag regression for 1d timeseries
calc_leadlagreg_pointwise   : (g) Compute pointwise lead lag regression using ufunc
calc_leadlag_regression_2d  : (g) Compute lead lag regression for 2D timseries (modeled after enso_lag_regression)
calc_spmm                   : (g) Calculate the South Pacific Meridional Mode
center_events_ninodict      : (c) Center Identified Events around a particular month
convolve_kernel_ccf         : (c) Convolve radiative kernel with ccf variable 
combine_events              : (g) Given identified events, combine similar events and get other traits (duration, etc)
fit_ctone_enso              : (g) Fit Combination Tone Model
generate_periods            : (c) Generate Periods for a sliding window of a selected length
get_early_late_periods      : (c) Get N earliest/latest periods for a DataArray
get_moving_segments         : (g) Subset timeseries into segments with a sliding window
get_rawpath_awi             : (A) Get rawpath for AWI output on niu
get_center_time             : (g) Get Center Time given min and max of time range
identify_stack_events       : (g) (WIP) Need to add this from lead_lag_composite_SEP_Box.ipynb)
init_tp_map                 : (v) initialize tropical Pacific plot 
init_global_map             : (v) Initialize a global map
init_sep_map                : (v) Initializw SEP Map
init_sliding_spectra_plot   : (v) Initialize sliding spectra plot...
load_ccf_kernel             : (l) load Kernels for each CCF (based on calc_ccf_radiation_byexp)
load_ensoid                 : (l) load enso indices calculated by calc_nino34.py
load_enso_eof               : (l) load ENSO indices from rotated EOF analysis (calc_EOF_enso<_sliding>.py)
load_land_mask_awi          : (l) load land mask from AWI-CM3, where land points are np.nan
load_regrid                 : (l) load regridded 1x1 output
load_scott2020_kernels      : (l) load radiative kernels from Scott et al. 2020 Paper
make_ninotime               : (A) Get Center Time from time range
mcsample                    : (g) Monte Carlo Sampler to repeat function
mlr                         : (g) single point multiple linear regression using Scipy
mlr_ccfs                    : (c) Multiple linear regrssion for cloud controlling factor analysis
preprocess_byperiod         : (g) Remove mean seasonal cycle and linear detrend by period
preprocess_enso             : (c) detrend (quadratic) and deseasonalize for ENSO calculations
remove_duplicate_times      : (g) Remove duplicate times from a DataArray
shift_time_monthstart       : (g) Shift times from center to start of month
sliding_spectra             : (g) Preprocess by period then compute sliding power spectra
swap_rename                 : (g) check if variable exists and rename if so
stack_events                : (g) For 1-d timeseries, stack events along specified leads/lags
stack_events_2d             : (g) Stack Events, but applied to 2D case with time x lat x lon....
standardize_names           : (A) uses swap_rename to replace variable and dimension names in AWI_CM3 output 
varcheck                    : (A) checks and converts variables for AWI-CM3

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
import glob
import pandas as pd
from sklearn.linear_model import LinearRegression
import sklearn
import scipy as sp

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

def band_avg_spectra(spectra_in,cutoff_periods_month=None,
                     debug=True,return_all=False,return_ds=True,band_sum=False):
    # Take Band-average around ENSO and Combination Mode Frequencies
    # spectra_in is a Dataset with spectra and freq dimensions (in 1/sec)
    # cutoff_periods is a list of bounds [[lower,upper]] in Months
    # Note, Debug only works with ENSO Cutoffs
    # Take from: sliding_spectra_ccf_region.ipynb
    # For Debugging band_sum, see bottom of: check_spectra_estimate_ctone.ipynb
    
    dtmon           = 3600*24*30
    bpnames         = ["Decadal to ENSO","ENSO Max Freq.","Difference Tone","Annual Band","Sum Tone","Sub-Annual"]
    bbcolors        = ["gray","red","blue","k","yellow","cyan"]
    
    
    if cutoff_periods_month == None:
        cutoff_periods_month         = [[5.5*12,10*12,],[2*12,5.5*12],[15,2*12],[10,15],[8,10],[0,8]] # Period in Months
    nthres = len(cutoff_periods_month)
    
    # Convert cutoffs to Seconds
    cutoffs_sec     = np.array(cutoff_periods_month) * dtmon
    # Convert cutoffs to Frequencies (1/sec)
    cutoffs_freq    = 1/cutoffs_sec
    
    # Select Frequencies
    bbn             = len(cutoff_periods_month)
    spec_bybb       = []
    specmean_bybb   = []
    for bb in range(nthres):
        freqsel = cutoffs_freq[bb]
        specsel = spectra_in.sel(freq=slice(freqsel[1],freqsel[0])).spectra
        spec_bybb.append(specsel)
        if band_sum:
            # Sum the Spectra (trapz)
            specmean_bybb.append(specsel.integrate(coord='freq'))
        else:
            specmean_bybb.append(specsel.mean('freq'))
    
    if debug:
        
        fig,ax,ax2=viz.init_specplot_enso(1,1)
    
        #for pp in range(2):
        specout = spectra_in
        lab     = "Spectra, Tcenter=%s" % specout.tcenter.data.item() 
        ax.loglog(specout.freq*dtmon,specout.spectra/dtmon,c="k",label=lab,lw=2)
        ax.loglog(specout.freq*dtmon,specout.CC.isel(clvl=1)/dtmon,c="k",ls='dotted',lw=1)

        
        ax.set_ylabel("Power")
        viz.add_ctones()
        
        # Plot the selected frequencies
        for bb in range(nthres):
        
            freqplot = cutoffs_freq[bb]
            spec2    = spec_bybb[bb]
            avgval   = specmean_bybb[bb].data.item()
        
            
            ax.loglog(spec2.freq*dtmon,spec2/dtmon,c=bbcolors[bb],ls='dashed',lw=4,alpha=0.8,label=lab)
            
        
            
            specm  = np.ones(len(spec2.freq)) * avgval
            lab    = "%s = %.2e" % (bpnames[bb],avgval)
            ax.plot(spec2.freq*dtmon,specm/dtmon,c='gray',ls='solid',lw=2,alpha=0.8,label="")
            
            

        ax.legend()
    
    if return_ds:
        da_values = xr.concat(specmean_bybb,dim='band')
        da_values['band'] = bpnames
        coords2   = dict(band=bpnames,bounds=['lower','upper'])
        da_freqs  = xr.DataArray(cutoff_periods_month,coords=coords2,dims=coords2,name='cutoffs')
        specmean_bybb = xr.merge([da_values,da_freqs])
    
    if return_all:
        return specmean_bybb,spec_bybb
    return specmean_bybb

def calc_enso_index(sstanom,ninoid_name='nino34',apply_movmean=False):
    """
    Calculate ENSO Index by taking area-weighted average of SST anomalies, with
    optional running mean.

    Inputs:
        ninoid_name (str)           : Name of Index/Box (nino12,nino3,nino34,nino4)
        sstanom     (xr.DataArray)  : Anomalized monthly SSTs with dimensions ('lat','lon','time')
        apply_movmean (bool)        : True to apply 5-month running mean
    Outputs:
        ensoid      (xr.DataArray)  : ENSO Index
    """
    
    # ENSO Index Boxes from UCAR's Climate Data Guide
    bbox_nino12       = [-90+360,-80+360,-10,0]      # Nino 1+2
    bbox_nino3        = [-150+360, -90+360 , -5, 5]  # Nino 3 Box
    bbox_nino34       = [-170+360,-120+360,-5,5]     # Nino3.4 Box
    bbox_nino4        = [ 160    ,-150+360,-5,5]     # Nino 4 Box
    bbox_cti        =   [-180+360,-90+360,-6 ,6]       # Cold Tongue Index
    if ninoid_name == "nino34":
        bbox = bbox_nino34
    elif ninoid_name == 'nino3':
        bbox = bbox_nino3
    elif ninoid_name == "nino4":
        bbox = bbox_nino4
    elif ninoid_name == "nino12":
        bbox = bbox_nino12
    elif ninoid_name == "cti":
        bbox = bbox_cti
    
    # Restrict to box and take area average
    sst_enso_box = proc.sel_region_xr(sstanom,bbox)
    ensoid       = proc.area_avg_cosweight(sst_enso_box)
    
    # Optionally Take Moving Mean
    if apply_movmean:
        movmean = lambda ds,win: np.convolve(ds.data,np.ones(win)/win,mode='same')
        ensoid  = movmean(ensoid,5)  # 5-month running mean
        coords  = {'time':sstanom.time}
        ensoid  = xr.DataArray(ensoid,coords=coords,dims=coords,name=ensoid.name)
    return ensoid

def calc_grad_centered(ds,latname='lat',lonname='lon'): # Copied from structure in calc_ekman_advection_htr
    
    if latname != 'lat' or lonname != 'lon':
        lldict          = {latname : 'lat', lonname : 'lon'}
        lldict_reverse  = {'lat' : latname, 'lon' : lonname}
        rename_flag=True
        ds = ds.rename(lldict)
    else:
        rename_flag=False
    
    dx,dy = proc.calc_dx_dy(ds.lon.values,ds.lat.values,centered=True)
    
    # Convert to DataArray
    daconv   = [dx,dy]
    llcoords = {'lat':ds.lat.values,'lon':ds.lon.values,}
    da_out   = [xr.DataArray(ingrad,coords=llcoords,dims=llcoords) for ingrad in daconv]
    dx,dy = da_out
    
    # Roll and Compute Gradients (centered difference)
    ddx = (ds.roll(lon=-1) - ds.roll(lon=1)) / dx
    ddy = (ds.roll(lat=-1) - ds.roll(lat=1)) / dy
    ddy.loc[dict(lat=ddy.lat.values[-1])] = 0 # Set top latitude to zero (since latitude is not periodic)
    
    if rename_flag:
        ddx,ddy =  [ds.rename(lldict_reverse) for ds in [ddx,ddy]]
    
    return ddx,ddy

def calc_lag_regression_1d(var1_lag,var2_base,lags,correlation=False): # CAn make 2d by mirroring calc_lag_covar_annn
    # Calculate the regression where
    # (+) lags indicate var1 lags  var2 (var 2 leads)
    # (-) lags indicate var1 leads var2 (var 1 leads)
    
    if np.any(np.isnan(var1_lag)) or np.any(np.isnan(var2_base)):
        #print("NaN detected. Returning NaN...")
        return np.nan * np.ones(len(lags*2))
    if np.any(var1_lag == 0.) or np.any(var2_base == 0.):
        return np.nan * np.ones(len(lags*2))
    
    ntime = len(var1_lag)
    betalag = []
    poslags = lags[lags >= 0]
    for l,lag in enumerate(poslags):
        varlag   = var1_lag[lag:]
        varbase  = var2_base[:(ntime-lag)]
        
        # Calculate correlation
        if correlation:
            beta = np.corrcoef(varbase,varlag)[0,1]
        else:
            beta = np.polyfit(varbase,varlag,1)[0]   
        betalag.append(beta.item())
    
    neglags = lags[lags < 0]
    neglags_sort = np.sort(np.abs(neglags)) # Sort from least to greatest #.sort
    betalead = []
    
    for l,lag in enumerate(neglags_sort):
        varlag   = var2_base[lag:] # Now Varlag is the base...
        varbase  = var1_lag[:(ntime-lag)]
        # Calculate correlation
        if correlation:
            beta = np.corrcoef(varlag,varbase)[0,1]
        else:
            beta = np.polyfit(varlag,varbase,1)[0]   
        betalead.append(beta.item())
    
    # Append Together
    return np.concatenate([np.flip(np.array(betalead)),np.array(betalag)])

# Copied from xrfunc in amv (based on leadlagcorr.
def calc_leadlagreg_pointwise(dslag,dsbase,leadlags):
    # Computes lead lag regression between ds1 and ds2 over dimension ['time']
    # over specified lead/lags [lags]. Loops over all other dimensions (lat,lon,depth,etc)
    # Copied over from  calculate_pointwise_sw_feedback.py
    # Might need to ensure there are no NaNs or everything is not Zero... 
    def calc_leadlag_regr(varlag,varbase):
        try:
            llout = calc_lag_regression_1d(varlag,varbase,leadlags)
            return llout
        except:
            return np.ones(len(leadlags)) * np.nan
    
    llreg = xr.apply_ufunc(
        calc_leadlag_regr,
        dslag,
        dsbase,
        input_core_dims=[['time'],['time']],
        output_core_dims=[['lags']],
        vectorize=True,
        )
    #leadlags     = np.concatenate([np.flip(-1*lags)[:-1],lags],) 
    llreg['lags'] = leadlags
    return llreg


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
            maxval     = np.nanmax(np.abs(amplitudes.data)).item() * np.sign(amplitudes.data[idmax])
            event_max.append(maxval)
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


def fit_ctone_enso(anomalies,ensoid,tmax=None,initial_guess=None,debug=True,
                   use_sine=False,linear_enso=False,verbose=False,fixed_params=None,
                   fix_frequency_offset=False):
    # Fit Idealized Combination Tone Model to Anomalie Timeserise [anomalies]
    # Main difference from [proc.fit_sinfunc] above is that we multiple sine function
    # By ENSO Timeseries [ensoid] and allow for flexibility in the period
    # Inputs:
    # [anomalies] (xr.DataArray)    : 1-d anomaly timeseries to fit
    # [ensoid] (xr.DataArray)       : 1-d ENSO Timeseries
    # [tmax]  (int)                 : max index of timeseries to subset...
    # [Initial Guess] (array)       : Guess of initial parameter values. None to use default guess.
    # [debug] (bool)                : True to Make Timeseries Plot
    # [linear_enso] (bool)          : True to add linear ENSO term (beta) Nino3.4
    # [use_sine] (bool)             : True to assume sin(a+b), False to use cos(a+b) (yields same answer...)
    # [fix_frequency_offset] (bool) : Use fixed annual cycle freq (2*pi/12) and fixed offset (0)
    #
    # Outputs: Dict Containing:
    # [Parameter Estimates] (int), (amplitude,frequency,phase,offset)
    # [prediction] (xr.DataArray) : Modeled 1-D anomaly timeseries
    # [fitfunc] (function)        : Function with parameters input
    # [period] (numeric)          : 1/frequency, multiply by 2*np.pi to convert
    # [corr] (numeric)            : Correlation between modeled + actual y
    # [r2] (numeric)              : r2between modeled + actual y
    # 
    # For Debugging Script, see `simple_mode_ctone.ipynb`
    
    # Get Length of Time
    if tmax is not None:
        y      = anomalies[:tmax]
        ensoin = ensoid[:tmax]
    else:
        y      = anomalies.copy()
        ensoin = ensoid.copy()
    ntime = len(y)
    x     = np.arange(ntime)
    
    # Make the function (based on angle sum identities)
    if use_sine: # Sin(a+b) = sin(a)cos(b) + cos(a)sin(b)
        if linear_enso:
            def sine_func_ctone(t,amplitude,frequency,phase,offset,beta):
                return beta * ensoin + ensoin * amplitude * ( (np.cos(phase) * np.sin(frequency*t)) + (np.sin(phase) * np.cos(frequency*t)) ) + offset
        else:
            def sine_func_ctone(t,amplitude,frequency,phase,offset):
                return ensoin * amplitude * ( (np.cos(phase) * np.sin(frequency*t)) + (np.sin(phase) * np.cos(frequency*t)) ) + offset
    else:        # Cos(a+b) = cos(a)cos(b) - sin(a)sin(b)
        if linear_enso:
            
            def sine_func_ctone(t,amplitude,frequency,phase,offset,beta): # Misnomer should be cos
                return beta*ensoin + ensoin * amplitude * ( (np.cos(phase) * np.cos(frequency*t)) - (np.sin(phase) * np.sin(frequency*t)) ) + offset
        else:
                        
            def sine_func_ctone(t,amplitude,frequency,phase,offset): # Misnomer should be cos
                return ensoin * amplitude * ( (np.cos(phase) * np.cos(frequency*t)) - (np.sin(phase) * np.sin(frequency*t)) )
    
    # Make initial Guess
    if initial_guess is None:
        initial_guess = [2*np.nanstd(y), # Amplitude ~ 2*Standard Deviation
                         np.pi*2/12    , # Frequency ~ Annual Cycle 
                         0             , # Assume no phase shift
                         np.nanmean(y) , # Offset = Mean of timeseries, ~0
                        ]
        if linear_enso: # First Guess is the slope between the two
            coeffs     = np.polyfit(ensoin.data,anomalies.data,1)
            guess_beta = coeffs[0]
            initial_guess.append(guess_beta) # Relationship to ENSO

    # Make Adjustments for Fixed Frequency
    if fix_frequency_offset:
        
        # Set the Period and Offset
        fix_frequency = 2*np.pi / 12 # Seasonal Cycle
        fix_offset = 0            # Assume anomalies centered on 0
        frequency = fix_frequency
        offset    =   fix_offset
        
        # Update the Function 
        if use_sine:
            if linear_enso:
                def sine_func_ctone(t,amplitude,phase,beta):
                    return beta * ensoin + ensoin * amplitude * ( (np.cos(phase) * np.sin(fix_frequency*t)) + (np.sin(phase) * np.cos(dix_frequency*t)) ) + fix_offset
            else:
                def sine_func_ctone(t,amplitude,phase):
                    return ensoin * amplitude * ( (np.cos(phase) * np.sin(fix_frequency*t)) + (np.sin(phase) * np.cos(dix_frequency*t)) ) + fix_offset
        else:
            if linear_enso:
                def sine_func_ctone(t,amplitude,phase,beta): # Misnomer should be cos
                    return beta*ensoin + ensoin * amplitude * ( (np.cos(phase) * np.cos(fix_frequency*t)) - (np.sin(phase) * np.sin(fix_frequency*t)) ) + fix_offset
            else:
                def sine_func_ctone(t,amplitude,phase): # Misnomer should be cos
                    return ensoin * amplitude * ( (np.cos(phase) * np.cos(fix_frequency*t)) - (np.sin(phase) * np.sin(fix_frequency*t)) ) + fix_offset
        
        # Update Initial Guess
        initial_guess.remove(fix_frequency)
        initial_guess.remove(np.nanmean(y))
    
    
    # Make the Fit
    params, covariance = sp.optimize.curve_fit(sine_func_ctone, x, y, p0=initial_guess)

    # Extract Parameters
    if fix_frequency_offset:
        if linear_enso:
            amplitude,phase,beta = params
            param_string         = f"Fitted parameters: Amplitude={amplitude}, Phase={phase}, Beta={beta}"
        else:
            amplitude,phase      = params
            param_string         = f"Fitted parameters: Amplitude={amplitude}, Phase={phase}"
    else:
        if linear_enso:
            amplitude, frequency, phase, offset, beta = params
            #param_string                              = f"Fitted parameters: Amplitude={amplitude}, Frequency={frequency}, Phase={phase}, Offset={offset}, Beta={beta}"
            #param_string_title                        = f"Fitted parameters: Amplitude ($A_A$)={amplitude:.2f}, Period ($\omega_A$, months)={(2*np.pi)/frequency:.2f}, Phase ($\phi_A$, months)={(2*np.pi)/phase:.2f}, Offset ($C_A$)={offset:.2f}, Beta={beta:2f}"
        else:
            amplitude, frequency, phase, offset       = params
            #param_string                              = f"Fitted parameters: Amplitude={amplitude}, Frequency={frequency}, Phase={phase}, Offset={offset}"
            #param_string_title                        = f"Fitted parameters: Amplitude ($A_A$)={amplitude:.2f}, Period ($\omega_A$, months)={(2*np.pi)/frequency:.2f}, Phase ($\phi_A$, months)={(2*np.pi)/phase:.2f}, Offset ($C_A$)={offset:.2f}"
    if linear_enso:
        param_string_title                        = f"Fitted parameters: Amplitude ($A_A$)={amplitude:.2f}, Period ($\omega_A$, months)={(2*np.pi)/frequency:.2f}, Phase ($\phi_A$, months)={(2*np.pi)/phase:.2f}, Offset ($C_A$)={offset:.2f}, Beta={beta:2f}"
    else:
        param_string_title                        = f"Fitted parameters: Amplitude ($A_A$)={amplitude:.2f}, Period ($\omega_A$, months)={(2*np.pi)/frequency:.2f}, Phase ($\phi_A$, months)={(2*np.pi)/phase:.2f}, Offset ($C_A$)={offset:.2f}"
    if verbose:
        print(param_string)
    
    # Create Model
    if fix_frequency_offset:
        if linear_enso:
            model = lambda t: sine_func_ctone(t,amplitude,phase,beta)
        else:
            model = lambda t: sine_func_ctone(t,amplitude,phase)
    else:
        if linear_enso:
            model = lambda t: sine_func_ctone(t,amplitude,frequency,phase,offset,beta)
        else:
            model = lambda t: sine_func_ctone(t,amplitude,frequency,phase,offset)
    y_fit  = model(x)
    r2     = np.corrcoef(y,y_fit)[0,1]**2
    
    # Plot Output
    if debug:
        fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(12.5,4))
        plotx  = y.time
        #l1 = ax.scatter(plotx,y,alpha=0.8,c='k',marker="o",s=5,label="Target Anomaly Timeseries")
        l1 = ax.plot(plotx,y,alpha=0.8,c='k',marker="o",markersize=1.5,linewidth=0.55,label="Target Anomaly Timeseries",ls='dashed')
        #ax.plot(plotx,y,alpha=0.8,c="k",lw=0.25)
        l2 = ax.plot(plotx,y_fit,c='red',lw=2.5,label="Combination Tone Fit, $r^2$=%.2f%%" % (r2*100),alpha=0.6)
    
        ax2 = ax.twinx()
        l3 = ax2.plot(plotx,ensoin,c='gray',lw=0.5,label="ENSO Index")
        ax2 = viz.change_axcol('right',c='gray',ax=ax2)
        ax2.set_ylabel(r"ENSO Index ($\degree C)$")
        #viz.add_axlines(xonly=True)
    
        #lns  = [l1,]+l2+l3 # For Scatter Case
        lns  = l1+l2+l3
        labs = [l.get_label() for l in lns]
        ax.legend(lns,labs,ncol=3)
        ax.set_xlabel("Time")
        ax.set_ylabel("Anomaly")
        ax.set_title(param_string_title)
        
    # Make Output Dictionary
    # Prepare Output
    outdict = dict(
        amplitude   = amplitude,
        frequency  = frequency,
        phase = phase,
        offset= offset,
        prediction = y_fit,
        fitfunc=model,
        period= 1/frequency,
        corr=np.corrcoef(y,y_fit),
        r2=r2,
        )
    if linear_enso:
        outdict['beta'] = beta
    return outdict

def fit_ctone_pointwise(ds,ensoid,return_ds=True,
                        linear_enso=False,fix_frequency_offset=False):
    
    # Pointwise application of fit_sinfunc
    def unpack_sinfit(target,ensoid,debug=False):
        t      = np.arange(len(target))
        try:
            if linear_enso:
                fitout = fit_ctone_enso(target,ensoid,debug=False,fix_frequency_offset=fix_frequency_offset,
                                        linear_enso=True)
                amp,freq,phase,offset,ypred,yfunc,period,corr,r2,beta = list(fitout.values())
            else:
            
                fitout = fit_ctone_enso(target,ensoid,debug=False,fix_frequency_offset=fix_frequency_offset)
                amp,freq,phase,offset,ypred,yfunc,period,corr,r2 = list(fitout.values())
        except:
            amp  = np.nan 
            freq = np.nan
            phase = np.nan
            offset = np.nan
            ypred  = t * np.nan
            r2    = np.nan
            beta = np.nan
        
        if linear_enso:
            return amp,freq,phase,offset,ypred,r2,beta
        return amp,freq,phase,offset,ypred,r2 # Just retain some of the output
    
    stxr = time.time()
    if linear_enso:
        output_core_dims = [[],[],[],[],['time'],[],[]]
    else:
        output_core_dims = [[],[],[],[],['time'],[]]
    dsout = xr.apply_ufunc(
        unpack_sinfit,
        ds,
        ensoid,
        input_core_dims=[['time'],['time']],
        output_core_dims=output_core_dims,
        vectorize=True,
        )
    print("Completed fit in %.2fs" % (time.time()-stxr))
    
    # Place in DataSet
    if return_ds:
        if linear_enso:
            amp,freq,phase,offset,ypred,r2,beta = dsout # dsout
        else:
            amp,freq,phase,offset,ypred,r2 = dsout # dsout
        dsout = xr.merge([amp.rename('amplitude'),
                          freq.rename('frequency'),
                         phase.rename('phase'),
                         offset.rename('offset'),
                         ypred.rename('ypred'),
                         r2.rename('r2')])
        if linear_enso:
            dsout['beta'] = beta.rename('beta')
    return dsout


def calc_leadlag_regression_2d(ensoid,dsvar,leadlags,sep_mon=False):
    # Based on routine in enso_lag_regression and global_mean_nino_regressions
    # Compute lead/lag regression on timeseries [ensoid] and 3D lagged variable [dsvar]
    
    adddim=False
    if len(dsvar.shape) == 1:
        dummylatlon = dict(lat=[1,],lon=[1,])
        dsvar = dsvar.expand_dims(dim=dummylatlon,axis=(0,1))
        adddim=True
    
    if dsvar.name is None:
        dsvar = dsvar.rename("regression_coefficient")
    
    # Check to make sure the time matches
    dsvar,ensoid    = proc.match_time_month(dsvar,ensoid)
    
    # Get Dimension Lengths
    dsvar           = dsvar.transpose('lon','lat','time')
    nlon,nlat,ntime = dsvar.shape
    
    if sep_mon is False: # Do for all months
        # Do the Leads (variable leads)
        leads       = np.abs(leadlags[leadlags <=0])
        nleads      = len(leads)
        beta_leads  = np.zeros((nlon,nlat,nleads)) * np.nan
        sig_leads   = beta_leads.copy()
        for ll in range(nleads):
            lag                = leads[ll] 
            ints               = ensoid.data[lag:]
            invar              = dsvar.data[:,:,:(ntime-lag)]
            rout               = proc.regress_ttest(invar,ints,verbose=False)
            beta_leads[:,:,ll] = rout['regression_coeff']
            sig_leads[:,:,ll]  = rout['sigmask']
            
        # Do the lags
        lags        = leadlags[leadlags > 0]
        nlags       = len(lags)
        beta_lags   = np.zeros((nlon,nlat,nlags)) * np.nan
        sig_lags    = beta_lags.copy()
        for ll in range(nlags):
            lag   = lags[ll] 
            ints  = ensoid.data[:(ntime-lag)]
            invar = dsvar.data[:,:,lag:]
            rout  = proc.regress_ttest(invar,ints,verbose=False)
            beta_lags[:,:,ll] = rout['regression_coeff']
            sig_lags[:,:,ll]  = rout['sigmask']
        
        # Concatenate
        betas = np.concatenate([beta_leads,beta_lags],axis=2)
        sigs  = np.concatenate([sig_leads,sig_lags],axis=2)
        
        # Replace into DataArray
        coords   = dict(lon=dsvar.lon,lat=dsvar.lat,lag=leadlags)
        da_betas = xr.DataArray(betas,coords=coords,dims=coords,name=dsvar.name)
        da_sigs  = xr.DataArray(sigs,coords=coords,dims=coords,name="sig")
        da_out   = [ds.transpose('lag','lat','lon') for ds in [da_betas,da_sigs]]
        
    else: # Calculate Separately by Month
        print("Warning: This is currently not working properly...")
        # Reshape the variables
        nyr          = int(ntime/12)
        ints_yrmon   = ensoid.data.reshape(nyr,12)
        invar_yrmon  = dsvar.data.reshape(nlon,nlat,nyr,12)
        
        # Calculate Leads (a bit silly to write it this way I know... should prob make a fuction)
        leads       = np.abs(leadlags[leadlags <=0])
        nleads      = len(leads)
        beta_leads  = np.zeros((12,nlon,nlat,nleads)) * np.nan
        sig_leads   = beta_leads.copy()
        for ll in range(nleads):
            lag                = leads[ll]  # This approach is incorrect I think, need to carefully apply the lag
            if (lag >= nyr):
                print("Cannot perform calculation since the lag (%i) exceeds the # of years %02i. Skipping" % (lag,nyr))
                continue
            for im in range(12):
                ints                    = ints_yrmon[lag:,im]
                invar                   = invar_yrmon[:,:,:(nyr-lag),im] #dsvar.data[:,:,:(ntime-lag)]
                rout                    = proc.regress_ttest(invar,ints,verbose=False)
                beta_leads[im,:,:,ll]   = rout['regression_coeff']
                sig_leads[im,:,:,ll]    = rout['sigmask']
                
        # Calculate Lags
        lags        = leadlags[leadlags > 0]
        nlags       = len(lags)
        beta_lags   = np.zeros((12,nlon,nlat,nlags)) * np.nan
        sig_lags    = beta_lags.copy()
        for ll in range(nlags):
            lag   = lags[ll] 
            if (lag >= nyr):
                print("Cannot perform calculation since the lag (%i) exceeds the # of years %02i. Skipping" % (lag,nyr))
                continue
            for im in range(12):
                ints  = ints_yrmon[:(nyr-lag),im]
                invar = invar_yrmon[:,:,lag:,im]
                rout  = proc.regress_ttest(invar,ints,verbose=False)
                beta_lags[im,:,:,ll] = rout['regression_coeff']
                sig_lags[im,:,:,ll]  = rout['sigmask']
              
        # Concatenate
        betas = np.concatenate([beta_leads,beta_lags],axis=3)
        sigs  = np.concatenate([sig_leads,sig_lags],axis=3)
        
        # Replace into DataArray
        coords   = dict(mon=np.arange(1,13,1),lon=dsvar.lon,lat=dsvar.lat,lag=leadlags)
        da_betas = xr.DataArray(betas,coords=coords,dims=coords,name=dsvar.name)
        da_sigs  = xr.DataArray(sigs,coords=coords,dims=coords,name="sig")
        da_out   = [ds.transpose('lag','mon','lat','lon') for ds in [da_betas,da_sigs]]
    
    # Merge Variables
    if adddim:
        da_out = [ds.squeeze() for ds in da_out]
    ds_out   = xr.merge(da_out)
    
    return ds_out

def calc_spmm(ds_sst,ds_u10,ds_v10,debug=False,regress_wind=False,verbose=True,standardize=True,ensoid=None,output_ds=True):
    """
    
    Calculate the South Pacific Meridional Mode (SPMM) given 3 global DataArrays of SST, 10m Zonal Wind, and 10m Meridional Winds.
    Input DataArrays must have matching [lat x lon x time] dimensions.
    ENSO is removed via regression (EOF-based indices are automatically calculated if [ensoid] is not supplied)
    Note that [verbose] triggers timing messages, while [debug] triggers verbose messages within functions as well
    Uses the [xmca] package. 
    
    Inputs
    ------
        ds_sst         (xr.DataArray) - Global SST with [lat x lon x time] dimensions
        ds_u10         (xr.DataArray) - Global 10m zonal winds with [lat x lon x time] dimensions
        ds_v10         (xr.DataArray) - Global 10m meridional winds with [lat x lon x time] dimensions
        debug          (Bool)         - True to make debugging plots and print verbose messages in internal functions. Default is False.
        regress_wind   (Bool)         - True to also perform regression on Wind PCs for EOF1. Default is False.
        verbose        (Bool)         - True to print steps and timing. Default is True.
        standardize    (Bool)         - True to standardize PCs prior to regression. Default is True
        ensoid         (xr.DataArray) - ENSO Index to remove via linear regression. Must be 1-D and have [time] dimension. Calculates EOF-based ENSO if not supplied.
        output_ds      (Bool)         - True to output Dataset. False outputs dictionary, but more variables are included in regression such as r2, residual, etc.
        
    Output : outdict containing...
    ------
        If [output_ds] is True, following is a merged Dataset. Otherwise, it is a dictionary with the corresponding keys.
        pc_sst   (xr.DataArray)  - SST Principle Components (Time x Mode)
        pc_wind  (xr.DataArray)  - Wind Principle Components (Time x Mode)
        eof_sst  (xr.DataArray)  - Regression of SST PC to variables. (Variable x Lat x Lon). 
                                    If output_ds=False, eof_sst is a nested list of [variable][ds] where ds contains full pointwise_polyfit output...
        eof_wind (xr.DataArray)  - (IF regress_wind is True) Regression of Wind PC to variables. (Variable x Lat x Lon)
                                    If output_ds=False, eof_sst is a nested list of [variable][ds] where ds contains full pointwise_polyfit output...
        varexp   (xr.DataArray)  - Variance Explained (Mode)

    WIP:
        Adding 3-month sliding window prior to EOF Analysis
        Option to reshape and directly output EOF from MCA SVD
        Monthly Implementation
        More Flexible ENSO Index Removal (or option to bypass)
        Support for large files (or providing non-global inputs)
        Use non xMCA package (due to lack of updates)
        See [niu:notebooks/ccfs/mca_test.ipynb] for debugging scrap...
    
    Last Updated 2026.06.17, by Glenn Liu
    
    """
    # 
    
    st_all            = time.time()
    # Bounding Boxs for Calculations
    bbox_tropical     = [140,-90+360,-20,20]         # EOF Bounding Box
    bbox_nino3        = [-150+360, -90+360 , -5, 5]  # Nino 3 Box for Checking Sign

    # Step 1. Preprocess (Deseason and Linear Detrend) ------------------------------------
    st      = time.time()
    dsraws  = [ds_sst,ds_u10,ds_v10]
    vnames  = ["sst","u10","v10"]
    dsraws  = [ds.transpose('time','lat','lon') for ds in dsraws]
    dsanoms = [proc.xrdeseason(ds,verbose=debug) for ds in dsraws]
    dsanoms = [proc.xrdetrend(ds,verbose=debug) for ds in dsanoms]
    if verbose:
        proc.printtime(st,"(1/7) Deseason/Detrended")
    
    # Step 2. Calculate ENSO Index using EOF Analysis ------------------------------------
    st                = time.time()
    if ensoid is None: # Calculate ENSO using EOF analysis
        N_mode            = 2
        sst_tropical      = proc.sel_region_xr(dsanoms[0],bbox_tropical)
        eofout            = proc.eof_time_ds(sst_tropical,N_mode,check_sign=[bbox_nino3,],flip_if_negative=True,verbose=debug)
        ninoid            = eofout.pcs.isel(mode=0) # Get ENSO Index (First Mode)
    else: # Check if ENSO Index matches the length of the anoms
        ntime_enso        = len(ensoid)
        for ii in range(3):
            ntime_var = len(dsanoms[ii].time)
            if ntime_enso != ntime_var:
                if verbose:
                    print("Warning, Time does not match between ENSO Index and %s" % (vnames[vv]))
                ensoid,var_adjust=proc.match_time_month(ensoid,dsanoms[ii])
                dsanoms[ii] = var_adjust
        ninoid = ensoid 
    if verbose:
        proc.printtime(st,"(2/7) Performed ENSO EOF")
    
    # Step 3. Remove ENSO influence via linear regression --------------------------------
    st                = time.time()
    dtouts            = [proc.detrend_by_regression(ds,ninoid,verbose=debug) for ds in dsanoms]
    if verbose:
        proc.printtime(st,"(3/7) Performed ENSO EOF")
    
    # Step 4. Restrict to SPMM Region and do further preprocessing -----------------------
    st           = time.time()
    bbox_spmm    = [-180+360,-70+360,-35,-20]
    dsanoms_reg  = [proc.sel_region_xr(dtouts[vv][vnames[vv]],bbox_spmm) for vv in range(3)]
    
    #         (add option to apply 3-month running mean here...)
    
    #         Convert to array
    dsanoms_reg_arr     = [ds.data for ds in dsanoms_reg]
    
    #         Apply Sqrt Cosine Latitude Weight 
    lon,lat             = dsanoms_reg[0].lon,dsanoms_reg[0].lat
    xx,yy               = np.meshgrid(lon,lat)
    dsanoms_weighted    = [ds * (np.cos(np.deg2rad(yy))**(0.5))[None,:,:] for ds in dsanoms_reg_arr]
    
    #         Reshape and merge surface wind arrays along space dimension
    ntime,nlat,nlon     = dsanoms_reg[0].shape
    dsanoms_reshaped    = [ds.reshape(ntime,nlat*nlon) for ds in dsanoms_weighted]
    ds_in               = [dsanoms_reshaped[0],np.concatenate([dsanoms_reshaped[1],dsanoms_reshaped[2]],axis=1)]
    if verbose:
        proc.printtime(st,"(4/7) Reshaped + Formatted for MCA")
    
    # Step 5. Perform MCA using xMCA -----------------------------------------------------
    st     = time.time()
    mca    = MCA(ds_in[0],ds_in[1])
    mca.solve()
    expvar = mca.explained_variance()
    eofs   = mca.eofs()
    pcs    = mca.pcs()
    
    #         Post-process the output (Array --> DataArray)
    sstpcs        = np.array(pcs['left'])  # [Time x Mode]
    windpcs       = np.array(pcs['right']) # [Time x Mode]
    coords_pc     = dict(time=dsanoms_reg[0].time,mode = np.arange(ntime))
    da_sstpcs     = xr.DataArray(sstpcs,coords=coords_pc,dims=coords_pc)
    da_windpcs    = xr.DataArray(windpcs,coords=coords_pc,dims=coords_pc)
    coords_exp    = dict(mode=np.arange(ntime))
    da_varexp     = xr.DataArray(expvar,coords=coords_exp,dims=coords_exp)
    if verbose:
        proc.printtime(st,"(5/7) Performed MCA + Postprocesed")
    
    # Step 6. Regress Standardized Values back to variables -----------------------------
    st              = time.time()
    spmm_index      = da_sstpcs.isel(mode=0)
    if standardize:
        spmm_index  = spmm_index / spmm_index.std('time')
    spmm_regression = [proc.pointwise_polyfit(spmm_index,ds,1) for ds in dsanoms]
    if regress_wind:
        spmm_index_wind = da_windpcs.isel(mode=0)
        if standardize:
            spmm_index_wind  = spmm_index_wind / spmm_index_wind.std('time')
            
        spmm_regression_wind = [proc.pointwise_polyfit(spmm_index_wind,ds,1) for ds in dsanoms]
    if verbose:
        proc.printtime(st,"(6/7) Obtained EOF Patterns via regression")
    
    # Step 7. Combine the Output ----------------------------------------------------------
    st = time.time()
    if output_ds:
        # Get SST Regression Patterns
        sst_regression_patterns = xr.concat([spmm_regression[0].coefficients_by_degree.isel(coeff=1).rename("sst"),
                                             spmm_regression[1].coefficients_by_degree.isel(coeff=1).rename("u10"),
                                             spmm_regression[2].coefficients_by_degree.isel(coeff=1).rename("v10")],
                                             dim='variable') #compat='no_conflicts')
        sst_regression_patterns['variable'] = vnames
        outdict = xr.merge([sst_regression_patterns.rename("eof1_sst"),
                            da_sstpcs.rename('pc_sst'),
                            da_windpcs.rename('pc_wind'),
                            da_varexp.rename('variance_explained'),
                           ])
        if regress_wind:
            wind_regression_patterns             = xr.concat([spmm_regression_wind[0].coefficients_by_degree.isel(coeff=1).rename("sst"),
                                                              spmm_regression_wind[1].coefficients_by_degree.isel(coeff=1).rename("u10"),
                                                              spmm_regression_wind[2].coefficients_by_degree.isel(coeff=1).rename("v10")],
                                                             dim='variable')
            wind_regression_patterns['variable'] = vnames
            outdict                              = xr.merge([outdict,wind_regression_patterns.rename("eof1_wind")],compat='no_conflicts') # Add to Dataset
        
            
        
        

    else: # Just Output Dictionaey (note that regression output has more detail then)
        outdict = dict(
            pc_sst   = da_sstpcs,
            pc_wind  = da_windpcs,
            eof_sst  = spmm_regression,
            varexp   = da_varexp,
            vnames   = vnames,
        )
        if regress_wind:
            outdict['eof_wind'] = spmm_regression_wind
    if verbose:
        proc.printtime(st,"(7/7) Prepared output")
        print("SPMM Indices and Patterns calculated in %.2fs" % (time.time()-st_all))
    return outdict

def center_events_ninodict(ensoin,ninodict_in,center_month,search_window,verbose=True,include_full_duration=True):
    
    # Subset events identified in "combine_events" to those peaking within [search_window] 
    # and center on a particular month [center_month].
    
    # (1) Get Indices of Events within [search_window]  ----
    # Indices are relative to the original event list.
    eventmonths = ninodict_in['eventmonths']
    keepevents  = np.where(np.isin(eventmonths,search_window),True,False)
    if include_full_duration:
        
        combined_event_ids = ninodict_in['event_combine']
        ensoidtimes        = ensoin.time
        
        running_count = 0 # Keep count of events
        for ii in range(len(keepevents)):
            keepevent = keepevents[ii] # Check if True or False
            if keepevent == False: # Check if any months within the combined event is within the event
                single_event_times = ensoidtimes[combined_event_ids[ii]].time.dt.month #  Get month of each step within event
                if np.any(np.isin(single_event_times,search_window)):
                    if verbose:
                        print("Event @t=%i has times within search window (%s). Adding to list." % (ninodict_in['center_ids'][ii],search_window))
                        print("\tEvent Times: %s" % single_event_times.data)
                    keepevents[ii] = True
                    running_count += 1
            else:
                if verbose:
                    print("Event Month %i is within search window %s" % (eventmonths[ii],search_window))
        
        print("Kept %i additional events based on full event duration falling within search window" % running_count)
    keepid      = np.where(keepevents)[0] #np.where(np.isin(eventmonths,search_window))[0]
    
    # (2) Center subsetted events on the [center_month]  ----
    # (2.1) Get Closest Distance to Center
    center_ids           = ninodict_in['center_ids']
    eventmonths_keep     = np.array(eventmonths)[keepid]
    dist_to_center_min   = np.array([proc.shortest_distance_mod12(cc,center_month,verbose=False).item() for cc in eventmonths_keep])
    # (2.2) Get Index (from original timeseries) of new event centers, centered on [center_month]
    new_center_ids       = np.array(center_ids)[keepid] 
    corrected_center_ids = new_center_ids + dist_to_center_min 
    # Set to last timestep if it passes
    corrected_center_ids = np.where(corrected_center_ids >= len(ensoin),
                                    len(ensoin)-1,
                                    corrected_center_ids)
    # Set to first timestep if it passes
    corrected_center_ids =  np.where(corrected_center_ids <= 0,0,corrected_center_ids)
    if verbose:
        print("%i of %i events will be kept." % (len(new_center_ids),len(center_ids)))
    
    # (3) Subset original entries of ninodict_in ----
    # (3.1) Subset Event Times
    event_times        = ninodict_in['event_time']
    event_times        = np.array([et.data for et in event_times]) # Extract and Convert to Array
    event_times_subset = event_times[keepid]
    
    # (3.2) Event Sequences (need to convert to object array due to uneven lengths)
    event_combine_subset = np.array(ninodict_in['event_combine'],dtype='object')[keepid]
    
    # (3.3) Build initial dictionary, including new center months + indices + information
    ninodict_subset = dict(
        event_time         = event_times_subset,
        center_ids         = new_center_ids,
        new_center_ids     = corrected_center_ids, # Corrected to index the center month
        center_correction  = dist_to_center_min, # Distance of correction
        event_combine      = event_combine_subset, # 
        center_month       = center_month,
        search_window      = search_window,
        )
    
    # (3.4) Subset other Events (convert to array then subset, add to dict)
    keyloop     = ["event_max","durations","eventmonths"]
    #subset_arrs = []
    for kk in range(len(keyloop)):
        key        = keyloop[kk]
        subset_arr = np.array(ninodict_in[key])[keepid]
        #subset_arrs.append(subset_arr)
        ninodict_subset[key] = subset_arr.copy()
    
    # (3.5) Optional Print to check things
    if verbose:
        keynew = list(ninodict_subset)
        for kk in keynew:
            print(kk)
            print(ninodict_subset[kk])
            print("\n")
                
    return ninodict_subset

def convolve_kernel_ccf(ccfvar,kernel,ccfname,seasonal=False):
    # Convolve radiative kernel with ccf variable 
    # Inputs
    # ccfvar : anomalized CCF              xr.DataArray [time x lat x lon]
    # kernel : coefficients                xr.DataArray [lat x lon], or 
    #                                      LIST of xr.DataArray [seasonal][lat x lon]
    # seasonal : True for seasonal Kernel  BOOL
    # 
    
    selmons_loop = [[12,1,2],[3,4,5],[6,7,8],[9,10,11]]
    
    # Standardize variable
    ccfvar_std = ccfvar / ccfvar.std('time')
    
    if seasonal:
        coeff_seasonal       = [ds.coeffs.sel(ccf=ccfname) for ds in kernel] # /dtday
        R_component_seasonal = []
        for ss in range(4):
            
            selmons    = selmons_loop[ss]
            varmon     = proc.selmon_ds(ccfvar_std,selmons)
            varmon_out = varmon * coeff_seasonal[ss]
            R_component_seasonal.append(varmon_out)
            
        R_component_seasonal = xr.concat(R_component_seasonal,dim='time')
        R_component_seasonal = R_component_seasonal.sortby('time')
        R_component_seasonal = R_component_seasonal.rename(ccfname)
        R_component_seasonal = R_component_seasonal
        return R_component_seasonal
    else:
        
        # Multiply by the Coefficient
        coeff_allmons  = kernel.coeffs.sel(ccf=ccfname)
        R_component    = coeff_allmons * ccfvar_std
        return R_component.rename(ccfname)

def generate_periods(ds,winlen):
    # Generate chunks of [winlen]-years along time dimension in ds
    # Note that winlen must be in YEARS
    tstart = ds.time[0].dt.year.data.item()
    tend   = ds.time[-1].dt.year.data.item()

    nperiods = tend-tstart-winlen+1
    #nperiods
    tranges = []
    subsets = []
    for ii in range(nperiods):
        trange = ["%i-01-01" % (tstart+ii), "%i-01-01" % (tstart+winlen+ii)]
        subset = ds.sel(time=slice(trange[0],trange[1]))
        tranges.append(trange)
        subsets.append(subset)
    return subsets,tranges

def get_early_late_periods(ds_byperiod,num_periods,avg=True):
    # Taken from `sliding_spectra_ccf.py`
    nper = len(ds_byperiod.period)
    # Get Early
    ds_early = ds_byperiod.isel(period=slice(0,num_periods))#.mean('period')
    # Get Late
    ds_late  = ds_byperiod.isel(period=slice(nper-num_periods,nper))#.mean('period')
    if avg:
        ds_early = ds_early.mean('period')
        ds_late  = ds_late.mean('period')
    return ds_early,ds_late

def get_moving_segments(ts,winlen):
    '''
    Divide timeseries into segments using a sliding window of [winlen] months.
    
    Inputs
        ts (xr.DataArray)      : Target timeseries, with time dimension
        winlen (int)           : Number of months in sliding window
    Outputs
        ts_segments (np.array) : Sliding segments with dimensions  [segment_number,winlen]
        center_time (np.array) : Center time of the sliding window [segment_number]

    '''
    # Get Timeseries length and compute maximum number of segments
    ntime = len(ts)
    nsegments = ntime - winlen + 1
    
    # Take Segments
    ts_segments = np.zeros((nsegments,winlen))
    center_time = []#np.zeros((nsegments))
    for i in range(nsegments):
        idseg = np.arange(i,i+winlen)
        ts_segments[i,:] = ts[idseg]
        center_time.append(ts.time.isel(time=int(i+ int(winlen/2))))
    # Concatenate center timesteps
    center_time = xr.concat(center_time,dim='time')
    return ts_segments,center_time

def get_rawpath_awi(expname,vname,ensnum=None):
    # TCo319_ctl1950d
    ctlpath0_31= "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ctl1950d/monthly/"
    ctlpath1_31= "//export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_control/"

    # TCo319_ssp585
    ssppath0_31="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/"
    ssppath0_31_tropics="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/tropics_only/"
    if ensnum is not None:
        ssppath0_31_ens = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585_ens%02i/"
    
    # TCo1279_DART-1950
    ctlpath0_09_atm = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/"
    ctlpath0_09_ocn = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/ocn/mon/"

    # TCo1279_DART-2090
    ssppath0_09     = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/"
    
    # TCo2559-DART-1950C
    ctlpath0_05_atm = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo2559-DART-1950C_5km/atm/"
    ctlpath0_05_ocn = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo2559-DART-1950C_5km/ocn/"


    
    notfound=False
    # ========================================================
    if expname == "TCo319_ctl1950d": # 31 km
        # First Try the Control Path (post 1950 crops)
        searchstr = "%s/*_%s_*.nc" % (ctlpath0_31,vname)
        nclist    = glob.glob(searchstr)
        if len(nclist) > 0:
            print("\tFound in [TCo319_ctl1950d] path!")
            #continue
        else:
            print("\tNot found in [TCo319_ctl1950d] path...")
            # Next Try the older Control Path
            searchstr = "%s/*_%s_*.nc" % (ctlpath1_31,vname)
            nclist    = glob.glob(searchstr)
            if len(nclist) > 0:
                print("\tFound in [TCo319_control] path!")
            else:
                notfound=True
    # ========================================================
    elif expname == "TCo319_ssp585":
        if ensnum is not None:
            searchstr = ssppath0_31_ens % (ensnum) + "*_%s_*.nc" % (vname)
            nclist    = glob.glob(searchstr)
            print("\tFound in [TCo319_ssp585_ens%02i]!" % ensnum)
            #continue
        else:
            print("\tVariable not available for ens members... looking in original folder")
        
        # First Try the Control Path (post 1950 crops)
        searchstr = "%s/*_%s_*.nc" % (ssppath0_31,vname)
        nclist    = glob.glob(searchstr)
        if len(nclist) > 0:
            print("\tFound in [TCo319_ssp585] path!")
            #continue
        else: # Otherwise Try Tropics Only
            print("\tNot found in [TCo319_ssp585] path...")
            # First Try the Control Path (post 1950 crops)
            searchstr = "%s/*_%s_*.nc" % (ssppath0_31_tropics,vname)
            nclist    = glob.glob(searchstr)
            if len(nclist) > 0:
                print("\tFound in [TCo319_ssp585_tropics_only] path!")
                #continue
            else:
                notfound=True
    # ========================================================
    elif expname == "TCo1279-DART-1950":
        # First Try the ATM Path (post 1950 crops)
        searchstr = "%s/*_%s_*.nc" % (ctlpath0_09_atm,vname)
        nclist    = glob.glob(searchstr)
        if len(nclist) > 0:
            print("\tFound in [TCo1279-DART-1950_atm] path!")
            #continue
        else:
            print("\tNot found in [TCo1279-DART-1950_atm] path...")
            # First Try the ATM Path (post 1950 crops)
            searchstr = "%s/*_%s_*.nc" % (ctlpath0_09_ocn,vname)
            nclist    = glob.glob(searchstr)
            if len(nclist) > 0:
                print("\tFound in [TCo1279-DART-1950_ocn] path!")
                #continue
            else:
                notfound=True
    # ========================================================
    elif expname == "TCo1279-DART-2090":
        # First Try the Control Path (post 1950 crops)
        searchstr = "%s/*_%s_*.nc" % (ssppath0_09,vname)
        nclist    = glob.glob(searchstr)
        if len(nclist) > 0:
            print("\tFound in [TCo319_ctl1950d] path!")
        else:
            notfound=True
    # ========================================================
    elif expname == "TCo2559-DART-1950C":
        # First Try the ATM Path (post 1950 crops)
        searchstr = "%s/*_%s_*.nc" % (ctlpath0_05_atm,vname)
        nclist    = glob.glob(searchstr)
        if len(nclist) > 0:
            print("\tFound in [TCo2559-DART-1950C_atm] path!")
            #continue
        else:
            print("\tNot found in [TCo2559-DART-1950C_atm] path...")
            # First Try the ATM Path (post 1950 crops)
            searchstr = "%s/*_%s_*.nc" % (ctlpath0_05_ocn,vname)
            nclist    = glob.glob(searchstr)
            if len(nclist) > 0:
                print("\tFound in [TCo2559-DART-1950C_ocn] path!")
                #continue
            else:
                notfound=True
    # ========================================================
    else:
        print("%s not recognized..." % expname)
        notfound=True
    if notfound:
        print("%s Not found for %s" % (vname,expname))
        return None
    print(nclist)
    return nclist

def get_center_time(trange):
    
    # Return Jan-01 for the middle year inbetween two dates
    
    # for YYYY-MM-DD
    ystart  = int(trange[0][:4])
    yend    = int(trange[1][:4])
    ymid    = int(np.round((ystart + yend) / 2))
    tcenter = "%04i-01-01" % (ymid)
    
    return tcenter

def get_tcenters(tranges):
    # Given list of [tstart,tend], get center times and year string for plotting
    tcenters = [get_center_time(tt) for tt in tranges]
    years   = [tt[:4] for tt in tcenters]
    return tcenters,years

def init_tp_map(nrow=1,ncol=1,figsize=(12.5,4.5),ax=None,latmax=20,lonbounds=[120,290],):
    bbplot = [lonbounds[0], lonbounds[1], -latmax, latmax]
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

def init_globalmap(nrow=1,ncol=1,figsize=(12,8),centlon=200):
    proj            = ccrs.Robinson(central_longitude=centlon)
    #bbox            = [-180,180,-90,90]
    fig,ax          = plt.subplots(nrow,ncol,subplot_kw={'projection':proj},figsize=figsize,constrained_layout=True)
    
    multiax = True
    if (type(ax) == mpl.axes._axes.Axes) or (type(ax) == cartopy.mpl.geoaxes.GeoAxes):
        ax = [ax,]
        multiax = False
    
    if type(ax) == tuple or (ncol+nrow > 2):
        ax = ax.flatten()
    for a in ax:
        a.coastlines(zorder=10,lw=0.75,transform=proj)
        a.gridlines(ls ='dotted',draw_labels=True)
        
    if multiax is False:
        ax = ax[0]
    return fig,ax

def init_sep_map(nrow=1,ncol=1,figsize=(6,8),expandx=10,expandy=10):
    bbox_sep    = [-90+360,-75+360,-40,-15] # Southeast Tropical Pacific Box from Kang et al. 2026
    bbox_plot   = proc.expand_bbox(bbox_sep,expandy,expandx)
    fig,axs     = plt.subplots(nrow,ncol,subplot_kw={'projection':ccrs.PlateCarree()},constrained_layout=True,figsize=figsize)

    # Set up the axes
    if nrow == 1 and ncol == 1:
        axs.set_extent(bbox_plot)
        axs = viz.add_coast_grid(axs,bbox=bbox_plot,line_color='lightgray',grid_color="w",ignore_error=True)
    else:
        for ax in axs:
            ax.set_extent(bbox_plot)
            ax = viz.add_coast_grid(ax,bbox=bbox_plot,line_color='lightgray',grid_color="w",ignore_error=True)
    return fig,axs


def init_sliding_spectra_plot(figsize=(8,6),xshift_text=0.9):
    # Taken from `sliding_spectra_other_simulations`
    # This version was copied from sliding_spectra_ccf_regional..ipynb
    fsz_ticks = 12
    fsz_axis  = 12
    
    fig,ax          = plt.subplots(1,1,constrained_layout=True,figsize=figsize)
    
    
    # Setup Parameters ------------------------
    
    xper            = np.array([20,10,5,2,1,0.5])
    xper_ticks      = 1 / (xper*12)
    
    # Plot Combination Tones (Sum/Difference)
    x_min, x_max    = ax.get_xlim()
    Mfreq_enso      = np.array([1/(2*12), 1/(5.5*12)])
    Mfreq_fplus     = 1/12 + Mfreq_enso
    Mfreq_fmins     = 1/12 - Mfreq_enso
    
    ax.axhline([1/(5.5*12)],label="",ls='dotted',c='lightgray')
    ax.axhline([1/(2*12)],label="",ls='dotted',c='hotpink')
    ax.axhline([1/(15)],label="",ls='dotted',c='hotpink')
    ax.axhline([1/(8)],label="",ls='dotted',c='yellow')
    ax.axhline([1/(10)],label="",ls='dotted',c='yellow')
    
    # Twin Y-Axis for Period Labels
    ax.set_ylabel("Frequency ($month^{-1}$)")
    #ax.set_yticks(1/xper_ticks)
    ax.set_ylim([xper_ticks[0],0.5])
    ax.set_yscale('log')
    
    ax2 = ax.twinx()
    ax2.set_ylim([xper_ticks[0],0.5])
    ax2.set_yscale('log')
    ax2.set_yticks(xper_ticks,labels=xper)
    ax2.set_ylabel("Period (Years)",fontsize=fsz_axis)
    
    #xshift_text = 12
    ax.text(x_min+xshift_text,np.mean(Mfreq_enso)*0.9,'$f_{ENSO}=$ 2-5.5 yrs', ha='center', va='top',c='w')
    ax.text(x_min+xshift_text,np.mean(Mfreq_fplus)+0.1/12, '$1+f_{ENSO}$ 8-10 mons', ha='center', va='top',c='yellow')
    ax.text(x_min+xshift_text,np.mean(Mfreq_fmins), '$1-f_{ENSO}$ =  15-24 mons', ha='center', va='top',c='hotpink')
    # ===========
    return fig,ax,ax2

def load_ccf_kernel(expname,flxname,customname=None,standardize=True,seasonal=False,kernel_path=None):
    # load Kernels for each CCF (based on calc_ccf_radiation_byexp)
    if kernel_path is None:
        kernel_path  = '/home/niu4/gliu8/projects/ccfs/kernels/regrid_1x1/%s/' % expname
    
    selmons_loop = [[12,1,2],[3,4,5],[6,7,8],[9,10,11]]
    
    if seasonal:
        ds_byseason = []
        for selmons in selmons_loop:
            
            if customname is not None:
                ncname_kernel         = "%s%s_%s_kernels_standardize%i.nc" % (kernel_path,flxname,customname,standardize)
            else:
                ncname_kernel         = "%s%s_kernels_standardize%i.nc" % (kernel_path,flxname,standardize)
            if selmons is not None:
                selmonstr       = proc.mon2str(np.array(selmons)-1)
                ncname_kernel   = proc.addstrtoext(ncname_kernel,"_"+selmonstr,adjust=-1)
            
            # Load Seasonal Outputs
            ds = xr.open_dataset(ncname_kernel).load()
            ds_byseason.append(ds)
        return ds_byseason
        
    else:
        #% Load the all months case
        if customname is not None:
            ncname_kernel         = "%s%s_%s_kernels_standardize%i.nc" % (kernel_path,flxname,customname,standardize)
        else:
            ncname_kernel         = "%s%s_kernels_standardize%i.nc" % (kernel_path,flxname,standardize)
        dsall                     = xr.open_dataset(ncname_kernel).load()
        
        return dsall

def load_ccf_radiation(expname,flxname,datpath=None,seasonal=False):
    # Load CCF Radiation. Taken from 'awi_cm3_toa_leadlag_analysis_area_avg_ccf_sliding.ipynb'
    
    if datpath is None:
        datpath = "/home/niu4/gliu8/projects/ccfs/radiative_components/regrid_1x1/"
    
    ccfs = ["sst","eis","Tadv","r700","w700","ws10"]
    
    dscomps = []
    for cc in range(6):
        
        ccf_var  = ccfs[cc]
        if seasonal:
            ncname = "%s%s/%s_%s_component_seasonal.nc" % (datpath,expname,flxname,ccf_var)
        else:
            ncname   = "%s%s/%s_%s_component.nc" % (datpath,expname,flxname,ccf_var)
        dsccf    = xr.open_dataset(ncname)[ccf_var].load()

        dscomps.append(dsccf)
    return dscomps

def load_ensoid(expname,ninoid_name='nino34',datpath=None,standardize=True):
    # Load Enso indices calculated with calc_nino34.py
    if datpath is None:
        datpath = "/home/niu4/gliu8/projects/scrap/nino34/"
    ncname = "%s%s_%s.nc" % (datpath,expname,ninoid_name)
    ds = xr.open_dataset(ncname).load()
    if standardize:
        return ds.sst
    return ds.sst * ds['std'].data.item()

def load_enso_eof(expname,datpath=None,apply_smoothing=True,sep_mon=False,by_period=False,winlen=None):
    # Load EOF-based ENSO computed via `calc_EOF_ENSO.py`
    # by_period loads for sliding period [winlen], currently does not support [sep_mon] `calc_EOF_ENSO_sliding.py`
    
    if by_period:
        print("Loading ENSO computed for sliding windows of %i-years" % winlen)
        if datpath is None:
            datpath       = "/home/niu4/gliu8/projects/ccfs/enso_eof/"
        ninonc        = "%s%s_ENSO_EOF_slidingwinlen%02i.nc" % (datpath,expname,winlen)

    else: # Load for calculation over whole timeseries
        if datpath is None:
            datpath = "/home/niu4/gliu8/projects/scrap/nino34/"
        ninonc        = "%s%s_enso_eof_rotated.nc" % (datpath,expname)
        if sep_mon:
            ninonc    = proc.addstrtoext(ninonc,"_sepmon",adjust=-1)
    dsnino        = xr.open_dataset(ninonc).load()
    
    # Apply 1-2-1 filter
    if apply_smoothing:
        
        filter_coeffs = [0.25,0.5,0.25]
        smooth121     = lambda timeseries: np.convolve(timeseries,filter_coeffs,mode='same')
        
        if by_period:
            
            ep = xr.apply_ufunc(smooth121, dsnino.ep, input_core_dims=[["timeindex"]],output_core_dims=[["timeindex"]],vectorize=True)
            cp = xr.apply_ufunc(smooth121, dsnino.cp, input_core_dims=[["timeindex"]],output_core_dims=[["timeindex"]],vectorize=True)
            ninotimes = [make_ninotime(trange,dsnino.timeindex) for trange in dsnino.trange] 
        else:
            ep     = np.convolve(dsnino.ep,filter_coeffs,mode='same')
            cp     = np.convolve(dsnino.cp,filter_coeffs,mode='same')
            tcoord = dict(time=dsnino.time)
            ep     = xr.DataArray(ep,coords=tcoord,dims=tcoord,name='ep')
            cp     = xr.DataArray(cp,coords=tcoord,dims=tcoord,name='cp')
        
    else:
        ep        = dsnino.ep
        cp        = dsnino.cp
        
    if by_period:
        ninotimes = [make_ninotime(trange,dsnino.timeindex) for trange in dsnino.trange] 
        return ep,cp,ninotimes
    return ep,cp

def load_land_mask_awi(expname,regrid=False,outpath=None):
    if outpath is None:
        outpath = '/home/niu4/gliu8/projects/scrap/awi_common/'
    if "TCo319" in expname: # 31km Simulations
        print("Loading for 31 km simulations...")
        if regrid: # Load 180x360
            dsmask = "TCo319_ctl1950d_r360x180_landmask.nc"
        else: # Load 640x1213
            dsmask = "TCo319_atmgrid_original_landmask.nc"#"TCo319_ssp585_atm_landmask.nc"
    elif "TCo1279" in expname: # 9km Simulations
        print("Loading for 9 km simulations...")
        if regrid:
            dsmask = "TCo1279_DART-1950_ocn_r3600x1800_landmask.nc"
        else:
            dsmask = "TCo1279_atmgrid_original_landmask.nc"#"TCo1279-DART-1950_atm_landmask.nc"
    elif "TCo2559" in expname: # 5km simulations
        print("Loading for 5 km simulations...")
        if regrid:
            dsmask = "TCo2559-DART-1950C_ocn_r3600x1800_landmask.nc"
        else:
            dsmask = "TCo2559_atmgrid_original_landmask.nc"
    else:
        print("Experiment not found")
        return np.nan
    return xr.open_dataset(outpath+dsmask).land_mask.load()

def loadregrid(expname,vname,reformat=False,bbox=None,check_var=True):
    
    # load regridded 1x1 product on niu (convenience function)
    # Check CERES First
    not_in_ceres_flag = True # True if variable has not been found in ceres
    if "CERES" in expname:
        datpath = "/home/niu4/gliu8/share/CERES/processed/"
        if "EBAF" in expname and np.isin(vname,['allsky','clearsky','cre','tscre','ttcre','tsr','tsrc','ttr','ttrc']): # EBAF (Fluxes Only)
            ncname = "%sCERES_EBAF_%s_2000-03_to_2025-08.nc" % (datpath,vname)
            not_in_ceres_flag = False
        elif "FBCT" in expname or vname in ['lcc','tcc','ucc']: # FBCT, includes cloud fraction
            ncname = "%sCERES_FBCT_%s_2002-07_to_2023-02.nc" % (datpath,vname)
            not_in_ceres_flag = False
        
    if "ERA5" in expname and not_in_ceres_flag:
        datpath = "/home/niu4/gliu8/projects/common_data/ERA5/regrid_1x1/"
        ncname  = "%s%s_1979_2024.nc" % (datpath,vname)
    elif not_in_ceres_flag:
        datpath = "/home/niu4/gliu8/projects/scrap/regrid_1x1/"
        ncname  = "%s%s_%s_regrid1x1.nc" % (datpath,expname,vname)
    
    # Open view of dataset
    ds      = xr.open_dataset(ncname)[vname]
    
    # Load and reformat, if option is set...
    if reformat:
        print("Loading and Reformatting DataSet...")
        if bbox is not None:
            ds = proc.sel_region_xr(ds,bbox) # Select a region
        ds = ds.load()
        ds = standardize_names(ds)
        if check_var is True: # Do additional unit conversions for AWI-CM3/Other Output
            ds = varcheck(ds,vname,expname)
        
        return ds
    return ds

def load_scott2020_kernels(vname="dRdxi",reformat=True):
    # Adapted from `visualize_kernels_ERA5_Scott_etal.ipynb`
    # Loading Inputs
    spath = "/home/niu4/gliu8/projects/ccfs/Scott_etal_2020/meteorological_cloud_radiative_kernels-main/obs/"
    ncnames = [
        "Scott_Myers_meteorological_kernels_CERES-FBCT_ERA5.nc",
        "Scott_Myers_meteorological_kernels_ISCCP-H_ERA5.nc",
        "Scott_Myers_meteorological_kernels_MODIS_ERA5.nc",
        "Scott_Myers_meteorological_kernels_PATMOS-x_ERA5.nc",
    ]
    ds2020_names = [
        "CERES-FBCT",
        "ISCCP-H",
        "MODIS",
        "PATMOS-x",  
    ]
    
    
    # Load Datasets
    ds2020 = [xr.open_dataset(spath+nc).load() for nc in ncnames]
    
    if reformat:
        # Get the Kernels
        ds2020_kernels = [ds[vname] for ds in ds2020] 
        # Adjust formatting
        def reformat_ds(ds):
            # Assign CCFs to Dimensions, Standardize Lat/Lon Names
            ccf_scott      = ["sst", "eis", "Tadv", "r700", "w700", "ws10"]
            ds             = ds.rename(dict(i='ccf'))
            ds['ccf']      = ccf_scott
            ds             = standardize_names(ds)
            return ds
        ds2020_kernels_rfmt = [reformat_ds(ds) for ds in ds2020_kernels]
        # Concatenate by Dataset
        ds2020_kernels_rfmt = xr.concat(ds2020_kernels_rfmt,dim='dataset')
        ds2020_kernels_rfmt['dataset'] = ds2020_names
        return ds2020_kernels_rfmt
    else:
        return ds2020,ds2020_names
    
def make_ninotime(trange,timeindex):
    ntime = len(timeindex)
    times = xr.date_range(start=trange[0].data.item(),end=trange[1].data.item(),periods=ntime)
    return times # Note this does not return the exact middle of the month

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

def mlr(X,y):
    # MLR fit using scipy 
    
    # Initialize Model and Fit
    model             = LinearRegression()
    model.fit(X,y)
    pred              = model.predict(X)
    # Calculate Error and other variables
    mlr_out           = {}
    mlr_out['pred']   = pred
    mlr_out['err']    = y - pred # Model Error # []
    mlr_out['coeffs'] = model.coef_
    mlr_out['r2']     = sklearn.metrics.r2_score(y,pred)
    
    return mlr_out

def mlr_ccfs(ccfs,flx,standardize=True,fill_value=0,verbose=False):
    # Perform MLR
    #    ccfs: LIST of DataArrays [variable x time]
    #    flx:  DataArray [time x 1]
    
    # Set up Predictors and Target (convert to Numpy Arrays)
    predictors = np.array([ds for ds in ccfs]) # [variable x time]
    if standardize:
        if verbose:
            print("Standardizing each variable")
        predictors = np.array([ds/np.nanstd(ds) for ds in list(predictors)])
    X = predictors.T
    y = flx.data
    
    # Replace NaN Values in Predictors
    if verbose:
        if np.any(np.isnan(X)):
            print("NaN values detected! Replace with %f" % fill_value)
    X = np.where(np.isnan(X),fill_value,X) # Set NaN to zero
    
    # Use sklearn for now (can try LSE manual later...)
    mlr_out = mlr(X,y)
    return mlr_out

def preprocess_byperiod(dswins,order=1,verbose=False,detrend=True):
    # Input: List of Arrays where elements are DataArrays for each period
    nwin    = len(dswins)
    dsanoms = []
    for nw in range(nwin):
        dsin   = dswins[nw].squeeze()
        dsanom = proc.xrdeseason(dsin,verbose=verbose)
        if detrend: # Optionally Detrend
            if len(dsanom.shape) > 1:
                dsanom = proc.xrdetrend_nd(dsanom,order,verbose=verbose)
            else:
                dsanom = proc.xrdetrend(dsanom,verbose=verbose)
        dsanoms.append(dsanom)
    return dsanoms

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

def remove_duplicate_times(ds,verbose=True,timename='time'):
    # From : https://stackoverflow.com/questions/51058379/drop-duplicate-times-in-xarray
    _, index = np.unique(ds[timename], return_index=True)
    print("Found %i duplicate times. Taking first entry." % (len(ds[timename]) - len(index)))
    return ds.isel({timename:index})

def shift_time_monthstart(dsin,timename='time'):
    # Shift dates from center of month to 1st of month (09-15 --> 09-01)
    # example: calculate radiative kernel.py
    oldtime         = dsin[timename]
    tstart          =  str(oldtime[0].data)[:7] + "-01"
    tend            =  str(oldtime[-1].data)[:7] + "-01"
    newtime         = pd.date_range(start=tstart,end=tend,freq="MS")
    dsin[timename]    = newtime
    print("New time dimension between %s and %s" % (dsin[timename][0].data,dsin[timename][-1].data))
    return dsin


def sliding_spectra(dsraw,nyr_window,nsmooth,detrend=False):
    # Adapted from sliding_pointwise_spectra
    # Take in dsraw, subset into periods, detrend/deseason, then compute spectra
    # [dsraw]           : raw 1-D data-array containing variable with dimenion 'time'
    # [nyr_window]      : sliding window length in years
    # [nsmooth]         : amount adjacent bands to smooth across 
    # [detrend]         : set to true to apply detrend
    # 
    # Uses following funcs from utils: 
    #   generate_periods, get_center_time,preprocess_byperiod
    
    subsets,tranges = generate_periods(dsraw,nyr_window)
    nperiods        = len(subsets) # Get Number of Periods
    tcenters        = [get_center_time(t) for t in tranges]
    
    # Preprocess by Period
    subsets_anom   = preprocess_byperiod(subsets,verbose=False,detrend=detrend) # Don't Detrend...
    
    # Calculate Power Spectra for each period
    spec_byperiod  = []
    for nw in range(nperiods):
        sample_in            = subsets_anom[nw].data
        specout              = proc.point_spectra(sample_in,nsmooth=nsmooth,return_conf=True)
        specout['tcenter']   = tcenters[nw]
        spec_byperiod.append(specout)
    
    # Concatenate By period
    spec_byperiod             = xr.concat(spec_byperiod,dim='period')
    
    return spec_byperiod

 
def standardize_names(ds):
    
    # Time Names
    ds = swap_rename(ds,'valid_time','time')
    ds = swap_rename(ds,'time_counter','time') # ORAS5
    ds = swap_rename(ds,"TIME_COUNTER",'time')
    
    # Longitude Names
    ds = swap_rename(ds,"longitude","lon")
    ds = swap_rename(ds,"LON","lon")
    ds = swap_rename(ds,"nav_lon","lon") # ORAS5
    
    # Latitude Names
    ds = swap_rename(ds,"latitude","lat")
    ds = swap_rename(ds,"LAT","lat")
    ds = swap_rename(ds,"LAT232_409","lat")
    ds = swap_rename(ds,"nav_lat","lat") # ORAS5
    
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
            
        elif iend  >= ntime:
            filler = np.zeros( (iend-ntime+1)) * np.nan
            subset = np.hstack([target_var[istart:],filler])
            stacked_events[ie,:] = subset
        elif istart < 0: # Note havent tested this
            filler  = np.zeros(np.abs(istart)) * np.nan
            subset  = np.hstack([filler,target_var[:(iend+1)],])
            stacked_events[ie,:] = subset
    return stacked_events

def stack_events_2d(invar,eventids,ibefore,iafter,times_da=None):
    
    invar           = invar.transpose('time','lat','lon')
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
        
        if iend >= ntime:
            print("iend is at %s" % (iend))
            insert_end = nlags + (ntime-iend) -1 #(lags*2+1) - (iend - ntime)
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
    
    # Performs the following conversions (for AWI-CM3 Output)
    # all fluxes: Convert to w/m2 (considering accumulation period)
    # sst: Convert to celsius
    # pr,lsp,cp: Convert to mm/day
    # msl: Convert to hPa
    # w700,w500: Convert to hPa/day
    
    if type(ds) == xr.Dataset:
        print("Converting to DataArray!")
        ds = ds[vname]
    if np.any(ds > 273) and vname == "sst": # Convert to Celsius
        print("Converting from Kelvin to Celsius for %s" % expname)
        ds = ds - 273.15
        
    
    # For Heat Fluxes, convert to W/m2 for AWI simulations considering the
    # accumulation time: 
    # 6h for 31km (TCo319)
    # 3h for 9km (TCo1279) and 5km (TCo2559)
    if vname in ['str','ssr',
                 'strc','ssrc',
                 'ttr','tsr',
                 'ttrc','tsrc',
                 'ttcre','tscre',
                 'allsky','clearsky','cre'
                 'sshf','slhf']: # Accumulation over 3h
        if "TCo" in expname:
            # Conversion for STR and SSR considering 3h Accumulation
            if "TCo319" in expname:
                print("Correction for accumulation over 6 hours for %s" % expname)
                accumulation_hr = 6
            elif ("TCo1279" in expname) or ("TCo2559" in expname):
                print("Correction for accumulation over 3 hours for %s" % expname )
                accumulation_hr = 3
            conversion  = 1/(3600 * accumulation_hr)  # 3 h accumulation time...? #1/(24*30*3600)
            # https://forum.ecmwf.int/t/surface-radiation-parameters-joule-m-2-to-watt-m-2/1588
            ds          = ds * conversion
        else:
            # No Conversion is Done
            print("%s unchanged ([TCo] not detected in [expname])" % vname)
            
    
    # Accounting for accumulation time as above in AWI-CM3
    # Convert from [meters/accumulation period] to [mm/day]
    if vname in ["cp","lsp","pr"]:
        if "TCo" in expname:
            if "TCo319" in expname:
                print("Correction for accumulation over 6 hours for %s" % expname)
                accumulation_hr = 6 
            elif ("TCo1279" in expname) or ("TCo2559" in expname):
                print("Correction for accumulation over 3 hours for %s" % expname )
                accumulation_hr = 3
            # From Discussion with Sun-Seon
            # the value is provided as the accumulated precipitation (in meters) over a N-hour period.
            # To convert to mm/day.... 
            # (1) Divide by accumulation period (nday * sec * min) to convert to precip rate (m/s)
            # (2) Multiply by day (m/day)
            # (3) Multiply to get milimeters (mm/day)
            #multiply by the number of seconds in a day (86400 s), and convert meters to millimeters (×1000):
            nsec_perday = 86400 
            conversion = (1/(accumulation_hr*60*60)) * nsec_perday * 1000#(24/accumulation_hr) * 1000
            ds         = ds * conversion
        else:
            # No Conversion is Done
            print("%s unchanged ([TCo] not detected in [expname])" % vname)
    
    if vname in ['msl']:
        if "TCo" in expname or np.isin(expname,["ERA5"]):
            print("Converting from Pa to hPa for %s" % vname)
            ds = ds/100 # Convert from Pa to hPa
    if vname in ['w700','w500']:
        if "TCo" in expname:
            print("Converting from Pa/sec to hPa/day for %s" % vname)
            ds = ds/100 * 3600*24 # Convert to hPa/day
    
    return ds



