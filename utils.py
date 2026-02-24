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
calc_grad_centered          : (g) Calculate centered-difference for spatial gradients
calc_lag_regression_1d      : (g) Compute lead lag regression for 1d timeseries
calc_leadlag_regression_2d  : (g) Compute lead lag regression for 2D timseries (modeled after enso_lag_regression)
combine_events              : (g) Given identified events, combine similar events and get other traits (duration, etc)
get_moving_segments         : (g) Subset timeseries into segments with a sliding window
get_rawpath_awi             : (A) Get rawpath for AWI output on niu
init_tp_map                 : (v) initialize tropical Pacific plot 
init_global_map             : (v) Initialize a global map
load_ensoid                 : (l) load enso indices calculated by calc_nino34.py
load_land_mask_awi          : (l) load land mask from AWI-CM3, where land points are np.nan
mcsample                    : (g) Monte Carlo Sampler to repeat function
mlr                         : (g) single point multiple linear regression using Scipy
mlr_ccfs                    : (c) Multiple linear regrssion for cloud controlling factor analysis
preprocess_enso             : (c) detrend (quadratic) and deseasonalize for ENSO calculations
remove_duplicate_times      : (g) Remove duplicate times from a DataArray
shift_time_monthstart       : (g) Shift times from center to start of month
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

def init_globalmap(nrow=1,ncol=1,figsize=(12,8)):
    proj            = ccrs.Robinson(central_longitude=-180)
    bbox            = [-180,180,-90,90]
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

def calc_lag_regression_1d(var1_lag,var2_base,lags,correlation=False): # CAn make 2d by mirroring calc_lag_covar_annn
    # Calculate the regression where
    # (+) lags indicate var1 lags  var2 (var 2 leads)
    # (-) lags indicate var1 leads var2 (var 1 leads)
    
    if np.any(np.isnan(var1_lag)) or np.any(np.isnan(var2_lag)):
        print("NaN detected. Returning NaN...")
        return np.nan * np.ones(len(lags*2)-1)
    
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

def load_ensoid(expname,ninoid_name='nino34',datpath=None,standardize=True):
    # Load Enso indices calculated with calc_nino34.py
    if datpath is None:
        datpath = "/home/niu4/gliu8/projects/scrap/nino34/"
    ncname = "%s%s_%s.nc" % (datpath,expname,ninoid_name)
    ds = xr.open_dataset(ncname).load()
    if standardize:
        return ds.sst
    return ds.sst * ds['std'].data.item()

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
 
def standardize_names(ds):
    
    ds = swap_rename(ds,'valid_time','time')
    ds = swap_rename(ds,'time_counter','time')
    ds = swap_rename(ds,"TIME_COUNTER",'time')
    ds = swap_rename(ds,"LON","lon")
    ds = swap_rename(ds,"LAT","lat")
    ds = swap_rename(ds,"longitude","lon")
    ds = swap_rename(ds,"latitude","lat")
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
    
    if type(ds) == xr.Dataset:
        print("Converting to DataArray!")
        ds = ds[vname]
    if np.any(ds > 273) and vname == "sst": # Convert to Celsius
        print("Converting from Kelvin to Celsius for %s" % expname)
        ds = ds - 273.15
    
    if vname in ['str','ssr',
                 'strc','ssrc',
                 'ttr','tsr',
                 'ttrc','tsrc',
                 'ttcre','tscre',
                 'allsky','clearsky','cre',
                 'sshf','slhf']: # Accumulation over 3h
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



