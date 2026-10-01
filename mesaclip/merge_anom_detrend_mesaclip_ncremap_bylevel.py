#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Merge across BHIST and RCP85 for Regridded MESACLIP Output.
Also remove mean seasonal cycle and the [deg]-order trend.

Loop Version by Level (2026.08.17)

Copied `merge_anom_detrend_mesaclip.py` but adopt to output from ncremap script
by Cuong

Created on Mon Aug  3 15:02:32 2026

@author: gliu

"""

import time
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import glob
import os
from tqdm import tqdm

#%% Helper Functions
def remove_duplicate_times(ds,verbose=True,timename='time'):
    # From : https://stackoverflow.com/questions/51058379/drop-duplicate-times-in-xarray
    _, index = np.unique(ds[timename], return_index=True)
    print("Found %i duplicate times. Taking first entry." % (len(ds[timename]) - len(index)))
    return ds.isel({timename:index})

def detrend_dim(da, dim="time", deg=1):
    # Function by Rohit Ghosh: https://github.com/rg568/EERIE_scripts/blob/main/FESOM/ENSO_Z500_DJF_teleconnection_IFS-FESOM.ipynb
    coeffs = da.polyfit(dim=dim, deg=deg)
    trend = xr.polyval(da[dim], coeffs.polyfit_coefficients)
    return da - trend

def dailyclim(ds):
    return ds.groupby("time.dayofyear").mean("time")

def deseason_daily(ds,clim=False,period=None,verbose=True):
    if period is not None:
        if verbose:
            print("Taking climatology over the provided period: %s" % period)
        dsperiod = ds.sel(time=slice(*period))
        dsclim   = dailyclim(dsperiod)
    else:
        dsclim   = dailyclim(ds)#ds.groupby("time.dayofyear").mean("time")
    dsanom       = ds.groupby("time.dayofyear") - dsclim
    dsanom       = dsanom.drop_vars('dayofyear')
    if clim:
        return dsanom,dsclim
    return dsanom

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



# Copied Path functions from crop_TP_MESACLIP
def get_paths(expname,ens,freq="month_1",scenario="BHIST",regrid=True,realm='atm'):
    # Need to add scenario to the data       
    # First get the path to the data
    # Scenarios Supported" BHIST, BRCP85
    if expname == "hires": # Hi-Res Runs in Campaign
        dpath = "/glade/campaign/collections/cmip/CMIP6/CESM-HR/RDA/%s" % scenario
        gridname="ne120_t12"
    elif expname == "lores":
        gridname="ne30_g16"
        if ens <= 10:  # First 10 Ens Lo-Res
            dpath = "/glade/campaign/collections/cmip/CMIP6/CESM-HR/RDA/lowres/%s" % scenario
        else:   # Ens 11-40 by Cuong
            dpath = "/glade/campaign/cgd/cas/scsan/CESMLR_New_Run"
    else:
        print("expname must be [hires] or [lores]")
    
    
    if regrid: # Set path to data regridded by ncremap
        if expname == "lores":
            dpath = "/glade/derecho/scratch/glennliu/MESACLIP/regridded/lores"
        else:
            dpath = "/glade/derecho/scratch/glennliu/MESACLIP/regridded"
    
    # Next, get the experiment string (based on scenario)
    if scenario == "BHIST":
        if ens == 1:
            ystr = "cesm-ihesp-sehires38-1850-2005"
        else: 
            # Hi-Res Has Specific Numbers
            if expname == 'hires':
                if (ens > 1) and (ens < 4): # 2-3
                    ystr = "cesm-ihesp-hires1.0.30-1920-2005"
                elif (ens >= 4) and (ens < 6): # 4-5
                    ystr = "cesm-ihesp-hires1.0.44-1920-2005"
                elif ens == 6: # 6
                    ystr = "cesm-ihesp-hires1.0.45-1920-2005"
                elif (ens > 6) or (ens < 11): # 6-10
                    ystr = "cesm-ihesp-hires1.0.46-1920-2005"
                else:
                    print("Warning! ens not recognized. For Hi Res. Ens is between 1 and 10")  
            elif expname == 'lores':  # Lo-Res is either 42 or 46 ====================
                if ens <= 10:
                    ystr = "cesm-ihesp-hires1.0.42-1920-2005"
                else:
                    ystr = "cesm-ihesp-hires1.0.46-1920-2005"
        
        # Now Combine
        if regrid:
            exppath = "%s/b.e13.BHISTC5.%s.%s.%03i/%s/%s/" % (dpath,gridname,ystr,ens,realm,freq)
        else:
            exppath = "%s/b.e13.BHISTC5.%s.%s.%03i/%s/proc/tseries/%s" % (dpath,gridname,ystr,ens,realm,freq)
    elif scenario == "BRCP85":
        if ens == 1:
            ystr = "cesm-ihesp-sehires38-2006-2100"
        else:
            if expname == 'hires':
                if ens == 2:   # (Ens 2)
                    ystr = "cesm-ihesp-hires1.0.30-2006-2100"
                elif ens == 3: # (Ens 3)
                    ystr = "cesm-ihesp-hires1.0.31-2006-2100"
                elif (ens >= 4) and (ens < 6): # (Ens 4-5)
                    ystr = "cesm-ihesp-hires1.0.44-2006-2100"
                else:          # (Ens 6-10)
                    ystr = "cesm-ihesp-hires1.0.46-2006-2100"
            elif expname == "lores":
                if ens <= 10: # (Ens 1-10) on gdex
                    ystr = "cesm-ihesp-hires1.0.42-2006-2100"
                else:
                    ystr = "cesm-ihesp-hires1.0.46-2006-2100"
        # Now Combine
        if regrid:
            exppath = "%s/b.e13.BRCP85C5.%s.%s.%03i/%s/%s/" % (dpath,gridname,ystr,ens,realm,freq)
        else:
            exppath = "%s/b.e13.BRCP85C5.%s.%s.%03i/%s/proc/tseries/%s" % (dpath,gridname,ystr,ens,realm,freq)
    return exppath

def get_nclist_ens(expname,vname,freq,scenario,realm='atm',regrid=False,debug=False):
    # Convenience Function to get all NetCDFs
    # Needs to be updated once permissions are changed for Lo_Res Run
    if expname == "hires":
        ensall = np.arange(1,11,1)
    else:
        ensall = np.arange(1,41,1) # np.arange(1,41,1)#
    
    ncall = []
    print("Searching for NetCDFs for %s" % expname)
    for ens in ensall:
        datpath  = get_paths(expname,ens,freq=freq,scenario=scenario,realm=realm,regrid=regrid)
        #print(datpath)
        ncsearch = "%s/*%s*.%s.*.nc" % (datpath,scenario,vname)
        if debug:
            print(ncsearch)
        nclist   = glob.glob(ncsearch)
        nclist.sort()
        nfiles   = len(nclist)
        print("\tFound %2i files for ens %03i..." % (nfiles,ens))
        #print(nclist)
        ncall.append(nclist)
    return ncall,ensall

# ===============================================================
#%% User Edits
# ===============================================================

tstart  = "1982-01-01"
tend    = "2025-12-31" #"2025-12-31"
outpath = "/glade/derecho/scratch/glennliu/MESACLIP/processed/"
expname = "lores"
vname   = "TEMP"
freq    = "month_1"
deg     = 1
climperiod = ['1985-01-01','2014-12-31'] # Set Period to Calculate Climatology

enslist_restrict = None #np.arange(2,11,1)

if expname == "lores":
    rawpath =  "/glade/derecho/scratch/glennliu/MESACLIP/regridded/lores" #% (scenario)
else:
    rawpath = "/glade/derecho/scratch/glennliu/MESACLIP/regridded" #% (scenario)

# NOTE From Here, Need to Add in the Special Characters Depending on Scenario...

# Set Level Settings
levels       = [50,100,200,500,1000]
#levels       = #[100,]#50,200,500,1000]
nlvls        = len(levels)

# ===============================================================
#%% Get Number of Ensemble Members, Make Output Directory
# ===============================================================

tstart       = tstart.replace('-','')#npdatetime_to_str(dsanom.time[0]).replace('-','')
tend         = tend.replace('-','')  #npdatetime_to_str(dsanom.time[-1]).replace('-','')
outpath_proc = "%s/anom_detrend%i_%s-%s/" % (outpath,deg,tstart,tend,)
if climperiod is not None:
    climname     = "_climatology%sto%s" % (climperiod[0][:4],climperiod[1][:4])
    outpath_proc = outpath_proc[:-1] + climname + "/"
makedir(outpath_proc)

if enslist_restrict is not None:
    print("Setting Manual Ens List")
    enslist = enslist_restrict
else:
    if expname == "lores":
        enslist = np.arange(1,41,1)
    else:
        enslist = np.arange(1,11,1)
        
nens = len(enslist)

for ll in range(nlvls):
    level        = levels[ll]
    vname_in     = "%s%s" % (vname,level)
    print("Processing for %s" % vname_in)
    
    #%% Looping for Each Ensemble Member
    
    # Get NetCDFs By Scenario
    nclists_byscenario = [] # [scenario][ens][file]
    for scenario in ["BHIST","BRCP85"]:
        nclists,ensall = get_nclist_ens(expname,vname_in,freq,scenario,realm='ocn',regrid=True,debug=True)
        nclists_byscenario.append(nclists)
    
    for e in tqdm(np.arange(nens)):#range(nens)):
        st  = time.time()
        ens = enslist[e]
        
        # Load DS by Scenario
        ss = 0
        ds_byscenario = []
        for ss in range(2):
            dsview_scenario = xr.open_mfdataset(nclists_byscenario[ss][e],concat_dim='time',combine='nested',decode_times=True)
            dsview_scenario = remove_duplicate_times(dsview_scenario,verbose=True,timename='time')
            dsslice         = dsview_scenario.sel(time=slice(tstart,tend))[vname_in]
            ds_byscenario.append(dsslice)
        
        # Merge In Time (Note that 2026-01-02 is missing), Takes ~ 10 Seconds
        dsmerge = xr.concat(ds_byscenario,dim='time')
        dsmerge = remove_duplicate_times(dsmerge)
        dsmerge = dsmerge.sortby('time')
        dsmerge = dsmerge.load()
            
        
        # ====================
        
        
        # Remove Seasonality (~9sec)
        dsanom = deseason_daily(dsmerge,period=climperiod)
        dsmerge.close()
        del dsmerge
            
        # Remove [deg]-order trend
        dsanom = detrend_dim(dsanom,dim='time',deg=deg)
            
        # Save Output
        outname      = "%smesaclip_%s_%s_%s_anom_ens%03i.nc" % (outpath_proc,expname,freq,vname_in,ens)
        dsanom.to_netcdf(outname)
        print("Processed Ens %03i in %.2fs" % (ens,time.time()-st))
        print("\tSaved to %s" % outname)
        dsanom.close()
        del dsanom
        
    
        
    

