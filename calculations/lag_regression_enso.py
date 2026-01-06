#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Perform Lag Regression for ENSO and selected variable

Created on Thu Oct  9 16:01:18 2025

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

#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% User Edits

# Process variables for BWJ Analysis
# expnames        = ["TCo319_ctl1950d","TCo319_ssp585",]#"glorys",]#["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
# expnames_long   = ["31km Control","31km SSP585"]#"GLORYS Reanalysis"]#["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"]
# ensoid_name     = "nino34"
# standardize     = False
# timecrops        = [[1950,2100],None]#[1993,2024]]
# datpath         = "/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/"#"/home/niu4/gliu8/projects/scrap/TP_crop/anom_detrend2/"
# outpath         = "/home/niu4/gliu8/projects/enso_feedbacks/lag_regressions/regrid_1x1/" #"/home/niu4/gliu8/projects/scrap/TP_crop/lag_regressions/"
# ncname_rep      = "%s%s_%s_regrid1x1.nc" #"%s%s_%s_anom.nc" (datpath,expnames[ex],vname)
# vnames          = ["temp75","vvel_ML50","uvel_ML50",'sst',"tx_sur","ssh"] #str','ssr','skt','ssh','lcc','tcc','ttr','ttrc','tsr','tsrc'] # 'sst',#

# Calculations with updated AWi-Cm3 output
expnames        = ["TCo1279-DART-2060","TCo319-DART-ctl1950d-gibbs-charn","TCo319-DART-ssp585d-gibbs-charn",]
expnames_long   = ["9km 2060","31km Control","31km SSP5.85"]#"GLORYS Reanalysis"]#["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"]
ensoid_name     = "nino34"
standardize     = False
timecrops       = [None,None,None]#[1993,2024]]
datpath         = "/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1/"#"/home/niu4/gliu8/projects/scrap/TP_crop/anom_detrend2/"
outpath         = "/home/niu4/gliu8/projects/enso_feedbacks/lag_regressions/" #"/home/niu4/gliu8/projects/scrap/TP_crop/lag_regressions/"
ncname_rep      = "%s%s_%s.nc" #"%s%s_%s_anom.nc" (datpath,expnames[ex],vname)
vnames          = ["sst" "msl"] #str','ssr','skt','ssh','lcc','tcc','ttr','ttrc','tsr','tsrc'] # 'sst',#


figpath         = "/home/niu4/gliu8/figures/bydate/2026-01-13/"
proc.makedir(figpath)


#%% Load ENSO ID

ensoids         = [ut.load_ensoid(expname,ensoid_name,standardize=standardize) for expname in expnames]

#%% Looping for each experiment

# Indicate Lead Lags
leadlags    = np.arange(-24,25,1)
#sep_mon     = True

#ex          = 0
nexps       = len(expnames)

for vname in vnames:
    for sep_mon in [False,]:#True]:
        for ex in tqdm.tqdm(range(nexps)):
            
            # Load the variable
            st              = time.time()
            #ncname          = #"%s%s_%s_anom.nc" % (datpath,expnames[ex],vname)
            ncname           = ncname_rep % (datpath,expnames[ex],vname)
            try:
                dsvar           = xr.open_dataset(ncname).load()[vname]
            except:
                print("Could not find %s for %s... skipping." % (vname,expnames_long[ex]))
                continue
            proc.printtime(st,print_str="Loaded")
            
            # Apply name check
            dsvar = ut.standardize_names(dsvar)
            
            # Crop time (mostly for control run, pre 1950)
            timecrop = timecrops[ex]
            if timecrop is not None:
                print("Cropping time for %s: %s to %s" % (expnames_long[ex],str(timecrop[0])+'-01-01',str(timecrop[1])+'-12-31'))
                dsvar = dsvar.sel(time=slice(str(timecrop[0])+'-01-01',str(timecrop[1])+'-12-31'))
            
            # Check to make sure the time matches
            ensoid          = ensoids[ex]
            dsvar,ensoid    = proc.match_time_month(dsvar,ensoid)
            
            # Get Dimension Lengths
            dsvar           = dsvar.squeeze().transpose('lon','lat','time')
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
                    lag                = leads[ll]
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
            ds_out   = xr.merge(da_out)
            
            #% Save the output
            if sep_mon is False:
                outname = "%sLagRegression_AllMonths_%s_%s_%s_standardize%i_lag%ito%i.nc" % (outpath,expnames[ex],vname,ensoid_name,standardize,leadlags[0],leadlags[-1])
            else:
                outname = "%sLagRegression_SeparateMonths_%s_%s_%s_standardize%i_lag%ito%i.nc" % (outpath,expnames[ex],vname,ensoid_name,standardize,leadlags[0],leadlags[-1])
            
            edict   = proc.make_encoding_dict(ds_out)
            ds_out.to_netcdf(outname,encoding=edict)



