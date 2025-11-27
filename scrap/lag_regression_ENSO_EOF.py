#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Lag Regression ENSO EOF

    Regress selected variable (pointwise) onto EP and CP indices

Created on Mon Nov 24 18:48:18 2025

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
#import climlab
import cartopy.feature as cfeature
#from tqdm import tqdm

import pandas as pd

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% Declare Experiment

# Regression Calculation Options
leadlags    = np.arange(-18,19,1)
standardize = True
sep_mon     = False

# target variable information (anomalized and detrended)
# --------------------------------------


# # .............................................................................
# # ERA5 1979 2024 Regridded
# expname     = "ERA5_1979_2024"
# vname       = "sst"
# datpath     = "/home/niu4/gliu8/projects/common_data/ERA5/regrid_1x1/anom_detrend1/"
# ncsearch    = datpath+"%s_1979_2024.nc" # nc search string, can also provide actual nc
# outpath     = "/home/niu4/gliu8/projects/ccfs/enso_regression_patterns/anom_detrend1/regrid_1x1/"
# # .............................................................................

# .............................................................................
# AWI-CM3 Regridded
expnames    = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
# Start Loop here
for expname in expnames:
    vname       = "sst"
    datpath     = "/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/"
    ncsearch    = datpath+expname+"_%s_regrid1x1.nc" # nc search string, can also provide actual nc
    outpath     = "/home/niu4/gliu8/projects/ccfs/enso_regression_patterns/anom_detrend1/regrid_1x1/"
    # .............................................................................
    
    
    
    
    # --------------------------------------
    
    # Set Output Name
    # <outpath>/<expname>_<vname>_standardize<0>_lags<%02i>-lags<%02i>.nc
    outname     = outpath+"%s_%s_standardize%i_lags%02i.nc"  % (expname,vname,standardize,leadlags[-1])
    
    # ENSO information
    ensopath    = "/home/niu4/gliu8/projects/scrap/nino34/"
    print("File will be saved to %s" % outname)
    
    #%% Load ENSO indices (rotated PCs)
    
    
    ncname      = "%s%s_enso_eof_rotated.nc" % (ensopath,expname)
    dsenso      = xr.open_dataset(ncname).load()
    
    cp          = dsenso.cp
    ep          = dsenso.ep
    X           = np.array([ep,cp])
    X_std       = np.array([ep/np.std(ep), cp/np.std(cp)])
    
    #%% Load variable
    
    st    = time.time()
    ncvar = ncsearch % vname
    ds    = xr.open_dataset(ncvar)[vname].load().squeeze()
    proc.printtime(st,"Loaded in ")
    
    ds    = ut.standardize_names(ds)
    
    #ds    = ds.transpose('time','lat','lon')
    
    # Check to make sure time matches
    ds,cp = proc.match_time_month(ds,cp)
    
    
    #%% Perform regression
    
    
    ensoids      = [cp,ep]
    ensoid_names = ["rEOF-CP","rEOF-EP"]
    dsvar        = ds # From [lag_regression_enso_global]
    
    for ei in range(2):
        
        # Select Variables
        ensoid = ensoids[ei]
        ensoid_name = ensoid_names[ei]
        
        
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
            for ll in tqdm.tqdm(range(nlags)):
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
            for ll in tqdm.tqdm(range(nlags)):
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
        outname1      = proc.addstrtoext(outname,"_%s" % ensoid_name,adjust=-1)
        if sep_mon is True:
            outname1 = proc.addstrtoext(outname1,"_SeparateMonths",adjust=-1)
        
        edict   = proc.make_encoding_dict(ds_out)
        ds_out.to_netcdf(outname1,encoding=edict)




#%%





#%%



