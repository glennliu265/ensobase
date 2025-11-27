#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Streamlined Version of calc_ccfs_regrid

Created on Fri Nov 14 10:54:27 2025

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
import importlib

from sklearn.linear_model import LinearRegression
import sklearn

#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% Other Functions


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

def reduce_time(ds,dsst):
    dsnew,dsst = proc.match_time_month(ds,dsst)
    return dsnew

#%% Load Land Mask


landmask = ut.load_land_mask_awi("TCo319",regrid=True)

# =================
#%% User Selections
# =================

# Path to Data and Experiments
expnames    = ["TCo1279-DART-1950","TCo2559-DART-1950C"]
datpath     = "/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/"#"/home/niu4/gliu8/projects/scrap/regrid_1x1/"
outpath     = "/home/niu4/gliu8/projects/ccfs/kernels/regrid_1x1/"

anomalize   = False # Kept for legacy. Input should be anomalized before using `anom_detrend1' shellscripts

# Variables
flxnames    = ['cre']#['allsky','clearsky','cre']  # Loop for fluxes
ccf_vars    = ["sst","eis","Tadv","r700","w700","ws10",]#"ucc"] 

selmons_loop = [[12,1,2],[3,4,5],[6,7,8],[9,10,11]] # Set to None to do 

# MLR Calculation Options
standardize = True # Set to True to standardize predictors before MLR
fill_value  = 0    # Replace NaN values with <fill_value>
add_ucc     = False # Set to True to include upper cloud concentration as a predictor

#%% Now Load each predictor variable to do CCFs calculation
# Note that some were computed in [calc_ccfs_regrid.py]
# Others were preprocessed using remapbil in cdo

dsbyexp = []
for ex in range(2):
    dsvars = []
    for v in range(len(ccf_vars)):
        
        ncname = "%s%s_%s_regrid1x1.nc" % (datpath,expnames[ex],ccf_vars[v])
        ds     = xr.open_dataset(ncname)[ccf_vars[v]].load()
        ds     = ut.standardize_names(ds)
        
        if expnames[ex] == 'TCo2559-DART-1950C':
            if ccf_vars[v] == "w700": # Duplicate a month for the missing variables
                w700data = ds.data.squeeze()
                w700data_duplicate_jan1950 = np.concatenate([w700data[[0],...],w700data],axis=0)
                newtime = dsvars[v-1].time
                coords  = dict(time=newtime,lat=ds.lat,lon=ds.lon)
                w700new = xr.DataArray(w700data_duplicate_jan1950,coords=coords,dims=coords,name='w700')
                ds = w700new
        
        print("%s, %s" % (expnames[ex],ccf_vars[v]))
        print(ds.shape)
        print("")
        dsvars.append(ds.squeeze())
    
    # SST is only 8 years, so reduce the time...
    if expnames[ex] == 'TCo2559-DART-1950C': 
        dsvars = [reduce_time(ds,dsvars[0]) for ds in dsvars]
    
    dsbyexp.append(dsvars)

#%% Anomalize and detrend variables

if anomalize:
    dsbyexp_anoms = []
    for ex in range(2):
        dsvars_anoms = []
        
        for v in tqdm.tqdm(range(len(ccf_vars))):
            
            dsin    = dsbyexp[ex][v]
            dsanoms = ut.preprocess_enso(dsin)
            dsvars_anoms.append(dsanoms)
        
        dsbyexp_anoms.append(dsvars_anoms)
else:
    dsbyexp_anoms = dsbyexp
    



# =====================================
#%% Part (2): Compute Radiative Kernels
# =====================================


for flxname in flxnames:
    
    # Load Fluxes for each experiment
    dsflxs  = []
    for ex in range(2):
        dsflx   = xr.open_dataset("%s%s_%s_regrid1x1.nc" % (datpath,expnames[ex],flxname)).load()
        if anomalize:
            dsflx   = ut.preprocess_enso(ut.standardize_names(dsflx[flxname]))
        dsflx   = ut.varcheck(dsflx,flxname,expnames[ex])
        dsflx   = ut.standardize_names(dsflx)
        dsflxs.append(dsflx)
    
    # Here might be the point to start selmons loop...
    for selmons in selmons_loop:
        
        
        # Perform MLR for each experiment
        for ex in range(2):
            st = time.time()
            
            dsexp_sel      = dsbyexp_anoms[ex]
            if add_ucc is False and "ucc" in ccf_vars:
                print("Reducing dimensions by 1")
                dsexp_sel = dsexp_sel[:-1]
            dsexp_flx      = dsflxs[ex]
            
            # Check time dimension
            ntimes_predictors = [len(ds.time) for ds in dsexp_sel]
            ntimes_flux       = len(dsexp_flx.time)
            if expnames[ex] == 'TCo2559-DART-1950C':
                print("Adjusting flux length")
                dsexp_flx,_ = proc.match_time_month(dsexp_flx,dsexp_sel[0])
            
            # Subset months
            if selmons is not None:
                dsexp_flx = proc.selmon_ds(dsexp_flx,selmons)
                dsexp_sel = [proc.selmon_ds(ds,selmons) for ds in dsexp_sel]
            else:
                print("Calculating for all months!")
            
            # Pre-allocate
            dsexp_flx       = dsexp_flx.transpose('lat','lon','time')
            lon             = dsexp_flx.lon.data
            lat             = dsexp_flx.lat.data
            nlat,nlon,ntime = dsexp_flx.shape
            nccfs           = len(dsexp_sel)
            coeffs          = np.zeros((nlat,nlon,nccfs)) * np.nan # [ Lat x Lon x CCFs ]
            ypred           = np.zeros(dsexp_flx.shape) * np.nan   # [ Lat x Lon x Time ]
            r2              = np.zeros((nlat,nlon)) * np.nan       # [ Lat x Lon ]
            
            # Do a silly loop (took 5 min 17 sec)
            for o in tqdm.tqdm(range(nlon)):
                lonf = lon[o]
                
                for a in range(nlat):
                    latf = lat[a]
                    
                    chkland = proc.selpt_ds(landmask,lonf,latf).data
                    if np.isnan(chkland):
                        continue
                    
                    # Check for NaN in predictor
                    dspts  = [proc.selpt_ds(ds,lonf,latf) for ds in dsexp_sel]
                    chknan = [np.any(np.isnan(ds.data)) for ds in dspts]
                    if np.any(chknan):
                        iinan = np.where(chknan)[0][0]
                        #print("NaN detected for variables %s, lon (%.2f), lat (%.2f)... skipping." % (chknan,lonf,latf))
                        continue
                    # Check for NaN in target
                    flxpt  = proc.selpt_ds(dsexp_flx,lonf,latf)
                    if np.any(np.isnan(flxpt.data)):
                        #print("NaN detected for Flux, lon (%.2f), lat (%.2f)... skipping." % (lonf,latf))
                        continue
                    
                    # Do calculations
                    mlr_out = ut.mlr_ccfs(dspts,flxpt,standardize=standardize,verbose=False)
                    
                    r2[a,o] = mlr_out['r2']
                    ypred[a,o,:] = mlr_out['pred']
                    coeffs[a,o,:] = mlr_out['coeffs']
            
            if add_ucc:
                ccfnames = ccf_vars
            else:
                if "ucc" in ccf_vars:
                    ccfnames = ccf_vars[:-1]
                ccfnames = ccf_vars
            
            coords_r2       = dict(lat=lat,lon=lon)
            coords_coeffs   = dict(lat=lat,lon=lon,ccf=ccfnames)
            coords_pred     = dict(lat=lat,lon=lon,time=dsexp_flx.time)
            
            da_r2           = xr.DataArray(r2,coords=coords_r2,dims=coords_r2,name='r2')
            da_coeffs       = xr.DataArray(coeffs,coords=coords_coeffs,dims=coords_coeffs,name='coeffs')
            da_pred         = xr.DataArray(ypred,coords=coords_pred,dims=coords_pred,name='ypred')
            ds_out          = xr.merge([da_r2,da_coeffs,da_pred])
            edict           = proc.make_encoding_dict(ds_out)
            outname         = "%s%s_%s_CCFs_Regression_standardize%i_adducc%i.nc" % (outpath,expnames[ex],flxname,standardize,add_ucc)
            if selmons is not None:
                selmonstr = proc.mon2str(np.array(selmons)-1)
                outname      = proc.addstrtoext(outname,"_"+selmonstr,adjust=-1)
            
            ds_out.to_netcdf(outname,encoding=edict)
            
            print("Completed CCF kernel calculation for %s (%s) in %.2fs" % (flxname,expnames[ex],time.time()-st))

    