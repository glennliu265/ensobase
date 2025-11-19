#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Copied calculate_radiative_kernels, but generalize... to any simulation

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
    y = flxpt.data
    
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

# =================
#%% User Selections
# =================

regrid1x1   = False
expname     = "ERA5_1979_2024"
datpath     = "/home/niu4/gliu8/projects/common_data/ERA5/anom_detrend1/"
outpath     = "/home/niu4/gliu8/projects/ccfs/"

if regrid1x1:
    datpath = "/home/niu4/gliu8/projects/common_data/ERA5/regrid_1x1/anom_detrend1/"
    outpath = "/home/niu4/gliu8/projects/ccfs/regrid_1x1/"

flxnames     = ['cre',]#'allsky','clearsky']# ['cre',]
ccf_vars     = ["sst","eis","Tadv","r700","w700","ws10",]#"ucc"] 
ncstr        = datpath + "%s_1979_2024.nc"  
selmons_loop = [None,]#[[12,1,2],[3,4,5],[6,7,8],[9,10,11]] # Set to None to do 


tstart   = '1979-01-01'
tend     = '2024-12-31'
timename = 'valid_time'
latname  = 'latitude'
lonname  = 'longitude'

#% Load Land Mask
landnc   = datpath + "mask_1979_2024.nc"
landmask = xr.open_dataset(landnc).mask



# MLR Calculation Options
standardize = True # Set to True to standardize predictors before MLR
fill_value  = 0    # Replace NaN values with <fill_value>
add_ucc     = False # Set to True to include upper cloud concentration as a predictor

#%% First, check preprocessed CCFs and fluxes

nvars  = len(ccf_vars)
nclist = []
for v in range(nvars):
    vname = ccf_vars[v]
    ncsearch = ncstr % (vname)
    foundnc = glob.glob(ncsearch)
    print("Found the following for %s:" % vname)
    print("\t"+ str(foundnc))
    print("\n")
    try:
        nclist.append(foundnc[0])
    except:
        nclist.append(None)
        print("Nothing found for %s" % vname)
print("First element of each will be taken")

#%% Load the predictors

# Load Predictors
dsvars = []
for v in tqdm.tqdm(range(nvars)):
    vname = ccf_vars[v]
    nc    = nclist[v]
    ds    = xr.open_dataset(nc)[vname].load()
    dsvars.append(ds)

# Manual Fix (need to standardize approach)
# Fixed dimensions wtih  [manual_replace_eis_coords]
#dsvars[1] = dsvars[1].rename(dict(time='valid_time'))

# Make sure they all have the right shape, dim names
def preprocess(ds,tstart,tend,timename,latname,lonname):
    print("\n now preprocessing")
    rename_dict = {
        timename : 'time',
        latname : 'lat',
        lonname : 'lon'
        }
    if 'time' not in ds.coords or timename not in ds.coords:
        del rename_dict[timename]
    if 'lat' in ds.coords:
        print("lat already found...")
        del rename_dict[latname]
        #rename_dict = rename_dict
    if 'lon' in ds.coords:
        print("lon already found...")
        del rename_dict[lonname]
    
    ds = ds.rename(rename_dict)
    ds = ut.standardize_names(ds)
    
    ds = ds.squeeze()
    
    
    if 'time' in ds.coords and len(ds.shape) >= 3:
        ds = ds.sel(time=slice(tstart,tend))
    return ds
dsvars_anoms = [preprocess(ds,tstart,tend,timename,latname,lonname) for ds in dsvars]

landmask = preprocess(landmask,tstart,tend,timename,latname,lonname)

#%%

# Looping for fluxes
#ff      = 0
for ff in range(len(flxnames)):
    flxname = flxnames[ff]
    
    # Load the flux
    ncsearch    = ncstr % (flxname)    
    foundnc     = glob.glob(ncsearch)
    print("Found the following for %s:" % flxname)
    print("\t"+ str(foundnc))
    try:
        dsflx = xr.open_dataset(foundnc[0])[flxname]
    except:
        print("No flux found... (%s)" % flxname)
    dsflx_anom = preprocess(dsflx,tstart,tend,timename,latname,lonname)
    
    #% ----------------------------------------------------------------
    
    # Here might be the point to start selmons loop...
    st = time.time()
    for selmons in selmons_loop:
        
        # Subset months
        if selmons is not None:
            dsin_flx  = proc.selmon_ds(dsflx_anom,selmons)
            dsin_vars = [proc.selmon_ds(ds,selmons) for ds in dsvars_anoms]
        else:
            dsin_flx  = dsflx_anom
            dsin_vars = dsvars_anoms
            print("Calculating for all months!")
            
        
        # Pre-allocate
        lon             = dsin_flx.lon.data
        lat             = dsin_flx.lat.data
        dsin_flx        = dsin_flx.transpose('lat','lon','time')
        nlat,nlon,ntime = dsin_flx.shape
        nccfs           = len(dsin_vars)
        coeffs          = np.zeros((nlat,nlon,nccfs)) * np.nan # [ Lat x Lon x CCFs ]
        ypred           = np.zeros(dsin_flx.shape) * np.nan   # [ Lat x Lon x Time ]
        r2              = np.zeros((nlat,nlon)) * np.nan       # [ Lat x Lon ]
        
        # Do a silly loop (took 5 min 17 sec)
        for o in tqdm.tqdm(range(nlon)):
            lonf = lon[o]
            
            for a in range(nlat):
                latf = lat[a]
                
                chkland  = proc.selpt_ds(landmask,lonf,latf).data
                if np.isnan(chkland):
                    continue
                
                
                # Check for NaN in predictor
                dspts  = [proc.selpt_ds(ds,lonf,latf) for ds in dsin_vars]
                chknan = [np.any(np.isnan(ds.data)) for ds in dspts]
                if np.any(chknan):
                    iinan = np.where(chknan)[0][0]
                    #print("NaN detected for variables %s, lon (%.2f), lat (%.2f)... skipping." % (chknan,lonf,latf))
                    continue
                # Check for NaN in target
                flxpt  = proc.selpt_ds(dsin_flx,lonf,latf)
                if np.any(np.isnan(flxpt.data)):
                    #print("NaN detected for Flux, lon (%.2f), lat (%.2f)... skipping." % (lonf,latf))
                    continue
                
                # Do calculations
                mlr_out = mlr_ccfs(dspts,flxpt,standardize=standardize,verbose=False)
                
                r2[a,o] = mlr_out['r2']
                ypred[a,o,:] = mlr_out['pred']
                coeffs[a,o,:] = mlr_out['coeffs']
            
    
        
        if add_ucc:
            ccfnames = ccf_vars
        elif "ucc" in ccf_vars:
            if ypred.shape[-1] != len(ccf_vars): # Check again please...
                ccfnames = ccf_vars[:-1]
            else:
                ccfnames = ccf_vars
        else:
            ccfnames = ccf_vars    
            
        
        
        coords_r2       = dict(lat=lat,lon=lon)
        coords_coeffs   = dict(lat=lat,lon=lon,ccf=ccfnames)
        coords_pred     = dict(lat=lat,lon=lon,time=dsflx_anom.time)
        
        da_r2           = xr.DataArray(r2,coords=coords_r2,dims=coords_r2,name='r2')
        da_coeffs       = xr.DataArray(coeffs,coords=coords_coeffs,dims=coords_coeffs,name='coeffs')
        da_pred         = xr.DataArray(ypred,coords=coords_pred,dims=coords_pred,name='ypred')
        ds_out          = xr.merge([da_r2,da_coeffs,da_pred])
        edict           = proc.make_encoding_dict(ds_out)
        #outname = "test123.nc"
        outname         = "%s%s_%s_CCFs_Regression_standardize%i_adducc%i.nc" % (outpath,expname,flxname,standardize,add_ucc)
        
        if selmons is not None:
            selmonstr = proc.mon2str(np.array(selmons)-1)
            outname      = proc.addstrtoext(outname,"_"+selmonstr,adjust=-1)
        
        if regrid1x1:
            outname      = proc.addstrtoext(outname,"_regrid1x1",adjust=-1)
            
        
        ds_out.to_netcdf(outname)
        
        print("Completed CCF kernel calculation for %s (%s) in %.2fs" % (flxname,expname,time.time()-st))
    
    
    

    
    
    
    
    
