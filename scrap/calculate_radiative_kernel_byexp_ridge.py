#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Radiative Kernels by Experiment using Ridge Regression

Created on Wed Mar  4 18:13:42 2026

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
import matplotlib as mpl
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

landmask = ut.load_land_mask_awi("TCo319",regrid=True)

# def correct_lon(lonval,verbose=True):
#     "Correct longitude values to between 0-360 or -180-180 (for get_box)"
#     if lonval > 180:
#         if verbose:
#             print("Degrees East Detected (0 to 360)")
#         if lonval >= 360:
#             if verbose:
#                 print("\tCrossing Prime Meridian")
#             lonval = lonval - 360        
#     else:
#         if verbose:
#             print("Degrees West Detected (-180 to 180)")
#         if lonval < -180:
#             if verbose:
#                 print("Crossing Date Line E to W")
#             lonval = 360 + lonval
#         elif lonval > 180:
#             if verbose:
#                 print("Crossing Date Line W to E")
#             lonval = lonval - 360
#     return lonval

# def correct_lat(latval,verbose=True):
#     "Correct latitude values that go beyond -90 or 90 (for get_box)"
#     if latval < -90:
#         if verbose:
#             print("Warning! Box extends past the poles, cropping to -90")
#         latval = -90
#     if latval > 90:
#         if verbose:
#             print("Warning! Box extends past the poles, cropping to 90")
#         latval = 90
    
#     return latval

# def get_box(lonc,latc,lonwin,latwin,verbose=False):
#     "Get box centered at lat/lon with specificed size"
#     bbsel  = np.array([lonc-lonwin,lonc+lonwin,latc-latwin,latc+latwin])
#     for ii in [0,1]:
#         bbsel[ii] = correct_lon(bbsel[ii],verbose=verbose)
#     for ii in [2,3]:
#         bbsel[ii] = correct_lat(bbsel[ii],verbose=verbose)
#     return bbsel

# def sel_box(ds,lonc,latc,lonwin,latwin,verbose=False,return_bbsel=False):
#     bbsel = get_box(lonc,latc,lonwin,latwin,verbose=verbose)
#     if return_bbsel:
#         return proc.sel_region_xr(ds,bbsel),bbsel
#     return proc.sel_region_xr(ds,bbsel,verbose=verbose)

def ridge(X,y,alpha):
    # Ridge Regression Function fit using scipy 
    # X is # [ time x predictor ]
    # y is # [ time x 1]
    
    # Initialize Model and Fit
    model             = sklearn.linear_model.Ridge(alpha=alpha)
    model.fit(X,y)
    pred              = model.predict(X)
    # Calculate Error and other variables
    ridge_out           = {}
    ridge_out['pred']   = pred
    ridge_out['err']    = y - pred # Model Error # []
    ridge_out['coeffs'] = model.coef_ # Coefficients
    ridge_out['r2']     = sklearn.metrics.r2_score(y,pred)
    
    return ridge_out


def ridge_ccfs_box(ccfs,flx,lonw,latw,alpha=1,standardize=True,verbose=False,
                   fill_value=0,std_mon=False):
    #    ccfs : LIST of DataArrays [ccf][time x lat x lon]
    #    flx  :  DataArray [time x lat x lon]
    #    lonw : Degrees of Longitude to consider for box
    #    latw : Degress of Latitude to consider for box
    # Skips NaN points based on first CCF...
    # Lat/Lon is Referenced to the flx
    st              = time.time()
    
    ccfs            = [ds.transpose('time','lat','lon') for ds in ccfs]
    flx             = flx.transpose('time','lat','lon')
    
    ntime,nlat,nlon = flx.shape # Get Shape
    nccfs           = len(ccfs)
    
    # Determine Lat/Lon Box
    lon             = flx.lon
    lat             = flx.lat
    dx              = int(np.ceil(lon.data[1:] - lon.data[:(-1)]).mean().item())
    dy              = int(np.ceil(lat.data[1:] - lat.data[:(-1)]).mean().item())
    nlatbox         = int(np.ceil(latw/dy*2+1))
    nlonbox         = int(np.ceil(lonw/dx*2+1))

    
    boxcoords_lon   = np.zeros((nlat,nlon,nlonbox))
    boxcoords_lat   = np.zeros((nlat,nlon,nlatbox))
    
    #box_coords      = np.zeros((nlat,nlon,nlatbox,nlonbox)) * np.nan # Matrix of boxes
    coeffs          = np.zeros((nccfs,nlat,nlon,nlatbox,nlonbox))    # MLR fit values
    r2              = np.zeros((nlat,nlon))       # r2 Fit
    ypred           = np.zeros((ntime,nlat,nlon)) # Predicted values
    yerr            = np.zeros((ntime,nlat,nlon)) # Error
    
    
    if std_mon:
        def std_mon(ds):
            ds = ds - ds.mean() # Center
            return ds.groupby('time.month') / ds.groupby('time.month').std('time') # Standardize
    else:
        def std_mon(ds): 
            ds = ds - ds.mean() # Center
            return ds / ds.std() # Standardize
    
    for o in tqdm.tqdm(range(nlon)):
        for a in range(nlat):
            lonf = lon[o]
            latf = lat[a]
            
            ccf_boxes       = [proc.sel_box(ds,lonf,latf,lonw,latw,verbose=False) for ds in ccfs]
            ccf_boxes_std   = [std_mon(ds) for ds in ccf_boxes] # Standardize along time dimension
            _,nlatbox_loop,nlonbox_loop=ccf_boxes[0].shape
            
            # Skip NaN points based on first CCF
            if np.any(np.isnan(proc.selpt_ds(ccfs[0],lonf,latf))):
                continue
            
            # Reshape to time x space, and concatenate along space dimension
            #_,nlatbox,nlonbox = ccf_boxes.shape
            ccf_predictors  = np.concatenate([ds.data.reshape(ntime,nlatbox_loop*nlonbox_loop) for ds in ccf_boxes_std],axis=1)
            X               = ccf_predictors.data # [time x predictor]
            y               = proc.selpt_ds(flx,lonf,latf).data#[:,None]
            
            # Replace Zeros
            X               = np.where(np.isnan(X),fill_value,X)
            y               = np.where(np.isnan(y),fill_value,y)
            
            # Perform Ridge Regression Fit
            ridge_out       = ridge(X,y,alpha)
            
            # Adjust lat box index
            latboxidx      = np.arange(nlatbox_loop)
            lonboxidx      = np.arange(nlonbox_loop)
            
            
            # Read variables into arrays
            r2[a,o]        = ridge_out['r2']
            yerr[:,a,o]    = ridge_out['err']
            ypred[:,a,o]   = ridge_out['pred']
            coeffspt       = ridge_out['coeffs']
            coeffspt       = coeffspt.reshape(nccfs,nlatbox_loop,nlonbox_loop)
            #coeffs[:,a,o,latboxidx,lonboxidx] = coeffspt.copy()
            coeffs[:,a,o,:nlatbox_loop,:nlonbox_loop] = coeffspt.copy()
            boxcoords_lat[a,o,:nlatbox_loop] = ccf_boxes[0].lat.data
            boxcoords_lon[a,o,:nlonbox_loop] = ccf_boxes[0].lon.data
    
    
    coords_coeff = dict(ccf=ccf_vars,lat=lat,lon=lon,latbox=np.arange(nlatbox),lonbox=np.arange(nlonbox))
    da_coeffs = xr.DataArray(coeffs,coords=coords_coeff,dims=coords_coeff,name="coeffs")
    
    
    coords_r2       = dict(lat=lat,lon=lon)
    coords_pred     = dict(time=flx.time,lat=lat,lon=lon)
    coords_lonbox   = dict(lat=lat,lon=lon,lonbox=np.arange(nlonbox))
    coords_latbox   = dict(lat=lat,lon=lon,latbox=np.arange(nlatbox))
    
    da_r2           = xr.DataArray(r2,coords=coords_r2,dims=coords_r2,name='r2')
    da_pred         = xr.DataArray(ypred,coords=coords_pred,dims=coords_pred,name='ypred')
    da_yerr         = xr.DataArray(yerr,coords=coords_pred,dims=coords_pred,name='yerr')
    
    da_lonbox       = xr.DataArray(boxcoords_lon,coords=coords_lonbox,dims=coords_lonbox,
                                   name="boxlon")
    
    da_latbox       = xr.DataArray(boxcoords_lat,coords=coords_latbox,dims=coords_latbox,
                                   name="boxlat")
    
    ds_out          = xr.merge([da_r2,da_coeffs,da_pred,da_yerr,da_lonbox,da_latbox])
    
    
    return ds_out

# =================
#%% User Selections
# =================

regrid_ver       = "regrid_1x1"

# # Redo CERES-EBAF, limit to 2001-2024
# expname         = "CERES_EBAF_ERA5_2001_2024"
# datpath         = "/home/niu4/gliu8/projects/ccfs/input_data/%s/%s/anom_detrend1/" % (regrid_ver,expname)#,#/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/"#"/home/niu4/gliu8/projects/scrap/regrid_1x1/"
# outpath         = "/home/niu4/gliu8/projects/ccfs/kernels/ridge_kernels/%s/" % regrid_ver#%s/" % expname #/home/niu4/gliu8/projects/ccfs/kernels/regrid_1x1/"
# anomalize       = False # Kept for legacy. Input should be anomalized before using `anom_detrend1' shellscripts
# tstart          = '2001-01-01'
# tend            = '2024-12-31'
# customname      = None

# Do for TCo319-1950-gibbs_charn
expname         = "TCo319-DART-ctl1950d-gibbs-charn"
datpath         = "/home/niu4/gliu8/projects/ccfs/input_data/%s/%s/anom_detrend1/" % (regrid_ver,expname) #" % expname
outpath         = "/home/niu4/gliu8/projects/ccfs/kernels/ridge_kernels/%s/" % regrid_ver #" # Expname Added later
anomalize       = False # Kept for legacy. Input should be anomalized before using `anom_detrend1' shellscripts
tstart          = None
tend            = None
customname      = None


# Variables
flxname         = 'cre' #['allsky','clearsky','cre']  # Loop for fluxes
ccf_vars        = ["sst","eis","Tadv","r700","w700","ws10",]#"ucc"] 
#ccf_vars        = ["sst","eis","MeanAdvTanom","AnomAdvTmean","r700","w700","ws10",]

selmons_loop    = [None,[12,1,2],[3,4,5],[6,7,8],[9,10,11]] #[[12,1,2],[3,4,5],[6,7,8],[9,10,11]] # [None,]# # Set to None to do 

# MLR Calculation Options
standardize     = True # Set to True to standardize predictors before MLR
fill_value      = 0    # Replace NaN values with <fill_value>
#add_ucc         = False # Set to True to include upper cloud concentration as a predictor

#%% Now Load each predictor variable to do CCFs calculation
# Note that some were computed in [calc_ccfs_regrid.py]
# Others were preprocessed using remapbil in cdo

"""
Searches for datasets in (or 5x5)
    .../ccfs/input_data/regrid_1x1/<expname>/anom_detrend1/

Outputs computed kernels in
    .../ccfs/kernels/regrid_1x1/<expname>/
"""

# Prepare Output Path
outpath_kernel = outpath + "%s/" % expname
proc.makedir(outpath_kernel)

# Load the CCFs
nccfs = len(ccf_vars)
dsvars_anom = []
for v in range(nccfs):
    
    ncname = "%s%s.nc" % (datpath,ccf_vars[v])
    ds     = xr.open_dataset(ncname)[ccf_vars[v]].load()
    
    # Standardize Names
    ds     = ut.standardize_names(ds)
    
    # <Need to add additional fixes/automatic checks>
    if expname == 'TCo2559-DART-1950C':
        if ccf_vars[v] == "w700": # Duplicate a month for the missing variables
            w700data = ds.data.squeeze()
            w700data_duplicate_jan1950 = np.concatenate([w700data[[0],...],w700data],axis=0)
            newtime = dsvars_anom[v-1].time
            coords  = dict(time=newtime,lat=ds.lat,lon=ds.lon)
            w700new = xr.DataArray(w700data_duplicate_jan1950,coords=coords,dims=coords,name='w700')
            ds = w700new
        
        print("%s, %s" % (expname,ccf_vars[v]))
        print(ds.shape)
        print("")
        dsvars_anom.append(ds.squeeze())
    else:
        dsvars_anom.append(ds.squeeze())
    
    # SST is only 8 years, so reduce the time...
    # Commented out this since it has changed
    #if expname == 'TCo2559-DART-1950C': 
        #dsvars_anom = [ut.reduce_time(ds,dsvars_anom[0]) for ds in dsvars_anom]

# Load the fluxes
dsflx   = xr.open_dataset("%s%s.nc" % (datpath,flxname))[flxname].load()
if "TCo" in expname: # Convert Fluxes
    print("Converting fluxes based on accumulation period for AWI-CM3")
    dsflx   = ut.varcheck(dsflx,flxname,expname)
elif expname == "ERA5_1979_2024":
    print("Considering 1-day flux accumulation period for ERA_1979_2024")
    dtday   = 3600*24
    dsflx   = dsflx / dtday
dsflx   = ut.standardize_names(dsflx)
dsflx   = ut.remove_duplicate_times(dsflx)

# Shift CERES data (move function to ut)
import pandas as pd
def shift_time_monthstart(dsin,timename='time'):
    oldtime         = dsin[timename]
    tstart          =  str(oldtime[0].data)[:7] + "-01"
    tend            =  str(oldtime[-1].data)[:7] + "-01"
    newtime         = pd.date_range(start=tstart,end=tend,freq="MS")
    dsin[timename]    = newtime
    print("New time dimension between %s and %s" % (dsin[timename][0].data,dsin[timename][-1].data))
    return dsin
if "CERES" in expname:
    dsflx = ut.shift_time_monthstart(dsflx,timename='time')
    
# Limit to time
if tstart or tend is not None:
    dsflx       = dsflx.sel(time=slice(tstart,tend))
    dsvars_anom = [ds.sel(time=slice(tstart,tend)) for ds in dsvars_anom]
else:
    print("No time range will be applied")
    
# for ii in range(nccfs): (note match time month does not play nice with pd.date_range, need to figure it out)
#     dsvars_anom[ii],_=proc.match_time_month(dsvars_anom[ii],dsflx)

# Check Time Dimension
ntimes_predictors = [len(ds.time) for ds in dsvars_anom]
ntimes_flux       = len(dsflx.time)


#%% Now Try to Do Ridge Regression

selmons = None
dsexp_sel  = dsvars_anom
dsexp_flx  = dsflx

if selmons is not None:
    dsexp_flx = proc.selmon_ds(dsexp_flx,selmons)
    dsexp_sel = [proc.selmon_ds(ds,selmons) for ds in dsexp_sel]
else:
    print("Calculating for all months!")

ccfs      = dsexp_sel
flx       = dsexp_flx
lonw      = 10*5#20
latw      = 5*5#20   
alpha     = 1e3

ds_out          = ridge_ccfs_box(ccfs,flx,lonw,latw,alpha=alpha,std_mon=False)
edict           = proc.make_encoding_dict(ds_out)
outpath_exp     = "%s%s/" % (outpath,expname)
proc.makedir(outpath_exp)
outname         = "%s%s_ridge_kernels_lonw%02i_latw%02i_alpha%.0e.nc" % (outpath_exp,flxname,lonw,latw,alpha)
ds_out.to_netcdf(outname)
            
            


 