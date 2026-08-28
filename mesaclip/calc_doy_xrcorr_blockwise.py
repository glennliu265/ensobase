#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Daily Autocorrelation using XRCorr

Created on Thu Aug 27 09:46:58 2026

@author: gliu
"""

# import sys
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
from datetime import datetime
import cftime

from tqdm.notebook import tqdm

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut


#%% Additional Functions


def flatten_ens_time(ds,vname='mergebyens'):
    # Concatenate each ensemble member along the time dimension
    # and assign a new time axes
    # Reshape to Specified Dimensions
    print(ds)
    ds         = ds.transpose('ens','time','lat','lon')
    nens,ntime,nlat,nlon = ds.shape
    
    # Calculate New Period Length and Make Time Axis
    nperiods  = len(ds.time) * len(ds.ens)
    print(nperiods)
    startdate = ds.time[0].data.item()
    newtime   = xr.date_range(start=startdate,periods=nperiods,freq="D",calendar='noleap')
    
    # Make into new DataArray
    coords_new = dict(time=newtime,lat=ds.lat,lon=ds.lon)
    dsrs       = xr.DataArray(ds.data.reshape(nperiods,nlat,nlon),coords=coords_new,dims=coords_new)

    # Also include starting indices of each ensemble member
    ens_istart = np.cumsum([len(ds.isel(ens=e).time) for e in range(nens)])
    ens_dict   = dict(ens=np.arange(nens)+1)
    ens_istart = xr.DataArray(ens_istart,coords=ens_dict,dims=ens_dict)
    
    dsrs       = xr.merge([dsrs.rename(vname),ens_istart.rename('ens_start')])
    return dsrs

def fix_mesaclip_lr_date(ds):
    # Create Dummy File with zeros for timeslop
    ds_dummy         = xr.zeros_like(ds.isel(time=0))
    ds_dummy['time'] = cftime.DatetimeNoLeap(2006, 1, 2, 0, 0, 0, 0, has_year_zero=True)
    ds_dummy         = xr.where(np.isnan(ds.isel(time=1)),np.nan,0)
    return xr.combine_by_coords([ds,ds_dummy])
    #return xr.concat([dsens,ds_dummy],dim='time').sortby('time')
        
    
# def convert_size(size_bytes):
#    # https://stackoverflow.com/questions/5194057/better-way-to-convert-file-sizes-in-python
#    if size_bytes == 0:
#        return "0B"
#    size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
#    i = int(math.floor(math.log(size_bytes, 1024)))
#    p = math.pow(1024, i)
#    s = round(size_bytes / p, 2)
#    return "%s %s" % (s, size_name[i])


#%% Calculation Settings

nlags         = 180
doy_istart    = np.array(proc.get_doy_monstart())
doy_monmid    = doy_istart + 14
doys_sel      = np.sort(np.hstack([doy_istart,doy_monmid]))
winsize       = 0

ndoys = len(doys_sel)
#%% Load View of Dataset

st            = time.time()

datpath       = "/home/niu4/gliu8/share/CESM1/MESACLIP/processed/anom_detrend2_19820101-20251231/"
expname       = "mesaclip_lores"
vname         = "SST"

xrname        = "__xarray_dataarray_variable__"
ncsearch      = "%s%s_day_1_%s_anom_ens*.nc" % (datpath,expname,vname,)
nclist,nens   = proc.get_nclist(ncsearch,verbose=True)

dsall = xr.open_mfdataset(ncsearch,concat_dim='ens',combine='nested')[xrname]
print(dsall)

print("View opened in %.2fs" % (time.time()-st))

outdir_base  = "/home/niu4/gliu8/projects/mesaclip/memory/%s/" % expname
proc.makedir(outdir_base)


#%% Set Up Blocksize and output Directory

bsz       = 20

# Set Both to same for simplicity
latsize   = bsz 
lonsize   = bsz
lonblocks = np.arange(0,361,lonsize).astype(int)
latblocks = np.arange(-90,91,latsize).astype(int)

nx        = len(lonblocks)
ny        = len(latblocks)
ntotal    = nx*ny

ii        = -5
bbsel     = [lonblocks[ii],lonblocks[ii+1],latblocks[ii],latblocks[ii+1]]
dsreg     = proc.sel_region_xr(dsall,bbsel)

print("Total Blocks: %i" % ntotal )
print("Sample Arr: %s" % (dsreg))

#%% Set Output Files

bsz          = 20

outdir_temp  ="%stemp_blocksize%02i_nblocks%i/" % (outdir_base,bsz,ntotal)
proc.makedir(outdir_temp)

#%% Now Start the loop

nb      = 0
dsblocs = []
for xx in tqdm(range(nx-1)):
    for yy in range(ny-1):
        
        # Set Up Functions
        bbsel     = [lonblocks[xx],lonblocks[xx+1],latblocks[yy],latblocks[yy+1]] 
        def preproc(ds):
            dsreg = proc.sel_region_xr(ds,bbsel)
            dsreg = dsreg[xrname]
            return dsreg
        
        # Load
        # Look across all ensemble members and Concat  ~ 2min 53s
        dsens = []
        for e in tqdm(range(nens)):
        
            # Load and restrict to Region
            ds = xr.open_dataset(nclist[e])
            ds = preproc(ds)    
            dsens.append(ds.load())
        dsens = xr.concat(dsens,dim='ens')
        
        
        # Do Preprocessing
        dsens = fix_mesaclip_lr_date(dsens) # Fill missing date for lr (2006-01-02)
        dsens = flatten_ens_time(dsens)     # Concatenate each ens member along time
        
        # Do calculations
        istartens = dsens.ens_start
        ds        = dsens.mergebyens
        
        years   = ds.time.dt.year
        ystart  = years[0].data.item()
        yend    = years[-1].data.item()
        
        # Looping for each doy
        if winsize == 0:
            for doybase in doys_sel:
                
                base_var    = ds.sel(time=ds.time.dt.dayofyear.isin(doybase))
                lagcorr_day = []
                for ll in tqdm(range(nlags)):
                    
                    # Select Lag Variable
                    doylag = doybase + ll
                    if doylag > 365:
                        doylag = doylag % 365
                        year_cross = True
                    else:
                        year_cross = False
                    
                    lag_var  = ds.sel(time=ds.time.dt.dayofyear.isin(doylag))
                    
                    if year_cross:
                        base_slice  = ["%04i-01-01" % (ystart), "%04i-12-31" % (yend-1)] 
                        lag_slice   = ["%04i-01-01" % (ystart+1), "%04i-12-31" % (yend)] 
                        
                        base_var_in = base_var.sel(time=slice(*base_slice))
                        lag_var_in  = lag_var.sel(time=slice(*lag_slice))
                    
                    else:
                        base_var_in = base_var
                        lag_var_in  = lag_var
                    
                    lag_var_in['time'] = base_var_in['time']
                    
                    corrout = xr.corr(base_var_in,lag_var_in,dim='time')
                    lagcorr_day.append(corrout)
                    
                lagcorr_day = xr.concat(lagcorr_day,dim='lag')
                
                outname = "%sdaily_acf_lagmax%i_winsize%02i_block%04i_doy%03i.nc" % (outdir_temp,nlags,winsize,nb,doybase)
                lagcorr_day.to_netcdf(outname)
                
        else:
            print("Other Winsize Not Supported")
        
        nb += 1
        # End y block loop
    # End x block loop
    
                
            
            
            
        
            
            
    
        
        
        
        
        
        
        

#%%





#%%


# #%% Open view of OISST (8GB for global)

# st     = time.time()
# ncpath = "/home/niu4/gliu8/share/OISST/mergetest/anom_detrend2_19820101-20251231/"
# ncname = "oisst_day_1_sst_anom.nc"
# dsview = xr.open_dataset(ncpath+ncname).load()
# xrname = '__xarray_dataarray_variable__'
# dsview = dsview.convert_calendar('noleap') # Remove Leap Year
# print("Loaded Data in %.2fs" % (time.time()-st))
# print(dsview)

# 

# winsize = 5

#%%


#doyall  = ds.time.dt.dayofyear
nlags   = 365
lags    = np.arange(nlags+1)
days    = np.arange(1,366)
ds      = dsview.squeeze()[xrname]

years   = ds.time.dt.year
ystart  = years[0].data.item()
yend    = years[-1].data.item()

if winsize == 0:
    
    #lagcorr_alldoys = []
    for doybase in tqdm(np.arange(1,366)):
        
        base_var    = ds.sel(time=ds.time.dt.dayofyear.isin(doybase))
        
        lagcorr_day = []
        for ll in range(nlags):
            
            # Select Lag Variable
            doylag = doybase + ll
            if doylag > 365:
                doylag = doylag % 365
                year_cross = True
            else:
                year_cross = False
            lag_var  = ds.sel(time=ds.time.dt.dayofyear.isin(doylag))
            
            if year_cross:
                base_slice  = ["%04i-01-01" % (ystart), "%04i-12-31" % (yend-1)] 
                lag_slice   = ["%04i-01-01" % (ystart+1), "%04i-12-31" % (yend)] 
                
                base_var_in = base_var.sel(time=slice(*base_slice))
                lag_var_in  = lag_var.sel(time=slice(*lag_slice))
            
            else:
                base_var_in = base_var
                lag_var_in  = lag_var
            
            lag_var_in['time'] = base_var_in['time']
            
            corrout = xr.corr(base_var_in,lag_var_in,dim='time')
            lagcorr_day.append(corrout)
            
        lagcorr_day = xr.concat(lagcorr_day,dim='lag')
        
        outname = "%sdaily_acf_lagmax%i_nowindow_doy%03i.nc" % (outdir,nlags,doybase)
        lagcorr_day.to_netcdf(outname)
else:
    print("Running script with window size")
    
    for doybase in tqdm(np.arange(1,366)):
        
        
        # Construct Base Window
        baseinfo    = proc.construct_window_doy(doybase,winsize,verbose=False)
        
        lagcorr_day = []
        for ll in tqdm(range(nlags)):
            
            # Select Lag Variable
            doylag = doybase + ll
            if doylag > 365:
                doylag = doylag % 365
                year_cross = True
            else:
                year_cross = False
            
            # Construct Lag Widnow
            laginfo = proc.construct_window_doy(doylag,winsize,verbose=False)
            
            # Perform Indexing
            basevar,lagvar = proc.index_year_base_lag_xr(ds,baseinfo,laginfo,ll,verbose=False)
            
            # Calcualte xrcorr
            lagvar['time'] = basevar['time']
            corrout = xr.corr(basevar,lagvar,dim='time')
            lagcorr_day.append(corrout)
        
        lagcorr_day = xr.concat(lagcorr_day,dim='lag')
        outname = "%sdaily_acf_lagmax%i_winsize%02i_doy%03i.nc" % (outdir,nlags,winsize,doybase)
        lagcorr_day.to_netcdf(outname)
        
        # End DOY Loop

print("Completed Script in %.2fs" % (time.time()-st_all))
    

# #%% Quick Test if I just use xrcorr

# doy_monstart = proc.get_doy_monstart()

# st          = time.time()
# ll      = 1

# doybase = ds.sel(time=ds.time.dt.dayofyear.isin(1))

# print("Caluclated Lag in %.2fs" % (time.time()-st))

