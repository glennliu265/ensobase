#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Script to Compute the NAO Index


(1) Open data and crop to region
(2) Deseason and Detrend
(3) Compute the NAO Index using NAO Analysis
(4) Save the output


Created on Tue Dec  9 15:19:30 2025

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

#%%

bbox     = [-90,40,20,80] # Assumes with degrees west #[-90+360, 40, 20, 80]
expnames = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090",]#"TCo2559-DART-1950C"]
ncnames  = {
    "TCo319_ctl1950d"   : "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_control/TCo319_ctl1950d_msl_1m_1850-2134.nc",
    "TCo319_ssp585"     : "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/TCo319_ssp585_msl_1m_2015-2114.nc",
    "TCo1279-DART-1950" : "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/TCo1279_DART-1950_msl_1m_1950-1969.nc",
    "TCo1279-DART-2090" : "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/TCo1279_DART-2090_msl_1m_2090-2099.nc",
    }

naocrop_path = "/home/niu4/gliu8/projects/scrap/nao_crop/"

    
    

#%% Open Views of each Dataset

dsall = []
for ex in range(4):
    
    expname = expnames[ex]
    ncname  = "%s%s_msl.nc" % (naocrop_path,expname)
    ds      = xr.open_dataset(ncname)['msl']
    #ds      = xr.open_dataset(ncnames[expname])['msl']
    
    if expname == "TCo319_ctl1950d":
        ds = ds.sel(time_counter=slice('1950-01-01','2100-12-31'))
        print(ds.time_counter)
    #ds = ut.get_rawpath_awi(expname,'msl',)
    
    ds  = proc.lon360to180_ds(ds,lonname='lon')
    
    # Remove seasonal cycle and detrend
    st     = time.time()
    ds     = ut.standardize_names(ds)
    dsa    = proc.xrdeseason(ds)
    dsa_dt = proc.xrdetrend_nd(dsa,2)
    ncout  = "%s%s_msl_anom.nc" % (naocrop_path,expname)
    dsa_dt.to_netcdf(ncout)
    print("Deseason/Detrend in %.2fs" % (time.time()-st))
    
    
    dsall.append(ds)
    
#%% Select NAO Bounding Box Region and Load

#dsreg = [proc.sel_region_xr(ds,bbox) for ds in dsall]


def calc_monthly_eof(daf,bboxeof,N_mode=None,concat_ens=True,mask=None,bbox_check=None):
    
    if 'ens' not in daf.dims: # Add dummy ens variable
        daf = daf.expand_dims(dim={'ens':[0,]},axis=1)
        print("Adding ens dim")
        
    
    flxa     = daf # [Time x Ens x Lat x Lon] # Anomalize variabless
    
    # Apply area weight
    wgt    = np.sqrt(np.cos(np.radians(daf.lat.values))) # [Lat]
    flxwgt = flxa * wgt[None,None,:,None]
    
    # Apply Max if needed
    if mask is not None:
        print("Applying provided mask...")
        flxwgt = flxwgt * mask
    
    # Select Region
    flxreg     = proc.sel_region_xr(flxwgt,bboxeof)
    
    flxout     = flxreg.values
    ntime,nens,nlatr,nlonr = flxout.shape
    
    if concat_ens:
        # IMPORTANT NOTE (implement fix later)
        # Variable must be stacked as [ens x time x otherdims]
        if flxout.shape[0] != nens:
            ens_reshape_flag = True
            print("Warning, since ensemble dimension is NOT first, temporarily permuting array to ens x time")
            flxout = flxout.transpose(1,0,2,3)
        else:
            ens_reshape_flag = False
        print("Stacking Dimensions")
        flxout = flxout.reshape(nens*ntime,1,nlatr,nlonr)
        ntime,nens,nlatr,nlonr = flxout.shape
    npts       = nlatr*nlonr
    nyr        = int(ntime/12)
    if N_mode is None: # Set EOFs to number of years
        N_mode=nyr
    
    # Repeat for full variable
    flxout_full= flxa.values
    _,_,nlat,nlon=flxout_full.shape
    if ens_reshape_flag:
        print("Permuting full variable")
        print("\tOriginal Shape %s" % str(flxout_full.shape))
        flxout_full = flxout_full.transpose(1,0,2,3)
        print("\tNew Shape %s" % str(flxout_full.shape))
    npts_full  = nlat*nlon
    if concat_ens:
        flxout_full = flxout_full.reshape(ntime,1,nlat,nlon)
    print("\tFinal Shape %s" % str(flxout_full.shape))
    
    # Check to see if N_mode exceeds nyrs
    if N_mode > nyr:
        print("Requested N_mode exists the maximum number of years, adjusting....")
        N_mode=nyr
    
    # Preallocate for EOF Analysis
    eofall    = np.zeros((N_mode,nens,nlat*nlon)) * np.nan
    pcall     = np.zeros((N_mode,nens,ntime)) * np.nan
    varexpall = np.zeros((N_mode,nens)) * np.nan
        
    # Loop for ensemble memmber
    for e in tqdm.tqdm(range(nens)):
        
        # Remove NaN Points
        flxens            = flxout[:,e,:,:].reshape(ntime,npts) #  Time x Space
        okdata,knan,okpts = proc.find_nan(flxens,0)
        _,npts_valid = okdata.shape
        
        # Repeat for full data
        flxens_full       = flxout_full[:,e,:,:].reshape(ntime,npts_full)
        okdataf,knanf,okptsf = proc.find_nan(flxens_full,0)
        _,npts_validf = okdataf.shape
        
        # Reshape to [yr x mon x pts]
        okdatar  = okdata.reshape(ntime,npts_valid)
        okdatarf = okdataf.reshape(ntime,npts_validf)
        
        # Calculate EOF by month
        #for im in range(12):
            
        # Compute EOF
        datain          = okdatar.T # --> [space x time]
        eofs,pcs,varexp = proc.eof_simple(datain,N_mode,1)
        
        # Standardize PCs
        pcstd = pcs / pcs.std(0)[None,:]
        
        # Regress back to dataset
        datainf = okdatarf[:,:].T
        eof,b = proc.regress_2d(pcstd.T,datainf.T) # [time x pts]
        
        
        # Save the data
        eofall[:,e,okptsf] = eof.copy()
        pcall[:,e,:] = pcs.T.copy()
        varexpall[:,e] = varexp.copy()
    
    # Reshape the variable
    eofall = eofall.reshape(N_mode,nens,nlat,nlon) # (86, 42, 96, 89)
    
    
    # Flip Signs
    if bbox_check is not None:
        print("Flipping boxes based on [bbox_check]")
        nmode_check = len(bbox_check)
        for N in tqdm.tqdm(range(nmode_check)):
            chkbox = bbox_check[N]
            for e in range(nens):
                for m in range(12):
                    
                    
                    sumflx = proc.sel_region(eofall[N,[m],e,:,:].transpose(2,1,0),flxa.lon.values,flxa.lat.values,chkbox,reg_avg=True)
                    #sumslp = proc.sel_region(eofslp[:,:,[m],N],lon,lat,chkbox,reg_avg=True)
                    
                    if sumflx > 0:
                        print("Flipping sign for NHFLX, mode %i month %i" % (N+1,m+1))
                        eofall[N,m,e,:,:]*=-1
                        pcall[N,m,e,:] *= -1
    else:
        print("Sign of EOF pattern will not be checked.")
    
    startyr   = daf.time.data[0]
    nyrs      = int(len(daf.time)/12)
    if concat_ens:
        tnew      = np.arange(0,int(ntime/12))
    else:
        tnew      = xr.cftime_range(start=startyr,periods=nyrs,freq="YS",calendar="noleap")

    # Make Dictionaries
    coordseof = dict(mode=np.arange(1,N_mode+1),mon=np.arange(1,13,1),ens=np.arange(1,nens+1,1),lat=flxa.lat,lon=flxa.lon)
    daeof     = xr.DataArray(eofall,coords=coordseof,dims=coordseof,name="eofs")

    coordspc  = dict(mode=np.arange(1,N_mode+1),mon=np.arange(1,13,1),ens=np.arange(1,nens+1,1),yr=tnew)
    dapcs     = xr.DataArray(pcall,coords=coordspc,dims=coordspc,name="pcs")

    coordsvar = dict(mode=np.arange(1,N_mode+1),mon=np.arange(1,13,1),)
    davarexp  = xr.DataArray(varexpall,coords=coordsvar,dims=coordsvar,name="varexp")
    
    ds_eof    = xr.merge([daeof,dapcs,davarexp])
    return ds_eof.squeeze()




#%% Preprocess and detrend quadratically

ds     = dsall[0]

# Load in dataset
st     = time.time()
ds     = ds.load()
print("Loaded dataset in %.2fs" % (time.time()-st))

# Remove seasonal cycle and detrend
st = time.time()
ds     = ut.standardize_names(ds)
dsa    = proc.xrdeseason(ds)
dsa_dt = proc.xrdetrend_nd(dsa,2)
print("Deseason/Detrend in %.2fs" % (time.time()-st))

#%%

# Subset Season
selmons   = [12,1,2]
selmonstr = proc.mon2str(np.array(selmons)-1) #"DJF"

# Perform EOF Analysis (NAO)
for ex in tqdm.tqdm(range(4)):
    expname = expnames[ex]
    ncin   = "%s%s_msl_anom.nc" % (naocrop_path,expname)
    dsa_dt = xr.open_dataset(ncin)['msl'].load()
    
    if selmons is not None:
        dsa_dt = proc.selmon_ds(dsa_dt,selmons)
    
    # Apply area weight
    dsa_dt = dsa_dt.transpose('time','lat','lon')
    wgt    = np.sqrt(np.cos(np.radians(dsa_dt.lat.values))) # [Lat]
    dswgt  = dsa_dt.data * wgt[None,:,None]
    
    # Preallocate for EOF Analysis
    N_mode    = 3
    ntime,nlat,nlon = dsa_dt.shape
    npts      = nlat*nlon
    eofall    = np.zeros((N_mode,nlat*nlon)) * np.nan
    pcall     = np.zeros((N_mode,ntime)) * np.nan
    varexpall = np.zeros((N_mode)) * np.nan
    
    # Remove NaN Points
    reshapedata       = dswgt.reshape(ntime,npts) #  Time x Space
    okdata,knan,okpts = proc.find_nan(reshapedata,0)
    _,npts_valid      = okdata.shape
    
    # Repeat for full data
    flxens_full       = dsa_dt.data.reshape(ntime,npts)
    okdataf,knanf,okptsf = proc.find_nan(flxens_full,0)
    _,npts_validf = okdataf.shape
    
    # Perform EOF
    datain            = okdata.T # --> [space x time]
    eofs,pcs,varexp   = proc.eof_simple(datain,N_mode,1)
    
    # Standardize PCs
    pcstd             = pcs / pcs.std(0)[None,:]
    
    # Regress back to dataset
    datainf = okdataf[:,:].T
    eof,b = proc.regress_2d(pcstd.T,datainf.T) # [time x pts]
    
    # # Standardize PCs
    
    # Reshape the variable
    eofall[:,okptsf] = eof
    eofall = eofall.reshape(N_mode,nlat,nlon) # (86, 42, 96, 89)
    pcall = pcs
    varexpall = varexp
    
    # Flip Signs
    spgbox     = [-60,20,40,80]
    eapbox     = [-60,20,40,60] # Shift Box west for EAP
    bbox_check = [spgbox,eapbox,]    
    print("Flipping boxes based on [bbox_check]")
    nmode_check = len(bbox_check)
    for N in tqdm.tqdm(range(nmode_check)):
        chkbox = bbox_check[N]
        
        sumflx = proc.sel_region(eofall[N,:,:].transpose(1,0),dsa_dt.lon.values,dsa_dt.lat.values,chkbox,reg_avg=True)
        #sumslp = proc.sel_region(eofslp[:,:,[m],N],lon,lat,chkbox,reg_avg=True)
        
        if sumflx > 0:
            print("Flipping sign for NHFLX, mode %i" % (N+1))
            eofall[N,:,:]*=-1
            pcall[N,:] *= -1
    
    
    # Make Dictionaries
    flxa      = dsa_dt
    coordseof = dict(mode=np.arange(1,N_mode+1),lat=flxa.lat,lon=flxa.lon)
    daeof     = xr.DataArray(eofall,coords=coordseof,dims=coordseof,name="eofs")
    
    coordspc  = dict(mode=np.arange(1,N_mode+1),time=dsa_dt.time,)
    dapcs     = xr.DataArray(pcall.T,coords=coordspc,dims=coordspc,name="pcs")
    
    coordsvar = dict(mode=np.arange(1,N_mode+1))
    davarexp  = xr.DataArray(varexpall,coords=coordsvar,dims=coordsvar,name="varexp")
    
    ds_eof    = xr.merge([daeof,dapcs,davarexp])
    
    ncout  = "%s%s_nao_index.nc" % (naocrop_path,expname)
    if selmons is not None:
        ncout = proc.addstrtoext(ncout,"_"+selmonstr,adjust=-1)
    ds_eof.to_netcdf(ncout)
    # pcstd = pcs / pcs.std(0)[None,:]
    
    
    

# Regress back to dataset
#datainf = okdatarf[:,:].T
#eof,b   = proc.regress_2d(pcstd.T,datainf.T) # [time x pts]

















