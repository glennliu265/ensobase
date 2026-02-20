#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Temperature Advection for the selected Dataset
copied from [calc_ccfs_regrid.py]

Inputs: sst, u10, v10
Output: Tadv

Notes:
    - For AWI-CM3 9km, preprocessed variables in [preprocess_awi_global_manual.sh]
    - 

Created on Wed Nov 26 11:22:49 2025

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

#%% User Edits

# Options
vnames   = ['sst','u10','v10']#["tx_sur","ty_sur","sst"] #Dont change order
regrid   = False

# Information by Experiment
# 9km Control
expname  = "TCo1279-DART-1950"
datpath  = "/home/niu4/gliu8/projects/scrap/processed_global/"
ncstr    = datpath + expname + "_%s.nc"

# ERA5
expname  = "ERA5_1979_2024"
datpath  = "/home/niu4/gliu8/share/ERA5/processed/"
ncstr    = datpath + "%s_1979_2024.nc"


# TCo319 Control Regridded
expname  = "TCo319_ctl1950d"
datpath  = "/home/niu4/gliu8/projects/scrap/regrid_1x1/"
ncstr    = datpath + expname + "_%s_regrid1x1.nc"


# # TCo319 Control Full
# expname  = "TCo319_ctl1950d"
# datpath  = "/home/niu4/gliu8/projects/scrap/processed_global/"
# ncstr    = datpath + expname + "_%s.nc"


# Calculation Options
calculate_total_Tadv      = False
calculate_components_Tadv = True # Compute ubar dot grad T' and uprime dot grad Tbar

# Load Land Mask
landmask = ut.load_land_mask_awi(expname,regrid=regrid)
print("Will search for %s..." % ncstr)

# =================
#%% Load variables
# =================
dsvars = []
for vv in tqdm.tqdm(range(3)): # 4min 30 for 2 variables, first one took 3 min, total 5min 16s
    vname  = vnames[vv]
    ncname = ncstr % vname
    ds     = xr.open_dataset(ncname)[vname].load()
    dsvars.append(ds)
 
dsall = xr.merge(dsvars)
dsall = ut.standardize_names(dsall)

# =================
#%% Calculate Tadv
# =================

u10       = dsall.u10
v10       = dsall.v10
sst       = dsall.sst
st        = time.time()

ddx,ddy   = ut.calc_grad_centered(sst,latname='lat',lonname='lon')

if calculate_total_Tadv: # Calculate Total Temperature Advection Term
    
    #ddx,ddy   = ut.calc_grad_centered(sst,latname='lat',lonname='lon')
    Tadv2     = - u10.data * ddx.data - v10.data * ddy.data 
    coords    = dict(time=sst.time,lat=sst.lat,lon=sst.lon)
    Tadv2     = xr.DataArray(Tadv2,coords=coords,dims=coords,name="Tadv")
    print("Calculated Tadv in %.2fs" % (time.time()-st))
    
    # =================
    #%% Save Output
    # =================
    
    st         = time.time()
    vname_out  = "Tadv"
    ncname_out = ncstr % vname_out #"%s%s_%s_regrid1x1.nc" % (datpath,expnames[ex],vnames_out[v])
    print("Saving as %s" % ncname_out)
    dsout      = Tadv2.rename(vname_out)
    
    dsout.to_netcdf(ncname_out)
    print("Saved in %.2fs" % (time.time()-st))

elif calculate_components_Tadv:
    
    
    #sst_bar   = sst.groupby('time.month').mean('time')
    #sst_prime = sst.groupby('time.month') - sst_bar
    
    st        = time.time()
    u10_bar   = u10.groupby('time.month').mean('time')
    u10_prime = u10.groupby('time.month') - u10_bar
   
    
    v10_bar   = v10.groupby('time.month').mean('time')
    v10_prime = v10.groupby('time.month') - v10_bar
    
    st        = time.time()
    
    ddx_bar   = ddx.groupby('time.month').mean('time')
    ddx_prime = ddx.groupby('time.month') - ddx_bar
    
    
    ddy_bar      = ddy.groupby('time.month').mean('time')
    ddy_prime    = ddy.groupby('time.month') - ddy_bar
    
    print("Computed in %.2fs" % (time.time()-st))
    
    
    # Mean Advection of Anomalous Gradient
    meanadv_anomgrad = -u10_bar * ddx_prime.groupby('time.month') - v10_bar * ddy_prime.groupby('time.month') 
    vname_out        = "MeanAdvTanom"
    meanadv_anomgrad = meanadv_anomgrad.rename(vname_out)
    ncname_out       = ncstr % vname_out #"%s%s_%s_regrid1x1.nc" % (datpath,expnames[ex],vnames_out[v])
    print("Saving as %s" % ncname_out)
    meanadv_anomgrad.to_netcdf(ncname_out)
    
    # Anomalous Advection of Mean Gradient
    anomadv_meangrad =  -1* (u10_prime.groupby('time.month') * ddx_bar) - v10_prime.groupby('time.month') * ddy_bar
    vname_out        = "AnomAdvTmean"
    anomadv_meangrad = anomadv_meangrad.rename(vname_out)
    ncname_out       = ncstr % vname_out #"%s%s_%s_regrid1x1.nc" % (datpath,expnames[ex],vnames_out[v])
    print("Saving as %s" % ncname_out)
    anomadv_meangrad.to_netcdf(ncname_out)
    
    # Eddying Terms # (Decided not to perform the subtraction, I can do this later...)
    anomadv_anomgrad = - u10_prime * ddx_prime - v10_prime * ddy_prime
    #anomadv_anomgrad = anomadv_anomgrad.groupby('time.month') + anomadv_anomgrad.groupby('time.month').mean('time')
    vname_out        = "AnomAdvTanom"
    anomadv_anomgrad = anomadv_anomgrad.rename(vname_out)
    ncname_out       = ncstr % vname_out #"%s%s_%s_regrid1x1.nc" % (datpath,expnames[ex],vnames_out[v])
    print("Saving as %s" % ncname_out)
    anomadv_anomgrad.to_netcdf(ncname_out)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


    
    


# #%% Debug Plot
# Tadv2mean = Tadv2.mean('time')
# dtday     = 3600*24
# (Tadv2mean*dtday).plot(vmin=-2.5,vmax=2.5,cmap='cmo.balance'),plt.show()

