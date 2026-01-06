#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Script to compute Nino3.4 Index

copied upper section of 

Created on Fri Sep 12 15:12:57 2025

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

#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"
sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% Indicate paths

figpath         = "/home/niu4/gliu8/figures/bydate/2026-01-13/"
proc.makedir(figpath)

outpath         = "/home/niu4/gliu8/projects/scrap/nino34/"

# # Calculate for AWI-CM3
# datpath         = "/home/niu4/gliu8/projects/scrap/TP_crop/"
# expnames        = ["TCo2559-DART-1950C"]#["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090"] # ["glorys"]#
# expnames_long   = ["5km 1950"]#["31km Control","31km SSP585","9km 1950","9km 2090"]
# vname           = "sst"
nclist            = None

# # Calculate for ERA5 as well
# expnames        = ["ERA5_1979_2024",]
# expnames_long   = ["ERA5 (1979-2024)",]
# datpath         = "/home/niu4/gliu8/share/ERA5/processed/"
# vname           = 'sst'
# nclist          = ["sst_1979_2024.nc",]

# Calculate for Updated AWI-CM3 Runs
datpath           = "/home/niu4/gliu8/projects/scrap/processed_global/"
expnames          = ["TCo1279-DART-2060",
                     "TCo319-DART-ctl1950d-gibbs-charn",
                     "TCo319-DART-ssp585d-gibbs-charn",]
vname             = "sst",
nclist            = ["%s_sst.nc" % (ex) for ex in expnames ]


nexps             = len(expnames)

ninoid_name       = 'nino3'#'nino34' # 

bbox_nino34       = [-170+360,-120+360,-5,5]
bbox_nino3        = [-150+360, -90+360 , -5, 5]  # Nino 3 Box: For SST, <tau_x>

if ninoid_name == "nino34":
    bbox = bbox_nino34
elif ninoid_name == 'nino3':
    bbox = bbox_nino3

#%%Functions

def rename_time_dart(ds):
    if 'time_counter' in list(ds.coords):
        print("Renaming [time_counter] to [time]")
        ds = ds.rename({'time_counter':'time'})
    if 'valid_time' in list(ds.coords):
        print("Renaming [valid_time] to [time]")
        ds = ds.rename({'valid_time':'time'})
        
    return ds

def movmean(ds,win):
    return np.convolve(ds.data,np.ones(win)/win,mode='same')

def remake_da(ds,dsref,name='sst'):
    coords = {'time':dsref.time}
    ds_new = xr.DataArray(ds,coords=coords,dims=coords,name=name)
    return ds_new

#%% Load the Dataset and preliminary processing

if nclist is not None:
    dsall         = [xr.open_dataset(datpath+nc).load() for nc in nclist]
else:
    dsall         = [xr.open_dataset("%s%s_%s.nc" % (datpath,ex,vname)).load() for ex in expnames]
dsall         = [rename_time_dart(ds) for ds in dsall]
dsall         = [ut.standardize_names(ds) for ds in dsall]

for dd,ds in enumerate(dsall):
    if ds.lat[1] < ds.lat[0]:
        print("Flipping lat for %s"  % expnames[dd])
        dsall[dd] = proc.fliplat(ds)
#sst_anom     = [proc.xrdeseason(ds) for ds in dsall]

#%% Calculate ENSO Index

apply_movmean = False


# Take Area-weighted average
dsall_reg     = [proc.sel_region_xr(ds,bbox_nino34) for ds in dsall]
nino_sim      = [proc.area_avg_cosweight(ds.sst) for ds in dsall_reg]

# Rename the Time dimension
nino_sim      = [rename_time_dart(ds) for ds in nino_sim]

# Remove seasonal cycle and detrend (quadratically)
nino_ds       = [proc.xrdeseason(ds) for ds in nino_sim] # Remove scycle (anomalize)
nino_ds       = [proc.xrdetrend_1d(ds,2) for ds in nino_ds] # REmove linear trend

# Option to apply a moving mean
if apply_movmean:
    nino34sim = [movmean(ds,5) for ds in nino_ds] # 5-month running mean
else:
    nino34sim = nino_ds

# Normalize
nino34_norm = [ds/np.nanstd(ds) for ds in nino34sim]
nino34_ds   = [remake_da(nino34_norm[ii],dsall[ii]) for ii in range(len(expnames))]


#%% Save Output


for ex in range(nexps):
    
    outname = "%s%s_%s.nc" % (outpath,expnames[ex],ninoid_name)
    print(outname)
    
    outds = nino34_ds[ex]
    outds['std'] = np.nanstd(nino34sim[ex])
    outdict = proc.make_encoding_dict(outds)
    outds.to_netcdf(outname,encoding=outdict)
    
    
    
    







