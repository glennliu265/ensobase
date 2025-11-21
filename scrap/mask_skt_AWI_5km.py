#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 18:38:29 2025

@author: gliu
"""



# Load and apply mask to skt


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

# ==============================================================================
#%% Do this for 5km simulations
datpath = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo2559-DART-1950C_5km/atm/"
nc      = datpath + "TCo2559-DART-1950C_atm_10256x5120_1m_skt_1950-1959.nc"

st      = time.time()
ds_sst  = xr.open_dataset(nc)#.load()
ds_sst  = ds_sst.chunk('auto')
print("Data Loaded in %.2fs" % (time.time()-st))


#%% Load  Mask

landmask = ut.load_land_mask_awi("TCo2559-DART-1950C")

latsst   = ds_sst.lat
lonsst   = ds_sst.lon

landmask['lon'] = lonsst
landmask['lat'] = latsst

st = time.time()
ds_sst_masked = ds_sst.skt * landmask
print("Masked in %.2fs" % (time.time()-st)) 

st            = time.time()
ds_sst_masked = ds_sst_masked.rename('sst')
edict         = proc.make_encoding_dict(ds_sst_masked)
outname       = "/home/niu4/gliu8/projects/scrap/processed_global/TCo2559-DART-1950C_sst.nc"
ds_sst_masked.to_netcdf(outname,encoding=edict)
print("Saved in %.2fs" % (time.time()-st)) 

#%% Umm check...

outname       = "/home/niu4/gliu8/projects/scrap/processed_global/TCo2559-DART-1950C_sst.nc"
ds_load       = xr.open_dataset(outname)

st = time.time()
ds_load = ds_load.load()

# =============================================================================
#%% Repeat Process for 31km simulation (Control)

# Control
# datpath         = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_control/"
# nc              = datpath + "ctl1950d_atm_remapped_1m_skt_1850-2134.nc"
# outname         = "/home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_sst.nc"

# Future
datpath         = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/"
nc              = datpath + "TCo319_ssp585_skt_1m_2015-2114_atmogrid.nc"
outname         = "/home/niu4/gliu8/projects/scrap/processed_global/TCo319_ssp585_sst.nc"



# Load Data
st      = time.time()
ds_sst  = xr.open_dataset(nc).load()
print("Data Loaded in %.2fs" % (time.time()-st))

# Load Mask and replace Lat/Lon
landmask = ut.load_land_mask_awi("TCo319")
latsst   = ds_sst.lat
lonsst   = ds_sst.lon
landmask['lon'] = lonsst
landmask['lat'] = latsst

# Apply Mask
st = time.time()
ds_sst_masked = ds_sst.skt * landmask
print("Masked in %.2fs" % (time.time()-st)) 

# Save Output
st            = time.time()
ds_sst_masked = ds_sst_masked.rename('sst')
edict         = proc.make_encoding_dict(ds_sst_masked)
ds_sst_masked.to_netcdf(outname,encoding=edict)
print("Saved in %.2fs" % (time.time()-st)) 

# =============================================================================
#%% Repeat for 31km Future





