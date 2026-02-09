#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Make ENSO Index as a Predictor

Broadcast Nino3.4 to a lon x lat x time map for easy input into the radiative
kernel calculation scripts

Created on Mon Feb  9 11:21:45 2026

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
import importlib


#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut


#%%


expname = "CERES_EBAF_ERA5_2001_2024"
# Path to Anomalized Data
anompath = "/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/%s/anom_detrend1/" % expname
refvar  = "sst"

ninopath = "/home/niu4/gliu8/projects/scrap/nino34/"
ninonc   = "ERA5_1979_2024_nino34.nc"
ninoname = "nino34"

#%% Load Reference Dataset

ds   = xr.open_dataset(anompath+"sst.nc").load()
ds   = ut.standardize_names(ds)
lon  = ds.lon
lat  = ds.lat



#%% Load ENSO

dsenso = xr.open_dataset(ninopath + ninonc).load()

ensoid = dsenso['sst'] * dsenso['std']


#%% Match time and Broadcast ENSO Index in Space
# Note: Land Points are also included...

ensoin,dsin = proc.match_time_month(ensoid,ds.sst)

ensomap = xr.ones_like(dsin) 
ensomap = ensomap * ensoid

#%% Save Output

outname = "%s%s.nc" % (anompath,ninoname,)
ensomap.to_netcdf(outname)

#%% Debug
fig,ax = plt.subplots(1,1)
ensomap.isel(lon=-22,lat=22).plot(ax=ax)
ensomap.isel(lon=-20,lat=-56).plot(ax=ax,ls='dashed')
plt.show()


