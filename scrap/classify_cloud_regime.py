#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Classify Cloud Regime based on Scott et al. 2020/Medeiros and Stevens 2011
using w700 and EIS

Created on Thu Jan 15 17:20:04 2026

@author: gliu
"""

import sys
import time
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import xarray as xr
import sys
import glob 
import scipy as sp
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
from scipy.io import loadmat
import matplotlib as mpl

import importlib
from tqdm import tqdm

#%%

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut


#%% Load Datasets

datpath = "/home/niu4/gliu8/share/ERA5/processed/"
vnames  = ["w700","eis"]

dsall = []
for vname in tqdm(vnames):
    ncname = "%s%s_1979_2024.nc" % (datpath,vname)
    ds = xr.open_dataset(ncname)[vname].load()
    ds = ut.standardize_names(ds) # Rename valid_time to time
    dsall.append(ds)
dsall = xr.merge(dsall)

#%% Take Time Mean and apply conversion to w700

dsmean  = dsall.mean('time')
dtday   = 3600*24
eis     = dsmean.eis # [K]
w700    = dsmean.w700.squeeze() * dtday # [hpa/day]


lons,lats = np.meshgrid(dsmean.lon,dsmean.lat,)

#%% Also Load Land Ice Mask

icenc = "/home/niu4/gliu8/projects/common_data/ERA5/era5_landmask_fromsst.nc"
dsmask = xr.open_dataset(icenc).load().mask

dsmask = ut.standardize_names(dsmask)

#%% Apply Masks

# Identify Marine Stratocumulus
mask_stratocumulus      = xr.where((w700 > 15) & (eis > 1), 1,0)

# Identify Trade Cumulus
mask_tradecumulus       = xr.where((w700 > 0)  & (eis < 1), 1,0)

# Identify Tropical Ascent
mask_tropicalascent     = xr.where((w700 < 0)  & (np.abs(lats)<25), 1,0)

# Identify Midlatitude
mask_midlatitude        = xr.where(((w700 < 0) | ((w700 < 15) & (eis > 1))) & (np.abs(lats)>25), 1,0)

#%% Apply Masks as Condition

strong_subsidence   = w700 > 15
strong_inversion    = eis  > 1

moderate_subsidence = w700 > 0
weak_inversion      = eis < 1

ascent              = w700 < 0
in_tropics          = np.abs(lats)<25
in_midlats          = np.abs(lats)>25

strong_ascent       = w700 < 15


# Identify Marine Stratocumulus
mask_stratocumulus      = xr.where(strong_subsidence & strong_inversion, 1,0)

# Identify Trade Cumulus
mask_tradecumulus       = xr.where(moderate_subsidence & weak_inversion, 1,0)

# Identify Tropical Ascent
mask_tropicalascent     = xr.where(ascent  & in_tropics, 1,0)

# Identify Midlatitude
mask_midlatitude        = xr.where((ascent | (strong_ascent & strong_inversion)) & in_midlats, 1,0)


#%% Plot each of the masks

inmasks   = [mask_stratocumulus,mask_tradecumulus,mask_tropicalascent,mask_midlatitude]
masknames = ["Marine Stratocumulus","Trade Cumulus","Tropical Ascent","Midlatitude"]



proj = ccrs.PlateCarree()

fig,axs = ut.init_globalmap(2,2,figsize=(12.5,6))


for ii in range(4):
    ax = axs[ii]
    plotvar = inmasks[ii] * dsmask
    ax.set_title(masknames[ii])
    pcm = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj)
    
cb = viz.hcbar(pcm,ax=axs.flatten())
plt.show()


#%% Lets save the output
masknames_var = ["marine_stratocumulus","trade_cumulus","tropical_ascent","midlatitude"]
outmasks = xr.merge([ (inmasks[ii] * dsmask).rename(masknames_var[ii]) for ii in range(4)])

outpath = "/home/niu4/gliu8/projects/ccfs/cloud_types/"
outname = outpath + "ERA5_cloud_regimes_1979_2024.nc"
edict   = proc.make_encoding_dict(outmasks)
outmasks.to_netcdf(outname,encoding=edict)




