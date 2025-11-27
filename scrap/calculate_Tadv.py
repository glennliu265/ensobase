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


#%% Debug Plot
Tadv2mean = Tadv2.mean('time')
dtday     = 3600*24
(Tadv2mean*dtday).plot(vmin=-2.5,vmax=2.5,cmap='cmo.balance'),plt.show()

