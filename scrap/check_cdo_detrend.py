#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Check how CDO is detrending

Created on Thu Nov 13 13:21:16 2025

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

#%% Import Custom Modules

amvpath     = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath    = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut


#%%

dpath       = "/home/niu4/gliu8/projects/scrap/regrid_1x1/"
cdopath     = "/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/"
spath       = "/home/niu4/gliu8/projects/scrap/regrid_1x1/scycle/"
ncname      = "TCo2559-DART-1950C_sst_regrid1x1.nc"

#%% Do a manual detrend + deseason

# Load Data
dsraw      = xr.open_dataset(dpath+ncname).load()
dsraw      = dsraw.sst

# Remove Seasonal Cycle
xrscycle   = dsraw.groupby('time.month').mean('time')
dsa        = proc.xrdeseason(dsraw)

# Linear Detrend
dsadt      = proc.xrdetrend(dsa)

#%% now compare with values computed from CDO


cdoscycle  = xr.open_dataset(spath+"TCo2559-DART-1950C_sst.nc").load().sst
cdodetrend = xr.open_dataset(cdopath+"TCo2559-DART-1950C_sst.nc").sst.load()


#%% First, compare the seasonal cycle computed


lonf = -30
latf = 55

mons3  = proc.get_monstr()
fig,ax = viz.init_monplot(1,1)

plotvar = proc.selpt_ds(xrscycle,lonf,latf)
ax.plot(mons3,plotvar,label="Python")
plotvar = proc.selpt_ds(cdoscycle,lonf,latf)
ax.plot(mons3,plotvar,label="cdo",ls='dashed')

ax.legend()
plt.show()


"""

It seems cdo scycle and python scycle are equivalent.

"""

#%% Next, compare the detrended timeseries

lonf    = -30
latf    = 50

fig,ax  = plt.subplots(1,1,constrained_layout=True,figsize=(12.5,4.5))

plotvar = proc.selpt_ds(dsadt,lonf,latf,)
ax.plot(plotvar.time,plotvar,label="Python")

plotvar = proc.selpt_ds(cdodetrend,lonf,latf,)
ax.plot(plotvar.time,plotvar,label="cdo",ls='dashed')

plt.show()

"""

Indeed again, it seems both are equivalent

"""
