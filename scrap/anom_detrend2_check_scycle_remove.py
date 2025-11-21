#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Scrap to Check if Scycle Removed

Created on Mon Nov  3 10:15:10 2025

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

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%%
lonf    = 330
latf    = 50

#datpath = "/home/niu4/gliu8/projects/scrap/global_anom_detrend2/"
datpath = "/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend2/"
nc      = "TCo319_ctl1950d_clearsky_anom.nc"
vname   = "clearsky"

ds      = xr.open_dataset(datpath + nc)
dspt    = proc.selpt_ds(ds,lonf,latf).load()

scycle  = dspt[vname].groupby('time_counter.month').mean('time_counter')


scycle.plot(),plt.show()


#%% Now check timeseries to see what is going on with all-sky calculations...

dpath  = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_control/"
lonf   = 330
latf   = 50
vnames = ['tsr','ttr']
tsr = xr.open_dataset(dpath+"ctl1950d_atm_remapped_1m_tsr_1850-2134.nc").tsr
ttr = xr.open_dataset(dpath+"ctl1950d_atm_remapped_1m_ttr_1850-2134.nc").ttr

dsall = [tsr,ttr]
dspts = [proc.selpt_ds(ds,lonf,latf).load() for ds in dsall]

#%% Plot the timeseries

fig,axs = plt.subplots(3,1,figsize=(12.5,6))

for vv in range(2):
    ax      = axs[vv]
    plotvar = dspts[vv]
    ax.plot(plotvar.time_counter,plotvar,label=vnames[vv])
    ax.legend()

tdict   = dict(time_counter='time')
ax      = axs[2]
plotvar = ut.remove_duplicate_times(dspts[1].rename(tdict)) + ut.remove_duplicate_times(dspts[0].rename(tdict))
ax.plot(plotvar.time,plotvar,label="allsky")
plt.show()'




#%% Check Seasonal Cycle in newly calculated variable

lonf = 330
latf = 50
ncname = "/home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_allsky.nc"
allsky = xr.open_dataset(ncname)
dspt   = proc.selpt_ds(allsky,lonf,latf).load()
dspta  = dspt.groupby('time.month').mean('time')
dspta  = dspt.groupby('time.month') - dspt.groupby('time.month').mean('time')

