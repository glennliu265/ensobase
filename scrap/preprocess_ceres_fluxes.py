#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

--- --- --- --- --- --- --- --- --- --- --- --- --- --- --
Preprocess Ceres Fluxes for Radiative Kernel Calculations
    copied sections from [compare_ERA5_CERES.ipynb]

Created on Mon Nov 24 16:05:11 2025

@author: gliu
--- --- --- --- --- --- --- --- --- --- --- --- --- --- --

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
#import climlab
import cartopy.feature as cfeature
from tqdm.notebook import tqdm

import pandas as pd

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut


#%%
cerespath   = "/home/niu4/gliu8/share/CERES/"
ncname      = "CERES_EBAF-TOA_Ed4.2.1_Subset_200003-202508.nc"
dsceres = xr.open_dataset(cerespath+ncname)

#%% Load fluxes

ceres_to_era_dict= {
    "toa_net_all_mon"   : "allsky",
    "toa_net_clr_c_mon" : "clearsky",
    "toa_sw_clr_c_mon"  : "tsrc",
    "toa_sw_all_mon"    : "tsr",
    "toa_lw_clr_c_mon"  : "ttrc",
    "toa_lw_all_mon"    : "ttr",
    }


#%%

outpath         = "/home/niu4/gliu8/share/CERES/processed/"
ceres_flxnames  =  list(ceres_to_era_dict.keys())
era5_flxnames   = [ceres_to_era_dict[vv] for vv in ceres_flxnames]

vv              = 0
for vv in tqdm(range(len(ceres_flxnames))):
    cname           = ceres_flxnames[vv]
    ename           = era5_flxnames[vv]
    dsin            = dsceres[cname].load()
    
    # Move timestep to first of month (to match ERA5)
    oldtime         = dsin['time']
    tstart          =  str(oldtime[0].data)[:7] + "-01"
    tend            =  str(oldtime[-1].data)[:7] + "-01"
    newtime         = pd.date_range(start=tstart,end=tend,freq="MS")
    dsin['time']    = newtime
    print("New time dimension between %s and %s" % (dsin.time[0].data,dsin.time[-1].data))
    
    # Rename
    dsin            = dsin.rename(ename)
    outname         = "%sCERES_EBAF_%s_2000-03_to_2025-08.nc" % (outpath,ename)
    dsin.to_netcdf(outname)
    print("Saved output to %s" % outname)

#%%

allsky_ceres   = dsceres.toa_net_all_mon.load()
clearsky_ceres = dsceres.toa_net_clr_c_mon.load()


