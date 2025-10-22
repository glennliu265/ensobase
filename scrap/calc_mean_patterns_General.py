#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

General version of "calc_mean_patterns_TP"
Written to work with regridded AWI-CM3 output...

things to compute:
    - Monthly Mean
    - Monthly Variance
    

Created on Tue Oct  7 16:26:49 2025

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


#%% 

# Simulation Names -----
expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"]

timecrops       = [[1950,2100],None,None,None,None]


expnames_all    = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
nexps_all       = len(expnames_all)

#%% Load Variables (processed in cdo by crop_TP_AWI_variableloop.sh)

outpath  = "/home/niu4/gliu8/projects/scrap/summary_stats/"
expname  = "TCo2559-DART-1950C" #"TCo319_ssp585"
vname    = "sst"
timecrop = None #[1950,2100]


#datpath = "/home/niu4/gliu8/projects/scrap/TP_crop/"
#vnames  = ["lsp"]#["cp","sshf","slhf"]#"ttr","ttrc","tsr","tsrc"]
#vname   = "tx_sur" #"lcc"


nclist = ut.get_rawpath_awi(expname,vname,ensnum=None)
#print(nclist)

ncname = nclist[0]

#%%

try:
    ds = xr.open_dataset(ncname).load()
except:
    print("Could not find %s in %s" % (vname,expname))

# Correct Upper Case Variable Names
if vname.upper() in list(ds.keys()):
    print("Renaming %s to %s for %s" % (vname.upper(),vname,expname))
    ds = ds.rename({vname.upper():vname})
    
# Standardize the names
ds = ut.standardize_names(ds)

# Crop time (mostly for control run, pre 1950)
if timecrop is not None:
    print("Cropping time for %s: %s to %s" % (expname,str(timecrop[0])+'-01-01',str(timecrop[1])+'-12-31'))
    ds = ds.sel(time=slice(str(timecrop[0])+'-01-01',str(timecrop[1])+'-12-31'))

startyr = ds.time.dt.year.isel(time=0).data.item()
endyr   = ds.time.dt.year.isel(time=-1).data.item()

timestr  = "%04ito%04i" % (startyr,endyr)
    
print(expname)
print("\tStart: %s" % ds.time.isel(time=0))
print("\tEnd : %s" % ds.time.isel(time=-1))


    
# Check the Shapes
print("%s shape=%s" % (expname,ds[vname].shape))

#%% Compute Means (both time and monthly)

ds_scycle   = ds.groupby('time.month').mean('time') 
ds_timemean = ds.mean('time')
ds_monvar   = ds.groupby('time.month').var('time')
ds_timevar  = ds.var('time')
outcalcs    = ["time_mean","scycle","variance","monvar"]

outall      = [ds_timemean,ds_scycle,ds_timevar,ds_monvar]

ncalcs      = len(outcalcs)
ds_out      = xr.merge([outall[nn][vname].rename(outcalcs[nn]) for nn in range(ncalcs)])


outname = "%s%s_%s_%s.nc" % (outpath,expname,vname,timestr)
edict   = proc.make_encoding_dict(ds_out)
ds_out.to_netcdf(outname)
print("Saved output to %s" % outname)


