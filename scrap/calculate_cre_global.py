#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Global Version of calculate_cre
Works with out preprocessed by preprocess_ENSO_Global

Created on FRI Oct 24

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

datpath = "/home/niu4/gliu8/projects/scrap/global_anom_detrend2/"
vnames  = ["ttr","tsr","ttrc","tsrc"]
chunkdict = dict(lat="auto",lon="auto")


#%% Compute Shortwave and Longwave CRE

def remove_duplicate_times(ds,verbose=True):
    # From : https://stackoverflow.com/questions/51058379/drop-duplicate-times-in-xarray
    _, index = np.unique(ds['time'], return_index=True)
    print("Found %i duplicate times. Taking first entry." % (len(ds.time) - len(index)))
    return ds.isel(time=index)

ex       = 3
#calcname = "ttcre" # 
expname  = expnames[ex]

for calcname in ['ttcre','tscre']:
    st = time.time()
            
    if calcname == "tscre":
        varids = [1,3]
    elif calcname == "ttcre":
        varids = [1,2]
    
    
    ds_all   = []
    foundvar = True
    for ii in range(2):
        v     = varids[ii]
        vname = vnames[v]
        
        ncname = "%s%s_%s_anom.nc"% (datpath,expname,vname)
        try:
            ds = xr.open_dataset(ncname)#.load()
            ds = ds.chunk(chunkdict)
        except:
            print("Could not find %s in %s" % (vname,expnames_long[ex]))
            foundvar = False
            ds_all.append(None)
            continue
        
        # Correct Upper Case Variable Names
        if vname.upper() in list(ds.keys()):
            print("Renaming %s to %s for %s" % (vname.upper(),vname,expnames_long[ex]))
            ds = ds.rename({vname.upper():vname})
         
        # Standardize the names
        ds = ut.standardize_names(ds)
        
        # Crop time (mostly for control run, pre 1950)
        timecrop = timecrops[ex]
        if timecrop is not None:
            print("Cropping time for %s: %s to %s" % (expnames_long[ex],str(timecrop[0])+'-01-01',str(timecrop[1])+'-12-31'))
            ds = ds.sel(time=slice(str(timecrop[0])+'-01-01',str(timecrop[1])+'-12-31'))
            
        print(expnames_long[ex])
        print("\tStart: %s" % ds.time.isel(time=0))
        print("\tEnd : %s" % ds.time.isel(time=-1))
        ds_all.append(ds[vname])
    
    
    if foundvar:
        ntimes = [len(ds.time) for ds in ds_all]
        if ntimes[0] != ntimes[1]:
            ds_all = [remove_duplicate_times(ds) for ds in ds_all]
        
        ds_cre  = ds_all[0] - ds_all[1]
        ds_cre  = ds_cre.rename(calcname)
        outname = "%s%s_%s.nc" % (datpath,expnames[ex],calcname)
        print("Saving output to %s" % outname)
        edict   = proc.make_encoding_dict(ds_cre)
        ds_cre.to_netcdf(outname,encoding=edict)
    print("Calculated %s in %.2fs" % (vname,time.time()-st))
    
#%%
