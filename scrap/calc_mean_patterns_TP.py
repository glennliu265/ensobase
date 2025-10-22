#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Mean Patterns of Different Variables for AWI-CM3 Simulations

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

datpath = "/home/niu4/gliu8/projects/scrap/TP_crop/"
vnames  = ["lsp"]#["cp","sshf","slhf"]#"ttr","ttrc","tsr","tsrc"]
#vname   = "tx_sur" #"lcc"

for vname in vnames:
    expnames = expnames_all.copy()
    nexps    = len(expnames)
    skipexps = []
    
    ds_all   = []
    for ex in tqdm.tqdm(range(nexps)):
        
        # if vname == "tx_sur": # Have to deal with this later (note I just manually renamed everything)
        #     file_vnames = ["tx_surf","tx_surf_1m"]
        
        ncname = "%s%s_%s.nc" % (datpath,expnames[ex],vname)
        try:
            ds = xr.open_dataset(ncname).load()
        except:
            print("Could not find %s in %s" % (vname,expnames_long[ex]))
            skipexps.append(expnames[ex])
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
            
        ds_all.append(ds)
    
    # Check the Shapes
    [print("%s shape=%s" % (expnames_long[ex],ds_all[ex][vname].shape)) for ex in range(nexps) if expnames[ex] not in skipexps]
    
    #%% Compute Means (both time and monthly)
    
    # Calculate Time Mean, Seasonal Cycle, and Monthly Variance
    ds_scycles  = [ds.groupby('time.month').mean('time') if ds is not None else None for ds in ds_all]
    ds_timemean = [ds.mean('time') if ds is not None else None for ds in ds_all]
    ds_monvar   = [ds.groupby('time.month').var('time') if ds is not None else None for ds in ds_all]
    ds_timevar  = [ds.var('time') if ds is not None else None for ds in ds_all]
    
    outpath     = "/home/niu4/gliu8/projects/scrap/TP_crop/summary/"
    
    outcalcs    = ["mean","scycle","var","monvar"]
    outall      = [ds_timemean,ds_scycles,ds_timevar,ds_monvar]
    # if expnames[0] in skipexps: # Add Dummy variable for None
    #     outall = [[None,] + ou for ou in outall] # This way the indexing is preserved
    
    ncalcs      = len(outcalcs)
    for nn in range(ncalcs):
        for ex in tqdm.tqdm(range(nexps)):
            if expnames[ex] in skipexps:
                print("Not saving file for %s" % (expnames[ex]))
                continue
            outname = "%s%s_%s_%s.nc" % (outpath,expnames[ex],vname,outcalcs[nn])
            #print(outname)
            ds_out  = outall[nn][ex]
            edict   = proc.make_encoding_dict(ds_out)
            ds_out.to_netcdf(outname)
            print("Saved output to %s" % outname)
    
    
    
# For Debugging
# for ex in range(nexps):
#     ds = ds_all[ex]
#     out = ds.groupby('time.month').var('time')



