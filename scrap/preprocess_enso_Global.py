#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Global Version of Preprocess_ENSO_Global

Created on Wed Oct 22 14:13:00 2025

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

# # Simulation Names -----
# expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
# expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"]
# timecrops       = [None,]#[[1950,2100],None,None,None,None]
# vnames          = ['sst',]#["ttcre","tscre"]#['sst','str','ssr','skt','ssh','lcc','tcc','ttr','ttrc','tsr','tsrc']
# nexps           = len(expnames)


timecrops = []

#%% Set data and outpath

summarypath     = "/home/niu4/gliu8/projects/scrap/summary_stats/"
outpath         = "/home/niu4/gliu8/projects/scrap/global_anom_detrend2/"

#%% Load Variable

expname         = "TCo2559-DART-1950C" #"TCo319_ssp585"
vname           = "sst"
timecrop        = None #[1950,2100]
nclist          = ut.get_rawpath_awi(expname,vname,ensnum=None)
print("Found the following:")
print("\tTaking the first found file: %s" % nclist[0])
ncname          = nclist[0]

#%% Load Seasonal Cycle for Deseasonalizing (calculated from calc_mean_patterns_General.py)

# ncsearch        = "%s%s_%s*.nc" % (summarypath,expname,vname)
# ncsummary       = glob.glob(ncsearch)
# print(ncsummary)
# print("Found the following:")
# ncsummary       = ncsummary[0]
# print("\tTaking the first found file: %s" % ncsummary)
# ds_scycle       = xr.open_dataset(ncsummary).scycle # Get seasonal scycle

#%% Open File and Chunk

dsvar           = xr.open_dataset(ncname)[vname]

# Standardize the names
dsvar = ut.standardize_names(dsvar)

# Crop time (mostly for control run, pre 1950)
if timecrop is not None:
    print("Cropping time for %s: %s to %s" % (expname,str(timecrop[0])+'-01-01',str(timecrop[1])+'-12-31'))
    dsvar = dsvar.sel(time=slice(str(timecrop[0])+'-01-01',str(timecrop[1])+'-12-31'))

# Set up Chunking
dsvar           = dsvar.chunk('auto')

#%% Set up function and preprocess

def preprocess_enso_point(ts):
    ntime   = len(ts)
    nyr     = int(ntime/12)
    x       = np.arange(ntime)
    tsrs    = ts.reshape((12,nyr))
    tsrsa   = tsrs - tsrs.mean(1)[:,None]
    tsa     = tsrsa.flatten()
    #tsa     = proc.deseason(ts)
    if np.any(np.isnan(ts)):
        return np.ones(tsa.shape)*np.nan
    ydetrended,model=proc.detrend_poly(x,tsa,2)
    return ydetrended

# Set Up Function
dsanom = xr.apply_ufunc(
    preprocess_enso_point,  # Pass the function
    dsvar,  # The inputs in order that is expected
    # Which dimensions to operate over for each argument...
    input_core_dims=[['time'],],
    output_core_dims=[['time'],],  # Output Dimension
    vectorize=True,  # True to loop over non-core dims
    dask='parallelized',
    output_dtypes='float32',
    )

# 520.45s, 492.45
st = time.time()
dsanom = dsanom.compute()
print("Preprocessed in %.2fs" % (time.time()-st))

#%% Save Output

st      = time.time()
outname = "%s%s_%s_anom.nc" % (outpath,expname,vname)
edict   = proc.make_encoding_dict(dsanom)
dsanom.to_netcdf(outname,encoding=edict)
proc.printtime(st,print_str="\tsaved")

# #%% Load variable for each case

# for vname in vnames:
#     for ex in tqdm.tqdm(range(nexps)):
        
#         print("Preprocessing %s (%s)" % (vname,expnames_long[ex]))
        
#         ncname  = "%s%s_%s.nc" % (datpath,expnames[ex],vname)
#         try:
#             ds      = xr.open_dataset(ncname).load()
#         except:
#             print("\tERROR: Could not find variable %s for %s... skipping" % (vname,expnames_long[ex]))
#             continue
        
#         if vname.upper() in list(ds.keys()):
#             print("Renaming %s to %s" % (vname.upper(),vname))
#             ds = ds.rename({vname.upper():vname})
        
#         ds      = ut.standardize_names(ds)
        
#         # Crop time (mostly for control run, pre 1950)
#         timecrop = timecrops[ex]
#         if timecrop is not None:
#             print("Cropping time for %s: %s to %s" % (expnames_long[ex],str(timecrop[0])+'-01-01',str(timecrop[1])+'-12-31'))
#             ds = ds.sel(time=slice(str(timecrop[0])+'-01-01',str(timecrop[1])+'-12-31'))
        
#         st  = time.time()
#         dsa = ut.preprocess_enso(ds)
#         proc.printtime(st,print_str="\tpreprocessed")
        
#         st  = time.time()
#         outname = "%s%s_%s_anom.nc" % (outpath,expnames[ex],vname)
#         edict   = proc.make_encoding_dict(dsa)
#         dsa.to_netcdf(outname,encoding=edict)
#         proc.printtime(st,print_str="\tsaved")


