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
expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
# expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"]
timecrops       = [[1950,2100],None,None,None,None]
# vnames          = ['sst',]#["ttcre","tscre"]#['sst','str','ssr','skt','ssh','lcc','tcc','ttr','ttrc','tsr','tsrc']
vnames = ["allsky",'cre','clearsky']
nexps           = len(expnames)

#timecrops       = #[]

#%% Set data and outpath

regrid_1x1      = True

if regrid_1x1:
    print("Using Regridded Output!")
    datpath         = "/home/niu4/gliu8/projects/scrap/regrid_1x1/"
    summarypath     = "/home/niu4/gliu8/projects/scrap/summary_stats/"
    outpath         = datpath + "global_anom_detrend2/"
else:
    print("Using original resolution")
    # datpath is found using ut.get_rawpath_awi
    summarypath     = "/home/niu4/gliu8/projects/scrap/summary_stats/"
    outpath         = "/home/niu4/gliu8/projects/scrap/global_anom_detrend2/"

#%% Load Variable

#expname         = "TCo319_ctl1950d" #"TCo2559-DART-1950C" #"TCo319_ctl1950d" #"TCo2559-DART-1950C" #"TCo1279-DART-1950" #"TCo319_ctl1950d" #"TCo1279-DART-1950" #"TCo2559-DART-1950C" #"TCo319_ssp585"
#timecrop        = None#[1950,2100] # None

for ex,expname in enumerate(expnames):
    
    timecrop = timecrops[ex]
    
    vnames          = ['eis','allsky','cre','clearsky',]#['ttr','ttrc','tsr','tsrc',]#'sst']
    #vname          = "sst"
    for vname in tqdm.tqdm(vnames):
        
        if not regrid_1x1:
            nclist          = ut.get_rawpath_awi(expname,vname,ensnum=None)
            print("Found the following:")
            print("\tTaking the first found file: %s" % nclist[0])
            ncname          = nclist[0]
            
            #% Load Seasonal Cycle for Deseasonalizing (calculated from calc_mean_patterns_General.py)
            # ncsearch        = "%s%s_%s*.nc" % (summarypath,expname,vname)
            # ncsummary       = glob.glob(ncsearch)
            # print(ncsummary)
            # print("Found the following:")
            # ncsummary       = ncsummary[0]
            # print("\tTaking the first found file: %s" % ncsummary)
            # ds_scycle       = xr.open_dataset(ncsummary).scycle # Get seasonal scycle
            
            #%% Open File and Chunk
            dsvar           = xr.open_dataset(ncname)[vname]#[vname]#[vname]
            
        else:
            ncname = "%s%s_%s_regrid1x1.nc" % (datpath,expname,vname)# TCo1279-DART-1950_tsrc_regrid1x1.nc
            try:
                
                dsvar  = xr.open_dataset(ncname)[vname]
            except:
                print("Could not find %s for %s... skipping" % (vname,expname))
                continue
            
        
        # Crop time (mostly for control run, pre 1950)
        if timecrop is not None:
            try:
                print("Cropping time for %s: %s to %s" % (expname,str(timecrop[0])+'-01-01',str(timecrop[1])+'-12-31'))
                dsvar = dsvar.sel(time_counter=slice(str(timecrop[0])+'-01-01',str(timecrop[1])+'-12-31'))
                timename    = "time_counter"
            except:
                print("Renaming variables first")
                # Standardize the names
                dsvar = ut.standardize_names(dsvar)
                print("Cropping time for %s: %s to %s" % (expname,str(timecrop[0])+'-01-01',str(timecrop[1])+'-12-31'))
                dsvar = dsvar.sel(time=slice(str(timecrop[0])+'-01-01',str(timecrop[1])+'-12-31'))
                timename    = "time"
        else:
            if "time_counter" in list(dsvar.coords):
                timename = "time_counter"
            else:
                dsvar = ut.standardize_names(dsvar)
                timename = "time"
        
        dsvar = ut.remove_duplicate_times(dsvar,timename=timename)
        
        # Set up Chunking
        st = time.time()
        if regrid_1x1:
            dsvar= dsvar.load()
        else:
            dsvar           = dsvar.chunk(dict(lat='auto',lon='auto'))
        
        #%% Set up function and preprocess
        def preprocess_enso_point(ts):
            ntime   = len(ts)
            nyr     = int(ntime/12)
            x       = np.arange(ntime)
            
            tsrs    = ts.reshape((nyr,12))
            #tsrs1   = ts.reshape((12,nyr))
            scycle    = tsrs.mean(0)[None,:]
            tsrsa     = tsrs - scycle
            #scyclechk = tsrsa.mean(0)
            
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
            input_core_dims=[[timename],],
            output_core_dims=[[timename],],  # Output Dimension
            vectorize=True,  # True to loop over non-core dims
            dask='parallelized',
            #output_dtypes='float32',
            )
        
        # 520.45s, 492.45, 5875.68s
        
        if not regrid_1x1:
            dsanom = dsanom.compute()
        print("Preprocessed in %.2fs" % (time.time()-st))
        
        #%% Save Output
        st      = time.time()
        outname = "%s%s_%s_anom.nc" % (outpath,expname,vname)
        edict   = proc.make_encoding_dict(dsanom)
        dsanom.to_netcdf(outname,encoding=edict)
        proc.printtime(st,print_str="\tsaved")
        
        dsanom.close()
        dsvar.close()
        


