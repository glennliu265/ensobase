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

# Simulation Names -----
expnames        = ['glorys',]#["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"]

timecrops       = [None,]#[[1950,2100],None,None,None,None]

vnames          = ['sst',]#["ttcre","tscre"]#['sst','str','ssr','skt','ssh','lcc','tcc','ttr','ttrc','tsr','tsrc']

nexps           = len(expnames)
#%% Set data and outpath

datpath         = "/home/niu4/gliu8/projects/scrap/TP_crop/"
outpath         = "/home/niu4/gliu8/projects/scrap/TP_crop/anom_detrend2/"
proc.makedir(outpath)


#%% Load variable for each case

for vname in vnames:
    for ex in tqdm.tqdm(range(nexps)):
        
        print("Preprocessing %s (%s)" % (vname,expnames_long[ex]))
        
        ncname  = "%s%s_%s.nc" % (datpath,expnames[ex],vname)
        try:
            ds      = xr.open_dataset(ncname).load()
        except:
            print("\tERROR: Could not find variable %s for %s... skipping" % (vname,expnames_long[ex]))
            continue
        
        if vname.upper() in list(ds.keys()):
            print("Renaming %s to %s" % (vname.upper(),vname))
            ds = ds.rename({vname.upper():vname})
        
        ds      = ut.standardize_names(ds)
        
        # Crop time (mostly for control run, pre 1950)
        timecrop = timecrops[ex]
        if timecrop is not None:
            print("Cropping time for %s: %s to %s" % (expnames_long[ex],str(timecrop[0])+'-01-01',str(timecrop[1])+'-12-31'))
            ds = ds.sel(time=slice(str(timecrop[0])+'-01-01',str(timecrop[1])+'-12-31'))
        
        st  = time.time()
        dsa = ut.preprocess_enso(ds)
        proc.printtime(st,print_str="\tpreprocessed")
        
        st  = time.time()
        outname = "%s%s_%s_anom.nc" % (outpath,expnames[ex],vname)
        edict   = proc.make_encoding_dict(dsa)
        dsa.to_netcdf(outname,encoding=edict)
        proc.printtime(st,print_str="\tsaved")


