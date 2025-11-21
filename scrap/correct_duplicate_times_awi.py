#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Correct Duplicate Times

Created on Thu Nov  6 11:54:35 2025

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

outpath = "/home/niu4/gliu8/projects/common_data/awi_cm3/"
datpath = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_control/"

vnames  = ["tsr","ttr","tsrc","ttrc"]

tstart  = '1950-01-01'
tend    = '2100-01-01'

#%%
# Took between 183.05s to 154.53 sec
for vname in tqdm.tqdm(vnames):
    st     = time.time()
    ncname = "%sctl1950d_atm_remapped_1m_%s_1850-2134.nc" % (datpath,vname)
    ds     = xr.open_dataset(ncname)
    ds     = ds.sel(time_counter=slice(tstart,tend))
    ds     = ds[vname].chunk()
    
    print(ds.shape)
    
    # Remove Duplicate Times and rename
    #ds     = ds.rename({'time_counter':'time'})
    #ds     = ut.remove_duplicate_times(ds)
    
    # Note accidentally made c lowercase
   # outname = "%sTco319_ctl1950d_%s_1950-2100.nc" % (outpath,vname)
    #edict   = proc.make_encoding_dict(ds)
    #ds.to_netcdf(outname,encoding=edict)
    #print("Processed variable in %.2fs" % (time.time()-st))
    
    
    
mv awi_cm3Tco319_ctl1950d_ttrc_1950-2100.nc Tco319_ctl1950d_ttrc_1950-2100.nc
mv awi_cm3Tco319_ctl1950d_ttr_1950-2100.nc Tco319_ctl1950d_ttr_1950-2100.nc



    