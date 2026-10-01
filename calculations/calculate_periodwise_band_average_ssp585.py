#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

For TCo319_ssp585 Runs, Calculate the Period-wise Band Average for Fluxes over a particular region

Copied from `compare_sst_cre_ssp585_awi.ipynb`

Created on Wed Sep  2 17:42:44 2026

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
import scipy as sp

from tqdm import tqdm

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

tbxpath = "/home/niu4/gliu8/scripts/commons/tbx"
sys.path.append(tbxpath)
import tbx as tbx

#%% Set Variable Properties
vnames  = ["cre","tscre","ttcre","sst"]
expnames  = ["TCo319-DART-ssp585d-gibbs-charn",
             "TCo319_ssp585",
             "TCo319_ssp585_ens01",
             "TCo319_ssp585_ens02",
             "TCo319_ssp585_ens03",
             'TCo319-DART-ctl1950d-gibbs-charn',
             'TCo319_ctl1950d']


expcolors = [
    "blue",
    "lightsalmon",
    "orange",
    "darkviolet",
    "gold",
    "k"
    ]

nexps   = len(expnames)
nvars   = len(vnames)

outpath = "/home/niu4/gliu8/projects/ccfs/metrics/regrid_1x1/TCo319_ssp585/SEP_band_avg/"


# Bounding Box Selection
bbox_sep     = [-90+360,-75+360,-40,-15]        # Southeast Tropical Pacific Box from Kang et al. 2026
bbox_sep_new = [-90+360,-75+360,-25,-15]        # Restricted Box based on previous Analysis
bbsel        = bbox_sep_new
bbname       = "SEPnew"

nyr_window   = 40
nsmooth      = 5

#%% Load Variables

varsbyexp = []
for ex in tqdm(range(nexps)):
    expname = expnames[ex]

    vbv = []
    for vv in range(nvars):
        vname = vnames[vv]
        ds    = ut.loadregrid(expname,vname,reformat=True)
        vbv.append(ds)
    varsbyexp.append(vbv)

varsawi = [xr.merge(ds) for ds in varsbyexp]



#%% Get Fluxes

#for vv in range(3):
vv    = 0

for vv in tqdm(range(3)):
    # Retrieve Fluxes
    vname        = vnames[vv]
    flxs         = [ds[vname] for ds in varsawi]
    
    # Preprocess
    flxsanom     = [proc.xrdeseason(ds) for ds in flxs]
    flxsanom     = [proc.xrdetrend_dim(ds,dim='time',deg=2) for ds in flxsanom]
    
    # Take area average
    flxaavg      = [proc.aavg(ds,bbsel) for ds in flxsanom]
    
    
    # Do Sliding Spectra Calculations
    specsbyexp   = [ut.sliding_spectra(ds,nyr_window,nsmooth,detrend=True) for ds in flxaavg]
    
    # Calculate Band Sum
    bandsum_byexp = [ut.band_avg_spectra(ds,debug=False,band_sum=True) for ds in specsbyexp]
    
    for ex in range(nexps):
        outname = "%sBand_Sum_%s_%s_%s_winlen%03i_nsmooth%02i.nc" % (outpath,expnames[ex],vname,bbname,nyr_window,nsmooth)
        bandsum_byexp[ex].to_netcdf(outname)
    
    





