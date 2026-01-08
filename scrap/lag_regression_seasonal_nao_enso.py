#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Lag Regression Between ENSO and selected variable

Created on Thu Dec 11 15:18:21 2025

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

#%% Load the Variable

anompath = "/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1/"
expname  = "TCo319_ssp585"
vname    = "pr"

expnames = ["TCo1279-DART-1950","TCo1279-DART-2090",] #["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090",]

for expname in expnames:
    
    ncname   = "%s%s_%s.nc" % (anompath,expname,vname)
    
    
    st       = time.time()
    dsanom   = xr.open_dataset(ncname)[vname].load()
    proc.printtime(st)
    
    dsanom   = ut.standardize_names(dsanom)
    
    #%% Load ENSO and NAO
    
    # Load ENSO
    standardize  = False
    ensoid_name  = "nino34"
    ensoid       = ut.load_ensoid(expname,ensoid_name,standardize=standardize)
    
    # Load NAO
    naopath      = "/home/niu4/gliu8/projects/scrap/nao_crop/"
    ncname_eof   = naopath + "%s_nao_index.nc"
    naoall       = xr.open_dataset(ncname_eof % (expname)).load()
    ncname_eof2  = naopath + "%s_nao_index_DJF.nc"
    naodjf       = xr.open_dataset(ncname_eof2 % (expname)).load()
    
    
    
    #%% Calculate Instantaneous Regression
    
    selmon       = [12,1,2]
    selmonstr    = proc.mon2str(np.array(selmon)-1)
    
    naoin        = naodjf.pcs.isel(mode=0)
    ensoin       = proc.selmon_ds(ensoid,[12,1,2]) # always Use Wintertime ENSO
    
    
    dsmon        = proc.selmon_ds(dsanom,selmon)
    
    
    lags         = np.array([-3,0,3])
    rrnao        = ut.calc_leadlag_regression_2d(naoin,dsmon,lags,sep_mon=False)
    rrenso       = ut.calc_leadlag_regression_2d(ensoin,dsmon,lags,sep_mon=False)
    
    outpath_regr = "/home/niu4/gliu8/projects/awi_hackathon/regression_patterns/"
    
    outnc_nao    = "%s%s_%s_regression_NAO_%s.nc" % (outpath_regr,expname,vname,selmonstr)
    outnc_enso   = "%s%s_%s_regression_ENSO_%s.nc" % (outpath_regr,expname,vname,selmonstr)
    
    rrnao.to_netcdf(outnc_nao)
    rrenso.to_netcdf(outnc_enso)













