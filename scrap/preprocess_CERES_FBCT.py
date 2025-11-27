#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Preprocess CERES Flux by Cloud Type

Split downloaded files and concatenate over time

Created on Wed Nov 26 13:13:01 2025

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
import importlib

from sklearn.linear_model import LinearRegression
import sklearn

#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%%

datpath = "/home/niu4/gliu8/share/CERES/FBCT/"
outpath = "/home/niu4/gliu8/share/CERES/processed/"

#%% Find files

ncsearch = "%sCERES_FluxByCldTyp-MON_Terra-Aqua-MODIS_Ed4.1_Subset_*.nc" % datpath
nclist   = glob.glob(ncsearch)
nclist.sort()

print(nclist)

#%% Open datasets

dsall = xr.open_mfdataset(nclist,concat_dim='time',combine='nested')

#%% Load for a variable

rename_dict = dict(
    cldarea_cldtyp_mon='tcc_cldtyp',# Cloud Area Fraction by Cloud Type [percent]
    #toa_lw_cldtyp_mon ='ttr_cldtyp',    # Observed TOA Longwave Flux - All-sky [W/m2]
    #toa_lw_clr_mon    ='ttrc',   # Observed TOA Longwave Flux - Clear-sky [W/m2]
    toa_sw_cldtyp_mon ='tsr_cldtyp',    # Observed TOA Shortwave Flux - All-sky [W/m2]
    toa_sw_clr_mon ='tsrc',    # Observed TOA Shortwave Flux - All-sky [W/m2]
    #toa_lw_all_mon ='ttr',    # Observed TOA Longwave Flux - All-sky [W/m2]
    #toa_sw_all_mon ='tsr',    # Observed TOA Shortwave Flux - All-sky [W/m2]
    )

ceresnames = list(rename_dict.keys())
era5names  = [rename_dict[nn] for nn in ceresnames]
nvar       = len(ceresnames)
 
for vv in tqdm.tqdm(range(nvar)):
    cname      = ceresnames[vv]
    ename      = era5names[vv]
    outname    = "%sCERES_FBCT_%s_2002-07_to_2023-02.nc" % (outpath,ename)
    st         = time.time()
    dsvar      = dsall[cname].load().rename(ename)
    dsvar.to_netcdf(outname)
    print("Saved [%s] in %.2fs" % (ename,time.time()-st))
    )



