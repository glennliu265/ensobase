#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Lead Lag regression of TOA Fluxes

- Copied from awi-cm3-toa-leadlag-analysis.ipynb


Created on Fri Feb 27 10:13:57 2026

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

from tqdm import tqdm

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% Plotting Parameters
proj    = ccrs.PlateCarree()
figpath = "/home/niu4/gliu8/figures/bydate/2026-03-03/"
proc.makedir(figpath)

#%%
expnames = [
    "TCo319_ctl1950d",
    "TCo319-DART-ctl1950d-gibbs-charn",
    "TCo319_ssp585",
    "TCo319-DART-ssp585d-gibbs-charn",
    "TCo1279-DART-1950",
    "TCo1279-DART-2060",
    "TCo1279-DART-2090",
    "TCo2559-DART-1950C",
    "CERES_EBAF_ERA5_2001_2024",]


# %%Some Functions

def load_enso_eof(expname,datpath=None,apply_smoothing=True,sep_mon=False):
    # Load EOF-based ENSO computed via `calc_EOF_ENSO.py`
    if datpath is None:
        datpath = "/home/niu4/gliu8/projects/scrap/nino34/"
    ninonc        = "%s%s_enso_eof_rotated.nc" % (datpath,expname)
    if sep_mon:
        ninonc = proc.addstrtoext(ninonc,"_sepmon",adjust=-1)
    dsnino        = xr.open_dataset(ninonc).load()
    
    # Apply 1-2-1 filter
    if apply_smoothing:
        filter_coeffs = [0.25,0.5,0.25]
        ep     = np.convolve(dsnino.ep,filter_coeffs,mode='same')
        cp     = np.convolve(dsnino.cp,filter_coeffs,mode='same')
        tcoord = dict(time=dsnino.time)
        ep     = xr.DataArray(ep,coords=tcoord,dims=tcoord,name='ep')
        cp     = xr.DataArray(cp,coords=tcoord,dims=tcoord,name='cp')
        
        
    else:
        ep     = dsnino.ep
        cp     = dsnino.cp
    return ep,cp

def load_input(vname,expname,anom=False,regrid=True):
    inpath = "/home/niu4/gliu8/projects/ccfs/input_data/"
    if regrid:
        inpath += "regrid_1x1"
    if anom:
        procname = "anom_detrend1"
    else:
        procname = "raw"
    ncname = "%s/%s/%s/%s.nc" % (inpath,expname,procname,vname)
    return xr.open_dataset(ncname)[vname]


#%% Load ENSO

ep_indices = []
cp_indices = []
for expname in expnames:
    if expname == "CERES_EBAF_ERA5_2001_2024":
        expname_in = "ERA5_1979_2024"
        eraflag    = True # Need to Crop to Dataset Timeperiod Later
    else:
        expname_in = expname
        eraflag    = False
    ep,cp = load_enso_eof(expname_in)

    ep_indices.append(ep)
    cp_indices.append(cp)
    
#%% Load Fluxes
flxnames    = ["cre","tscre","ttcre"]


flx_byexp = []
for ex,expname in tqdm(enumerate(expnames)):
    dsflxs = []
    for flxname in flxnames:
        dsvar = load_input(flxname,expname,anom=True,regrid=True).load()
        dsvar = ut.standardize_names(dsvar)
        if "TCo" in expname:
            dsvar = ut.varcheck(dsvar,flxname,expname)
        dsvar = ut.remove_duplicate_times(dsvar)
        
        dsvar,_ = proc.match_time_month(dsvar,ep_indices[ex])
        dsflxs.append(dsvar)

        
    #dsflxs = xr.merge(dsflxs)
    flx_byexp.append(dsflxs)
        
#%% Calculate Lead Lag regression with ENSO
outpath = "/home/niu4/gliu8/projects/ccfs/metrics/regrid_1x1/scrap/"

nexp       = len(expnames)
leadlags   = np.arange(-24,25,1)

#flxts_in   = gmeans # [experiment][variable]
#llbynino    = [] # [nino][experiment][variable]

ninos_in    = [cp_indices,ep_indices]
nino_names  = ["CP","EP"]

#llpatbynino = []
for nn in range(2):
    
    #llpatbyexp = []
    for ex in range(nexp):
        # Get Nino Index
        ninoin = ninos_in[nn][ex]
        
        #llpatbyflx = []
        for vv in tqdm(range(3)):
            st = time.time()
            # Get Flux
            flxinpat = flx_byexp[ex][vv]
            flxinpat,ninoin   = proc.match_time_month(flxinpat,ninoin)
            
            # Compute Lead Lag regression
            lloutpat          = ut.calc_leadlag_regression_2d(ninoin,flxinpat,leadlags)
            
            outname = "%sLag_Regression_%s_%s_%sENSO_Lags%02i.nc" % (outpath,expnames[ex],flxnames[vv],nino_names[nn],leadlags[-1])
            
            lloutpat.to_netcdf(outname)
            print("Completed %s ENSO, %s, %s in %.2fs" % (nino_names[nn],expnames[ex],flxnames[vv],time.time()-st))
            
    #         llpatbyflx.append(lloutpat)
        

    #     llpatbyexp.append(llpatbyflx)
    # llpatbynino.append(llpatbyexp)

    #         continue
    #     continue
    # continue






