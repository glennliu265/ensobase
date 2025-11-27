#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Modified from calculate_ccf_radiation_ERA5

Using radiative kernels computed through [calculate_radiative_kernels] _general],
calculate the radiative effect associated with each variable...


Copied upper section of [visualize_kernels_ERA5]

Created on Fri Nov 21 15:44:00 2025

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


#%%

# Kernel Information
expname         = "TCo2559-DART-1950C"
datpath         = "/home/niu4/gliu8/projects/ccfs/kernels/regrid_1x1/"
anompath        = "/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/"
component_path = "/home/niu4/gliu8/projects/ccfs/radiative_components/"
flxname         = "cre"
add_ucc         = False
selmons_loop    = [[12,1,2],[3,4,5],[6,7,8],[9,10,11]] # [[]]#.None
regrid1x1       = True
standardize     = True
ccf_vars        = ["sst","eis","Tadv","r700","w700","ws10",]#"ucc"] 
ccf_nc          = anompath + expname + "_%s_regrid1x1.nc" # anompath + "%s_1979_2024.nc"

figpath     = "/home/niu4/gliu8/figures/bydate/2025-11-25/"
proc.makedir(figpath)

#%% Load the all months case
outname         = "%s%s_%s_CCFs_Regression_standardize%i_adducc%i.nc" % (datpath,expname,flxname,standardize,add_ucc)
dsall           = xr.open_dataset(outname).load()


#%% Load each Month
# Load the Kernels (Seasonal)

ds_byseason = []
for selmons in tqdm(selmons_loop):
    
    outname         = "%s%s_%s_CCFs_Regression_standardize%i_adducc%i.nc" % (datpath,expname,flxname,standardize,add_ucc)
    
    if selmons is not None:
        selmonstr = proc.mon2str(np.array(selmons)-1)
        outname      = proc.addstrtoext(outname,"_"+selmonstr,adjust=-1)
    
    
    # Load Seasonal Outputs
    ds = xr.open_dataset(outname).load()
    
    ds_byseason.append(ds)
    
#%% Load each variable

dsvars_anom = []

for ccf in tqdm(ccf_vars):
    
    ncname = ccf_nc % (ccf)
    ds = xr.open_dataset(ncname)[ccf].load()
    dsvars_anom.append(ds)

#%% For each variable, do the multiplication

outpath = "%sregrid_1x1/%s/" % (component_path,expname) #"/home/niu4/gliu8/projects/ccfs/regrid_1x1/%s_components/" % flxname
proc.makedir(outpath)

vv    = 0
nccf  = len(ccf_vars)
dtday = 3600*60
for vv in tqdm(range(nccf)):
    
    # Get Variable
    ccfname        = ccf_vars[vv]
    varanom        = dsvars_anom[vv]
    varanom        = ut.standardize_names(varanom)
    
    varanom_std     = varanom / varanom.std('time') # Standardize Variable
    
    # Multiple by the Coefficient
    coeff_allmons  = dsall.coeffs.sel(ccf=ccfname)/dtday
    R_component    = coeff_allmons * varanom_std
    vname_new      = ccfname
    R_component    = R_component.rename(vname_new)
    rname_out      = "%s%s_component.nc" % (outpath,vname_new,)
    R_component.to_netcdf(rname_out)
    
    coeff_seasonal = [ds.coeffs.sel(ccf=ccfname)/dtday for ds in ds_byseason]
    
    R_component_seasonal = []
    for ss in range(4):
        
        selmons = selmons_loop[ss]
        
        varmon     = proc.selmon_ds(varanom_std,selmons)
        varmon_out = varmon * coeff_seasonal[ss]
        R_component_seasonal.append(varmon_out)
    
    R_component_seasonal = xr.concat(R_component_seasonal,dim='time')
    R_component_seasonal = R_component_seasonal.rename(vname_new)
    rname_out      = "%s%s_%s_component_seasonal.nc" % (outpath,flxname,vname_new,)
    R_component_seasonal.to_netcdf(rname_out)
        
    
