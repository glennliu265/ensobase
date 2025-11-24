#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Using radiative kernels computed through [calculate_radiative_kernels_general],
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

#%% User Edits

# Kernel Information
datpath      = "/home/niu4/gliu8/projects/ccfs/regrid_1x1/"
ccf_vars     = ["sst","eis","Tadv","r700","w700","ws10"]
flxname      = "cre"
standardize  = True
add_ucc      = False
selmons_loop = [[12,1,2],[3,4,5],[6,7,8],[9,10,11]] # [[]]#.None
regrid1x1    = True
expname      = "ERA5_1979_2024"

figpath     = "/home/niu4/gliu8/figures/bydate/2025-11-25/"
proc.makedir(figpath)

#%% Load the all months case
outname         = "%s%s_%s_CCFs_Regression_standardize%i_adducc%i.nc" % (datpath,expname,flxname,standardize,add_ucc)
if regrid1x1:
    outname      = proc.addstrtoext(outname,"_regrid1x1",adjust=-1)
dsall = xr.open_dataset(outname).load()


#%% Load each Month
# Load the Kernels (Seasonal)

ds_byseason = []
for selmons in tqdm(selmons_loop):
    
    outname         = "%s%s_%s_CCFs_Regression_standardize%i_adducc%i.nc" % (datpath,expname,flxname,standardize,add_ucc)
    
    if selmons is not None:
        selmonstr = proc.mon2str(np.array(selmons)-1)
        outname      = proc.addstrtoext(outname,"_"+selmonstr,adjust=-1)
    
    if regrid1x1:
        outname      = proc.addstrtoext(outname,"_regrid1x1",adjust=-1)
    
    # Load Seasonal Outputs
    ds = xr.open_dataset(outname).load()
    
    ds_byseason.append(ds)
    
#%% Load each variable

anompath = "/home/niu4/gliu8/projects/common_data/ERA5/regrid_1x1/anom_detrend1/"

dsvars_anom = []

for ccf in tqdm(ccf_vars):
    
    ncname = "%s%s_1979_2024.nc" % (anompath,ccf)
    ds = xr.open_dataset(ncname)[ccf].load()
    dsvars_anom.append(ds)

#%% For each variable, do the multiplication

outpath = "/home/niu4/gliu8/projects/ccfs/regrid_1x1/%s_components/" % flxname
proc.makedir(outpath)


vv    = 0
nccf  = len(ccf_vars)
dtday = 3600*60
for vv in tqdm(range(nccf)):
    
    # Get Variable
    ccfname        = ccf_vars[vv]
    varanom        = dsvars_anom[vv]
    varnom_std     = varanom / varanom.std('time')
    
    # Multiple by the Coefficient
    coeff_allmons  = dsall.coeffs.sel(ccf=ccfname)/dtday
    R_component    = coeff_allmons * varanom
    vname_new      = flxname + "_" + ccfname
    R_component    = R_component.rename(vname_new)
    rname_out      = "%s%s_component.nc" % (outpath,vname_new,)
    R_component.to_netcdf(rname_out)
    
    coeff_seasonal = [ds.coeffs.sel(ccf=ccfname)/dtday for ds in ds_byseason]
    
    R_component_seasonal = []
    for ss in range(4):
        
        selmons = selmons_loop[ss]
        
        varmon     = proc.selmon_ds(varanom.rename({'valid_time':'time'}),selmons)
        varmon_out = varmon * coeff_seasonal[ss]
        R_component_seasonal.append(varmon_out)
        
    R_component_seasonal = xr.concat(R_component_seasonal,dim='time')
    R_component_seasonal = R_component_seasonal.rename(vname_new)
    R_component_seasonal = R_component_seasonal.rename({'time':'valid_time'})
    rname_out      = "%s%s_component_seasonal.nc" % (outpath,vname_new,)
    R_component_seasonal.to_netcdf(rname_out)
        
    
