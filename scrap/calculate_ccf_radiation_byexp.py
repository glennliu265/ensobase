#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Copied from [calculate_ccf_radiation_era5]

Compute radiative component for each predictor using the kernels

Created on Fri Nov 28 17:32:24 2025

@author: gliu
"""

import sys
import time
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import xarray as xr
import sys
import glob 
import scipy as sp
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
from scipy.io import loadmat
import matplotlib as mpl
import climlab

import importlib
from tqdm import tqdm

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% User Edits

# Experiment Information
expname      = "CERES_FBCT_ERA5" #"CERES_EBAF_ERA5"#"ERA5_1979_2024"
flxname      = "tscrelowcloud"#"creln"
ccf_vars     = ["sst","eis","Tadv","r700","w700","ws10"]
tstart       = '2002-07-01'
tend         = '2023-02-01'

# Set Paths
kernel_path  = '/home/niu4/gliu8/projects/ccfs/kernels/regrid_1x1/%s/' % expname
ccf_path     = '/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/%s/anom_detrend1/' % expname
outpath      = "/home/niu4/gliu8/projects/ccfs/radiative_components/regrid_1x1/%s/" % expname
proc.makedir(outpath)

# Kernel Information
standardize  = True

calc_seasonal = True
selmons_loop  = [[12,1,2],[3,4,5],[6,7,8],[9,10,11]] # No Need to put None, All Months is loaded separately


#add_ucc      = False

#regrid1x1    = True


figpath     = "/home/niu4/gliu8/figures/bydate/2025-12-02/"
proc.makedir(figpath)

#%% Load the all months case

ncname_kernel   = "%s%s_kernels_standardize%i.nc" % (kernel_path,flxname,standardize)
dsall           = xr.open_dataset(ncname_kernel).load()





#%% Load each Month
# Load the Kernels (Seasonal)
if calc_seasonal:
        
    ds_byseason = []
    for selmons in tqdm(selmons_loop):
        
        ncname_kernel   = "%s%s_kernels_standardize%i.nc" % (kernel_path,flxname,standardize)

        if selmons is not None:
            selmonstr = proc.mon2str(np.array(selmons)-1)
            ncname_kernel   = proc.addstrtoext(ncname_kernel,"_"+selmonstr,adjust=-1)
        
        # Load Seasonal Outputs
        ds = xr.open_dataset(ncname_kernel).load()
        
        ds_byseason.append(ds)
 
#%% Load each variable


dsvars_anom = []

for ccf in tqdm(ccf_vars):
    
    ncname = "%s%s.nc" % (ccf_path,ccf)
    ds     = xr.open_dataset(ncname)[ccf].load()
    

    
    # Do some cropping (copied from calculate_radiative_kernels_byexp)
    tstart = str(dsall.time[0].data)[:7]
    tend   = str(dsall.time[-1].data)[:7]
    ds = ds.sel(valid_time=slice(tstart,tend))
    
    
    dsvars_anom.append(ds)
    
    

#%% For each variable, do the multiplication

# save_files = False # Set to True to Output NetCDFs
# print("Save Files is set to [%s]" % save_files)

proc.makedir(outpath)

vv    = 0
nccf  = len(ccf_vars)
dtday = 3600*60
for vv in tqdm(range(nccf)):
    
    # Get Variable
    ccfname        = ccf_vars[vv]
    varanom        = dsvars_anom[vv]
    
    varanom_std     = varanom / varanom.std('valid_time') # Standardize Variable
    
    # Multiple by the Coefficient
    coeff_allmons  = dsall.coeffs.sel(ccf=ccfname)#/dtday
    R_component    = coeff_allmons * varanom_std
    vname_new      = ccfname
    R_component    = R_component.rename(vname_new)
    rname_out      = "%s%s_%s_component.nc" % (outpath,flxname,vname_new,)
    R_component.to_netcdf(rname_out)
    
    if calc_seasonal:
        coeff_seasonal = [ds.coeffs.sel(ccf=ccfname) for ds in ds_byseason] # /dtday
        R_component_seasonal = []
        for ss in range(4):
            
            selmons    = selmons_loop[ss]
            
            varmon     = proc.selmon_ds(varanom_std.rename({'valid_time':'time'}),selmons)
            varmon_out = varmon * coeff_seasonal[ss]
            R_component_seasonal.append(varmon_out)
            
        R_component_seasonal = xr.concat(R_component_seasonal,dim='time')
        R_component_seasonal = R_component_seasonal.sortby('time')
        R_component_seasonal = R_component_seasonal.rename(vname_new)
        R_component_seasonal = R_component_seasonal.rename({'time':'valid_time'})
        rname_out      = "%s%s_%s_component_seasonal.nc" % (outpath,flxname,vname_new,)
        R_component_seasonal.to_netcdf(rname_out)


# #%% More Debugging

# latf = 50
# lonf = 330

# raw = proc.selpt_ds(varanom)
# allmon = proc.selpt_ds(R_component,lonf,latf)
# varymon = proc.selpt_ds(R_component_seasonal,lonf,latf)


# #%%

# fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(12.5,4))

# #ax.plot(proc.selpt_ds(varmon,lonf,latf),label="Variable")
# #ax.plot(proc.selpt_ds(varmon_out,lonf,latf),label="Multiplied",ls='dashed',c='k')


# ax.plot(proc.selpt_ds(varanom,lonf,latf),label="Variable")

# ax.plot(proc.selpt_ds(R_component,lonf,latf),label="Multiplied Constant",ls='solid',c='gray')
# ax.plot(proc.selpt_ds(R_component_seasonal,lonf,latf),label="Multiplied Seasonal",ls='dashed',c='k')

# plt.show()


# #%% Coeffs
# raw      = proc.selpt_ds(varanom,lonf,latf)
# coeffall = proc.selpt_ds(coeff_allmons,lonf,latf)
# coeffmon = [proc.selpt_ds(ss,lonf,latf) for ss in coeff_seasonal]