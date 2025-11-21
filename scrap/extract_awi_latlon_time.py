#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Extract Latitude, Longitude, and Time dimensions from AWI simulations

Created on Tue Nov  4 16:45:03 2025

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

#%%

outpath = "/home/niu4/gliu8/projects/scrap/awi_common/"

#%% 31 km Control

dpath           = '/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ctl1950d/monthly/'
nc              = "awicm3_tco319_ctl1950d_sst_fesom_r360x180_1950-2100.nc"
dsfesom         = xr.open_dataset(dpath+nc)

nc1             = "awicm3_tco319_ctl1950d_sst_fesom_r360x180_1950-2100.nc"
dsremap         = xr.open_dataset(dpath+nc1)

landmask_remap  = xr.where(np.isnan(dsremap.sst.isel(time=0)),np.nan,1)


outname         = outpath + "Tco319_ctl1950d_r360x180_landmask.nc"
landmask_remap  = landmask_remap.rename('land_mask')
landmask_remap.to_netcdf(outname,)

landmask_31c = xr.open_dataset(outname).load()

# ===============
#%% 31 km Future
# ===============

"""

Takeaway:
    - No land mask available...  (need to find one)
    - 
    
"""


dpath31s = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/"

# Check Regridded Output
# (time_counter: 200, bnds: 2, lat: 180, lon: 360), 311MB
nc      = "TCo319_ssp585_sst_1m_2015-2114_1x1regrid.nc"
dsremap = xr.open_dataset(dpath31s+nc)

# It seems Lat/Lon Matches with 31km Control grid
print(np.all(landmask_31c.lon == dsremap.lon))
print(np.all(landmask_31c.lat == dsremap.lat))

# Now Check (original?) atmo grid
# (time_counter: 1200, bnds: 2, lat: 640, lon: 1312), 4GB
ncatmo   = "TCo319_ssp585_lcc_1m_2015-2114_atmogrid.nc" 
dsatmo31 = xr.open_dataset(dpath31s+ncatmo) 

dstest   = dsatmo31.isel(time_counter=0).lcc.load()# 3MB
dstest.plot()

dstest        = dstest.rename('grid')
dsatmo31_grid = xr.zeros_like(dstest)
outname31s    = outpath + "TCo319_ssp585_atmogrid.nc"
dsatmo31_grid.to_netcdf(outname31s)

# Try SST for ENs Member 1
# (time_counter: 552, lat: 640, lon: 1312)> Size: 2GB
ncpath = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585_ens01/"
ncname           = "TCo319_ssp585_ens01_sst_1m_2055-2100_atmogrid.nc"
dsens1           = xr.open_dataset(ncpath+ncname).sst

dstest           = dsens1.isel(time_counter=0)
dstest1          = dsens1.isel(time_counter=1)
landmask_31km    = xr.where(dstest==dstest1,np.nan,1)
landmask_31km    = landmask_31km.rename("land_mask")
outname          = outpath + "TCo319_ssp585_atm_landmask.nc"
landmask_31km.to_netcdf(outname)

# ===============
#%% DART 1950 9km 
# ===============

# Land 
dpath_9c        = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/"
ncsst           = "TCo1279-DART-1950_atm_remapped_1m_sst_240months.nc"
ds9kmatm        = xr.open_dataset(dpath_9c+ncsst) # 13GB


#ds9kmatm_skt = xr.open_dataset(dpath_9c+"TCo1279-DART-1950_atm_remapped_1m_skt_240months.nc") # 13GB
#dstest_skt   = ds9kmatm_skt.isel(time_counter=0).skt.load()

dstest          = ds9kmatm.isel(time_counter=0).sst.load() # 53mb
dstest_1        = ds9kmatm.isel(time_counter=1).sst.load() # 53mb

# Save the Land Mask
landmask_9km    = xr.where(dstest==dstest_1,np.nan,1)
landmask_9km    = landmask_9km.rename('land_mask')

# Save 
outname         = outpath + "TCo1279-DART-1950_atm_landmask.nc"
landmask_9km.to_netcdf(outname)


# Try to find the land mask using skt and sst differences
# diff = dstest -dstest_skt
# dlist = diff.data.flatten()

# # Try to find threshold
# checkpct = lambda thres: (np.abs(dlist)<thres).sum()/len(dlist) * 100
# thres    = np.hstack([np.array([0,1e-5,1e-4,1e-3,1e-2,]),np.arange(0.1,10.1,0.1)])
# pctthres = np.array([checkpct(th) for th in thres])

# plt.plot(thres,pctthres,marker='d'),plt.show()

# mask_9c_atm = xr.where(diff>0.3,np.nan,1)
# mask_9c_atm.plot(),plt.show()

# #histout = plt.hist(np.abs(diff.data.flatten()),bins=[0,1e-5,1e-4,1e-3,1e-2,1e-1,1,10])
# #sktsst = dstest == dstest_skt
# sktsst = dstest == dstest_skt

# dlist = dstest.data.flatten()
# m     = sp.stats.mode(dlist)
 
# mask_9c_atm = xr.where(dstest==271.46,np.nan,1)
# mask_9c_atm.plot(),plt.show()

# ============= Just Save the Regridded Ocean <Land Mask>
dpath_9c_land       = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/ocn/mon/"
ncssh               = "TCo1279_1950ctl_ssh_1m_1950-1969.nc"
ds9kmocn_regrid     = xr.open_dataset(dpath_9c_land + ncssh)
ds9kmocn_regrid     = ds9kmocn_regrid.ssh.isel(time=0)

maskocn_9kmregrid   = xr.where(np.isnan(ds9kmocn_regrid),np.nan,1)
maskocn_9kmregrid.plot(),plt.show()

maskocn_9kmregrid    = maskocn_9kmregrid.rename('land_mask')
outname9kmregrid_ocn = outpath + "TCo1279_DART-1950_ocn_r3600x1800_landmask.nc"
maskocn_9kmregrid.to_netcdf(outname9kmregrid_ocn )

# ================
#%% DART 2090 9km
# ================

dpath_9kms = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/"
nc2        = "TCo1279_2090slice_aleph_sst_1m_2090-2099.nc"

ds9kms = xr.open_dataset(dpath_9kms + nc2).sst#.isel(time=0).load()

dstests0         = ds9kms.isel(time_counter=0).load()
dstests1         = ds9kms.isel(time_counter=1).load()
landmask_9kms    = xr.where(dstests0==dstests1,np.nan,1)

# It appears there are some differences... maybe relating to sea ice?
check9km_diff    = xr.where(np.isnan(landmask_9km) == np.isnan(landmask_9kms),1,0)

# Save 
landmask_9kms = landmask_9kms.rename('land_mask')
outname = outpath + "TCo1279-DART-2090_atm_landmask.nc"
landmask_9kms.to_netcdf(outname)


# ===========
#%% 5km 1950
# ===========

dpath_5km_atm = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo2559-DART-1950C_5km/atm/"
dpath_5km_ocn = dpath_5km_atm + "../ocn/"

# Make Regridded ocean landmask
ocnnc         = "TCo2559-DART-1950C_ocn_remapped_3600x1800_1m_sst_1950-1958.nc"
dsocn_5km     = xr.open_dataset(dpath_5km_ocn+ocnnc).sst

dstest        = dsocn_5km.isel(time=0)

maskocn_5kmregrid = xr.where(np.isnan(dstest),np.nan,1)
landmask_5km  = maskocn_5kmregrid.rename('land_mask')
outname       = outpath + "TCo2559-DART-1950C_ocn_r3600x1800_landmask.nc"
landmask_5km.to_netcdf(outname)

# Make Regridded atm mask? 
# (time_counter: 120, lat: 5120, lon: 10256) 24 GB
ncatm         = "TCo2559-DART-1950C_atm_10256x5120_1m_slhf_1950-1959.nc"
dsatm_5km     = xr.open_dataset(dpath_5km_atm + ncatm).slhf
# (lat: 5120, lon: 10256) 210MB
slhf          = dsatm_5km.isel(time_counter=0).load()

maskwinter    = xr.where(slhf > 0,np.nan,1)

grid5kmatm = xr.zeros_like(slhf).rename('grid')
outname = outpath + "TCo2559-DART-1950C_atm_grid.nc"
grid5kmatm.to_netcdf(outname)


#%% Check the files

datpath     = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/LSM_and_TOPO/"

resolutions = ["TCo319","TCo1279","TCo2559"]

dsall = []
for rr in resolutions:

    nc          = "%s%s.LSM.nc" % (datpath,rr)
    ds = xr.open_dataset(nc).load().lsm
    dsall.append(ds)
    

maskall = [xr.where(ds==0,1,np.nan) for ds in dsall]
maskall = [ds.rename('land_mask') for ds in maskall]
for ii in range(3):
    outname = "%s%s_atmgrid_original_landmask.nc" % (outpath,resolutions[ii])
    maskall[ii].to_netcdf(outname)

#ds      = xr.open_dataset(datpath+nc)


#%%

def load_land_mask_awi(expname,regrid=False,outpath=None):
    if outpath is None:
        outpath = '/home/niu4/gliu8/projects/scrap/awi_common/'
    if "TCo319" in expname: # 31km Simulations
        print("Loading for 31 km simulations...")
        if regrid: # Load 180x360
            dsmask = "Tco319_ctl1950d_r360x180_landmask.nc"
        else: # Load 640x1213
            dsmask = "TCo319_atmgrid_original_landmask.nc"#"TCo319_ssp585_atm_landmask.nc"
    elif "TCo1279" in expname: # 9km Simulations
        print("Loading for 9 km simulations...")
        if regrid:
            dsmask = "TCo1279_DART-1950_ocn_r3600x1800_landmask.nc"
        else:
            dsmask = "TCo1279_atmgrid_original_landmask.nc"#"TCo1279-DART-1950_atm_landmask.nc"
    elif "TCo2559" in expname: # 5km simulations
        print("Loading for 5 km simulations...")
        if regrid:
            dsmask = "TCo2559-DART-1950C_ocn_r3600x1800_landmask.nc"
        else:
            dsmask = "TCo2559_atmgrid_original_landmask.nc"
    else:
        print("Experiment not found")
        return np.nan
    return xr.open_dataset(dsmask).land_mask.load()
    
expnames            = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
msktest = [load_land_mask_awi(expname) for expname in expnames]

for ii in [0,2,4]:
    msktest[ii].plot()
    plt.show()








