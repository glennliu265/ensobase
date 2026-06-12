#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Perform Sliding Spectra for AWI-CM3 Output

based on SEP_NEP Spectra Analysis
Current works with regridded output 
(need to do lazy evaluation and better memory management for higher resolution...)

Can also rewrite to loop by # of periods, which may be faster


Created on Mon May  4 09:56:30 2026

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

proj    = ccrs.PlateCarree()
figpath = "/home/niu4/gliu8/figures/bydate/2026-05-05/"
proc.makedir(figpath)
mons3   = proc.get_monstr()

#%% Indicate Experiment and Flux Name

stall = time.time()

# Indicate Experiment and Flux Name
expname                = "TCo319-DART-ssp585d-gibbs-charn"#"TCo319_ssp585"
flxname                = "tscre"#"eis"
regrid                 = True 

# Indicate Sliding Window Options
nsmooth                = 5  # Smoothing over Adjacent Bands
nyr_window             = 24 # Length of Sliding Window

#%%
# Open/Load the files
if regrid:
    rawpath            = "/home/niu4/gliu8/projects/scrap/regrid_1x1/"
    ncname             = "%s_%s_regrid1x1.nc" % (expname,flxname)
else:
    rawpath            = "/home/niu4/gliu8/projects/scrap/processed_global/"
    ncname             = "%s_%s.nc" % (expname,flxname)
dsflx_raw_awi      = xr.open_dataset(rawpath + ncname)[flxname].load()

# Do some preprocessing, get dimension sizes
dsflx_raw_awi      = ut.standardize_names(dsflx_raw_awi)        # Standardize Dimension Names
dsflx_raw_awi      = dsflx_raw_awi.transpose('time','lat','lon')# Transepose
dsflx_raw_awi      = ut.varcheck(dsflx_raw_awi,flxname,expname) # Correct for accumulation period
ntime,nlat,nlon    = dsflx_raw_awi.shape

# ========================
create_arr = 0 # To Ensure that the arrays are created on 1st loop
for o in tqdm(range(nlon)):
    for a in range(nlat):
        
        # Get Point and Subset into Periods
        dspt            = dsflx_raw_awi.isel(lon=o,lat=a)
        if np.any(np.isnan(dspt)):
            #print("Skipping point a=%i,o=%i" % (o,a))
            continue
        
        subsets,tranges = ut.generate_periods(dspt,nyr_window)
        nperiods        = len(subsets) # Get Number of Periods
        tcenters        = [ut.get_center_time(t) for t in tranges]
        
        # Preprocess by Period
        subsets_anom   = ut.preprocess_byperiod(subsets,verbose=False) # Don't Detrend...
        
        
        
        # Calculate Power Spectra for each period
        spec_byperiod  = []
        for nw in range(nperiods):
            sample_in            = subsets_anom[nw].data
            specout              = proc.point_spectra(sample_in,nsmooth=nsmooth,return_conf=True)
            spec_byperiod.append(specout)
            
        # Concatenate By period
        spec_byperiod            = xr.concat(spec_byperiod,dim='period')
        nfreq = len(spec_byperiod.freq)
        
        # Preallocate on first loop (Takes ~5 min)
        
        if create_arr == 0: #(o == 0) and (a == 0): # Make the Arrays
            create_arr += 1
            print("Creating Array...")
            st1 = time.time()
            spectra_all = np.zeros((nperiods,nfreq,nlat,nlon)) * np.nan    # [Period x Freq x Lat x Lon]
            CC_all      = np.zeros((nperiods,nfreq,nlat,nlon,2)) * np.nan  # [Period x Freq x Lat x Lon x Conf Level]
            proc.printtime(st1,print_str="Preallocated in...")
        
        # Save to Array
        spectra_all[:,:,a,o] = spec_byperiod.spectra.data
        CC_all[:,:,a,o,:]    = spec_byperiod.CC.data
        
        
        
        # --- End Lat Loop
    # --- End Lon Loop

# Set Up DataArray Coords
coords_spectra    = dict(period=tcenters,
                      freq=spec_byperiod.freq,
                      lat=dsflx_raw_awi.lat,
                      lon=dsflx_raw_awi.lon)
coords_CC         = coords_spectra.copy()
coords_CC['clvl'] = spec_byperiod.clvl

# Make into Data Arrays and Merge
ds_spectra_all    = xr.DataArray(spectra_all,
                              coords=coords_spectra,
                              dims=coords_spectra,
                              name="spectra",
                              )
ds_CC_all         = xr.DataArray(CC_all,
                                 coords=coords_CC,
                                 dims=coords_CC,
                                 name="CC",
                                 )
ds_out            = xr.merge([ds_spectra_all,ds_CC_all])

# Save Output
outpath = "/home/niu4/gliu8/projects/ccfs/metrics/regrid_1x1/scrap/sliding_spectra/"
outname = "%sSliding_Spectra_%s_%s_winsize%2i_nsmooth%02i.nc" % (outpath,expname,flxname,nyr_window,nsmooth)
edict   = proc.make_encoding_dict(ds_out)
ds_out.to_netcdf(outname,encoding=edict)

proc.printtime(stall,"Completed Script in")

# Completed Script in in 121184.41s (took ~33 hours...)


