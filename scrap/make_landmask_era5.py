#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Make Landmask for ERA5 from SST

Created on Tue Nov 18 10:18:18 2025

@author: gliu

"""

import xarray as xr
import numpy as np

datpath = "/Users/gliu/Downloads/02_Research/01_Projects/05_SMIO/01_Data/"
ncfile  = "sst_1979_2024.nc"
ds      = xr.open_dataset(datpath+ncfile)

dsmask  = xr.where(np.isnan(ds.isel(valid_time=1).sst),np.nan,1)
dsmask  = dsmask.rename('mask')

dsmask.to_netcdf(datpath + "era5_landmask_fromsst.nc")
