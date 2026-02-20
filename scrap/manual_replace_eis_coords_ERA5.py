#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Having some issues with regridding, had to change the variables

Created on Tue Nov 18 10:57:41 2025

@author: gliu

"""

import numpy as np
import xarray as xr

#nc1 = #"/home/niu4/gliu8/projects/common_data/ERA5/anom_detrend1/oldfiles/eis_1979_2024.nc"
# nc1_out ="/home/niu4/gliu8/projects/common_data/ERA5/anom_detrend1/eis_1979_2024.nc"
nc2 = "/home/niu4/gliu8/projects/common_data/ERA5/anom_detrend1/sst_1979_2024.nc"
nc1 = "/home/niu4/gliu8/share/ERA5/processed/oldfiles/eis_1979_2024.nc"

nc1_out = "/home/niu4/gliu8/share/ERA5/processed/eis_1979_2024.nc"

ds1     = xr.open_dataset(nc1).load()

ds2     = xr.open_dataset(nc2)

times   = ds2.valid_time
lon     = ds2.longitude
lat     = ds2.latitude

dsnew   = xr.zeros_like(ds2.sst)


dsmerge = xr.merge([dsnew,ds1.eis.rename(dict(time='valid_time'))])

dsmerge = dsmerge.drop_vars('sst')

dsmerge.to_netcdf(nc1_out)



# dsnew['eis'] = ds1.eis.rename(dict(time='valid_time'))
# #dsnew['eis'] = ds1
# dsnew = dsnew.drop('sst')

