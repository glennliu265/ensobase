#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Scrap to check the skew calculations using cdo [compute_skew.sh]

Created on Tue Jul 14 13:33:38 2026

@author: gliu

"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import scipy as sp

def dailyclim(ds):
    return ds.groupby("time.dayofyear").mean("time")

def deseason_daily(ds,clim=False):
    dsclim       = dailyclim(ds)#ds.groupby("time.dayofyear").mean("time")
    dsanom       = ds.groupby("time.dayofyear") - dsclim
    if clim:
        return dsanom,dsclim
    return dsanom

def loadcdo(filename):
    ds = xr.open_dataset(dpath+"/"+filename)
    return ds.TS.isel(ncol=imin).load()


dpath = '/home/niu4/gliu8/projects/mesaclip/region_crops'

dsref = xr.open_dataset(dpath+"/HiRes_BHISTC5_TS_Pacific_ens001.nc")

lon   = dsref.lon.load()
lat   = dsref.lat.load()

skewnc  = "%s/skewtest2.nc" % dpath
cdoskew = xr.open_dataset(skewnc).load()

# Check Skew Computed Via CDO....
plt.scatter(lon,lat,c=cdoskew.TS.squeeze())
plt.colorbar()
plt.show()


# Load Raw Timeseries from File
lonf = 260
latf = -10
imin = np.argmin( (lon.data-lonf)**2 + (lat.data-latf)**2).item()
print(lon.isel(ncol=imin))
print(lat.isel(ncol=imin))

dsptraw = dsref.isel(ncol=imin).load().TS

# Part (1) Check Preprocessing (Ok) -------------------------------------------
dsptanom    = deseason_daily(dsptraw)
dsptanom_dt = xr.ones_like(dsptanom) * sp.signal.detrend(dsptanom.data,axis=0,type='linear')

# Load CDO-Processed Version
#ncanom  = "%s/HiRes_BHISTC5_TS_Pacific_ens001_anom_detrended.nc" % dpath
#cdoanom = xr.open_dataset(ncanom).TS.isel(ncol=imin).load()
cdoanom = loadcdo("HiRes_BHISTC5_TS_Pacific_ens001_anom_detrended.nc")

# Check
fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(12.5,4.5))

plotvar = cdoanom
ax.plot(plotvar.time,plotvar,label="cdo",c='gray')

plotvar = dsptanom_dt
ax.plot(plotvar.time,plotvar,label="cdo",c='red',lw=.75)

ax.legend()
plt.show()

# Check Max (it was 0.00012493134) 
print(np.nanmax(np.abs(cdoanom-dsptanom_dt)))


# Part (2) Check Numerator Calculation (Ok)

# Close enough?
mu    = np.nanmean(dsptanom_dt) # np.float32(-7.777656e-09)
cdomu = loadcdo('mu.nc')        # array([3.7999437e-11], dtype=float32)

# Pretty much there
numer    = ((dsptanom_dt - mu)**3).sum('time') # array(4985.423, dtype=float32)
cdonumer = loadcdo('skew_numer.nc')         # array([4985.424], dtype=float32)

# Part (3) Check Denominator Calculation (Ok)

# Check Denominator
denom   = dsptanom_dt.std('time')**3  # array(0.24899009, dtype=float32)
cdodenom = loadcdo('skew_denom.nc')   # array([0.24899009], dtype=float32)


# Part (4) Finally Check Skew
skew    = 1/len(dsptanom_dt.time) * numer/denom
cdoskew = loadcdo('skewtest2.nc')




