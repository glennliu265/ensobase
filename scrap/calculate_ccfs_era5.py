#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate CCFs using ERA5 datasets

 - Wind Speed
 - Temperature Advection (SST)

Created on Mon Nov 17 15:04:10 2025

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

from sklearn.linear_model import LinearRegression
import sklearn

#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% Additional Functions

def calc_grad_centered(ds,latname='lat',lonname='lon'): # Copied from structure in calc_ekman_advection_htr
    
    if latname != 'lat' or lonname != 'lon':
        lldict          = {latname : 'lat', lonname : 'lon'}
        lldict_reverse  = {'lat' : latname, 'lon' : lonname}
        rename_flag=True
        ds = ds.rename(lldict)
    else:
        rename_flag=False
        
    
    dx,dy = proc.calc_dx_dy(ds.lon.values,ds.lat.values,centered=True)
    
    # Convert to DataArray
    daconv   = [dx,dy]
    llcoords = {'lat':ds.lat.values,'lon':ds.lon.values,}
    da_out   = [xr.DataArray(ingrad,coords=llcoords,dims=llcoords) for ingrad in daconv]
    dx,dy = da_out
    
    # Roll and Compute Gradients (centered difference)
    ddx = (ds.roll(lon=-1) - ds.roll(lon=1)) / dx
    ddy = (ds.roll(lat=-1) - ds.roll(lat=1)) / dy
    ddy.loc[dict(lat=ddy.lat.values[-1])] = 0 # Set top latitude to zero (since latitude is not periodic)
    
    if rename_flag:
        ddx,ddy =  [ds.rename(lldict_reverse) for ds in [ddx,ddy]]
    
    return ddx,ddy



#%% Load the Datasets

datpath="/home/niu4/gliu8/share/ERA5/processed/"

vnames = ["u10","v10","sst"]

ds_all = []
for vv in tqdm.tqdm(range(len(vnames))):
    
    ncname  = "%s%s_1979_2024.nc" % (datpath,vnames[vv])
    ds      = xr.open_dataset(ncname)[vnames[vv]].load()
    ds_all.append(ds)
    
# ==============================    
#%% Part 1 (compute wind speed)
# ==============================   

u10         = ds_all[0]
v10         = ds_all[1]
dsws        = np.sqrt(u10**2 + v10**2)#**(1/2)

# def calc_ws(ds):
#     return ((ds[0])**2 + (ds[1])**2)**0.5

#dsws = calc_ws(ds_all)

#%% Check the mean pattern... (again it appears weird)

wsmean  = dsws.mean('valid_time') # .sel(valid_time=slice('1979-01-01','2010-01-01'))
wsmean.plot(vmin=3,vmax=13,cmap='cmo.balance'),plt.show()

#dsws = [calc_ws(ds) for ds in ds_all]

#%% Try loading explicitly downlaoded wind speed

ncw10 = "ws10_1979_2024.nc"
dsw10 = xr.open_dataset(datpath+ncw10)['ws10'].load()

ws10mean = dsw10.mean('valid_time')
ws10mean.plot(vmin=3,vmax=13,cmap='cmo.balance'),plt.show()
# This is correct, but somehow my wind speed computations are... wrong?


#%% Part 2 (compute temperature advection)

u10       = ds_all[0]
v10       = ds_all[1]
sst       = ds_all[-1]
RE        = 6375e3 # Earth Radius, in meters

ddx,ddy   = calc_grad_centered(sst,latname='latitude',lonname='longitude')

Tadv2     = - u10 * ddx - v10 * ddy
Tadv2mean = Tadv2.mean('valid_time')
dtday     = 3600*24

(Tadv2mean*dtday).plot(vmin=-2.5,vmax=2.5,cmap='cmo.balance'),plt.show()
 
#%% Save the output (Tadv)

Tadv2 = Tadv2.rename("Tadv")

ncname = "%sTadv_1979_2024.nc" % (datpath)
edict  = proc.make_encoding_dict(Tadv2)
Tadv2.to_netcdf(ncname)



# #%%

# xrddx   = sst.differentiate('longitude')
# xrddy   = sst.differentiate('latitude')

# fig,axs = plt.subplots(2,1)
# (xrddx.isel(valid_time=0)/RE).plot(ax=axs[0])
# ddx.isel(valid_time=0).plot(ax=axs[1]),plt.show()

#test = (xrddx.isel(valid_time=0)/ddx.isel(valid_time=0) )

#%% Old Scrap Below, where xr.differentiate was not working....
# Still need to figure out in the future
#
# lat = sst.latitude
# lon = sst.longitude
# xx,yy = np.meshgrid(lon,lat)
# phi       = np.radians(yy)#np.radians(sst.latitude)
# lbd       = np.radians(xx)#np.radians(sst.longitude)

# st        = time.time()
# dSST_dlbd = sst.differentiate('longitude')
# dSST_dphi = sst.differentiate('latitude')
# print("\tDifferentiated SST along Lon/Lat in %.2fs" % (time.time()-st))

# st        = time.time()
# Term1     = u10/np.cos(phi) * dSST_dlbd
# Term2     = v10 * dSST_dphi
# print("Calculated Term in %.2fs" % (time.time()-st))


# Tadv      = (-Term1 - Term2) / RE

# #%%
# Tadvmean = Tadv.mean('valid_time')

# dtday    = 3600*24
# (Tadvmean*dtday).plot(vmin=-2.5,vmax=2.5,cmap='cmo.balance'),plt.show()


# #%% Do a test for temperature advection






# lontest = np.arange(0,370,10)
# lattest = np.arange(-90,100,10)

# xt,yt = np.meshgrid(lontest,lattest)

# coords = dict(lat=lattest,lon=lontest,)

# lont = xr.DataArray(xt,coords=coords,dims=coords,name='xx')
# latt = xr.DataArray(yt,coords=coords,dims=coords,name='yy')
# #dummyvar = 


# #%% Try method borrwed from calc_geostrophic advection


# lldict_reverse  = dict(lon='longitude',lat='latitude')
# lldict          = dict(longitude='lon',latitude='lat')
                      
                      
# sstnew = sst.rename(dict(longitude='lon',latitude='lat'))




# ddx,ddy = calc_grad_centered(sstnew)

# ddx,ddy = [ds.rename(lldict_reverse) for ds in [ddx,ddy]]

# #%% Recompaute gradient

# Tadv2     = - u10 * ddx - v10 * ddy
# Tadv2mean = Tadv2.mean('valid_time')
# (Tadv2mean*dtday).plot(vmin=-2.5,vmax=2.5,cmap='cmo.balance'),plt.show()

# #test = ddx.isel(valid_time=0)



