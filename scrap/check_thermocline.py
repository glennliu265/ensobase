#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Check Thermocline Depth in AWI Simulations

Created on Thu Sep 11 14:05:30 2025

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


#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

#%% Indicate paths

figpath         = "/home/niu4/gliu8/figures/bydate/2025-09-16/"
proc.makedir(figpath)

datpath         = "/home/niu4/gliu8/projects/scrap/TP_crop/"

expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090"]

vname           = "temp_1m"


nexps = len(expnames)
#%% Load the thermocline data

dsall = []
for ex in tqdm.tqdm(range(nexps)):
    ncname = "%s%s_%s.nc" % (datpath,expnames[ex],vname)
    ds  = xr.open_dataset(ncname).load()
    dsall.append(ds)

#%% Find the Maximum Gradient in temperature

# Find the vertical gradient
dTdz_all = [ds.temp.differentiate('nz1') for ds in dsall]

# Locate the depth of maximum
Dmaxgrad = [np.abs(ds).idxmax('nz1') for ds in dTdz_all]

# More Debugging
# test = dTdz_all[-1].isel(time=0,lat=0)
# testmax = np.abs(test).idxmax('nz1')

# fig,ax = plt.subplots(1,1)
# test.plot(ax=ax)
# testmax.plot(ax=ax)
# plt.show()

# Tried to make it into function before I discovered idxmax....
#dTdz = dTdz_all[-1].temp
# def get_zmax(dTdz):
      # Apparently don't need a function because it is all handleed by this
#     # First, replace all-NaN slices with zero
#     dzfilt  = xr.where(np.isnan(dTdz),-99,dTdz)
#     #nanmask = xr.where(np.isnan(dTdz.sum('nz1',skipna=False)),np.nan,1)
#     # Find Argument of Max
#     #zlvlmax = dzfilt.idxmax('nz1')
#     return zargmax

#%% Find the 20C Isoterm

def find_D20(ds):
    # Detect if in Kelvin (rather coarse approach)
    if np.any(ds > 270):
        ds = ds - 273.15
    diff20 = ds - 20
    return np.abs(diff20).idxmin('nz1',skipna=True)


D20 = [find_D20(ds.temp) for ds in dsall]

# Debugging (turns out I forgot to add a maximum)
# test     = dsall[-1].temp.isel(time=0,lat=0)
# testdiff = test - 20
# zmin     = np.abs(testdiff).idxmin("nz1")

# fig,ax = plt.subplots(1,1)
# testdiff.plot(ax=ax)
# ax.contour(testdiff.lon,testdiff.nz1,testdiff,levels=[0,])
# zmin.plot(ax=ax)
# plt.show()


#%% Save the output for some visualization

for vv in range(2):
    
    if vv == 0:
        dsin  = Dmaxgrad
        dname = "Dmaxgrad"
    else:
        dsin = D20
        dname = "D20"
    
    for ex in tqdm.tqdm(range(4)):
        
        ncout = "%s%s_%s.nc" % (datpath,expnames[ex],dname)
        edict = proc.make_encoding_dict(dsin[ex])
        dsin[ex].to_netcdf(ncout,encoding=edict)

#%% Try to look at the meridional mean (or a latitude slice)

# Meridional Mean
# D20_ymean  = [ds.mean('lat') for ds in D20]
# Dmax_ymean = [ds.mean('lat') for ds in Dmaxgrad]
# temp_ymean = [ds.temp.mean('lat') for ds in dsall]

# Latitude Slice
ilat       = 20
D20_ymean  = [ds.isel(lat=ilat) for ds in D20]
Dmax_ymean = [ds.isel(lat=ilat) for ds in Dmaxgrad]
temp_ymean = [ds.isel(lat=ilat) for ds in dsall]

#%% Make some plots of Longitude vs Depth

cints = np.arange(0,39,3)
ex    = -4

itime = 0

fig,ax = plt.subplots(1,1,figsize=(12.5,6),constrained_layout=True)

plotvar = temp_ymean[ex].isel(time=itime).temp
#pcm = ax.pcolormesh(plotvar.lon,plotvar.nz1,plotvar)
pcm = ax.contourf(plotvar.lon,plotvar.nz1,plotvar,levels=cints,cmap='cmo.thermal')

cl = ax.contour(plotvar.lon,plotvar.nz1,plotvar,levels=[20,],colors='gray')#,ls='dotted',linewidths=.75)

ax.clabel(cl)

plotline = D20_ymean[ex].isel(time=itime)
ax.plot(plotline.lon,plotline,color="cyan",label="$D_{20}$")

plotline = Dmax_ymean[ex].isel(time=itime)
ax.plot(plotline.lon,plotline,color="k",label="$D_{GradMax}$")

ax.set_title("%s Zonal Section @ Lat = %s, Time=%s" % (expnames_long[ex],plotvar.lat.data,plotvar.time.data))

ax.set_ylim([0,1000])
ax.invert_yaxis()

fig.colorbar(pcm,ax=ax)
ax.legend()
plt.show()

#%% Take Seasonal Mean and Examine the Pattern

D20_scycle = [ds.groupby('time.season').mean('time') for ds in D20]
Dmax_scycle = [ds.groupby('time.season').mean('time') for ds in Dmaxgrad]

#%% Check The Seasonal Differences in Pattern

ex  = -2
sid = 0

for ex in range(4):
    for sid in range(4):
        proj   = ccrs.PlateCarree(central_longitude=180)
        projd  = ccrs.PlateCarree()
        def init_tp_map():
            bbplot = [120, 290, -20, 20]
            proj   = ccrs.PlateCarree(central_longitude=180)
            projd  = ccrs.PlateCarree()
            fig,ax = plt.subplots(1,1,figsize=(12.5,4.5),subplot_kw={'projection':proj})
            #ax     = viz.add_coast_grid(ax,bbox=bbplot,fill_color='k',proj=ccrs.PlateCarree())
            return fig,ax
        
        
        fig,ax  = init_tp_map()
        ax.coastlines()
        plotvar = D20_scycle[ex].isel(season=sid) - Dmax_scycle[ex].isel(season=sid)
        pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,vmin=-50,vmax=50,
                                cmap='cmo.balance',transform=projd)
        
        fig.colorbar(pcm,ax=ax,fraction=0.010)
        ax.set_title("D20 - Dmax (%s, %s) [meters]" % (plotvar.season.data,expnames_long[ex]))
        
        
        figname = "%sD20_v_DMax_%s_%s.png" % (figpath,expnames[ex],plotvar.season.data.item())
        plt.savefig(figname,dpi=150,bbox_inches='tight')
        print(figname)

#%% Other Scrap below



#%% Debugging for the Gradient CAlculation (checking centrede difference )
    




# test        = dsall[-1]
# testfake    = test.copy()
# z           = test.nz1

# # Differentiate
# #test['fakez'] = np.arange(len(z))
# testdz      = test.differentiate('nz1')

# # Make fake and differentiate as well
# testfake['nz1'] =  np.arange(len(z)) #= test.differentiate('fakez')
# testdzfake  = testfake.differentiate('nz1')


"""

.isel(time=0,lat=22,lon=22)
array([32.877617 , 32.815716 , 32.777912 , 32.71042  , 32.502182 ,
       32.35916  , 32.215977 , 32.095222 , 32.001102 , 31.881626 ,
       31.622875 , 30.803707 , 27.862547 , 23.577906 , 19.758133 ,
       16.70301  , 15.07375  , 13.778761 , 12.944878 , 12.581959 ,
       12.312361 , 11.99538  , 11.666052 , 11.2631855, 10.938492 ,
       10.637009 , 10.339666 , 10.070934 ,  9.737741 ,  9.4206505,
        9.0543785,  8.709652 ,  8.357594 ,  7.955512 ,  7.570633 ,
        7.268164 ,  7.0460424,  6.7783823,  6.583274 ,  6.283485 ,
        5.961723 ,  5.607806 ,  5.3441696,  4.9998035,  4.5800705,
        4.0821886,  3.6650605,  3.24057  ,  2.9211931,  2.6663957,
        2.4386392,  2.2417204,  2.076672 ,  1.909014 ,  1.7753792,
        1.6697878,  1.562441 ,  1.4584618,  1.3741848,
        
array([-0.01238022, -0.00944433, -0.00577285, -0.01378651, -0.01756287,
       -0.01431026, -0.01319695, -0.01074371, -0.01067982, -0.01891136,
       -0.04901544, -0.14974867, -0.26866686, -0.28123477, -0.22916317,
       -0.15614611, -0.09747498, -0.0709624 , -0.0398934 , -0.02013813,
       -0.01482916, -0.01365465, -0.0146439 , -0.0145512 , -0.01252354,
       -0.01197651, -0.0108834 , -0.01041024, -0.01083806, -0.01138938,
       -0.01184998, -0.01161308, -0.01162117, -0.01061721, -0.00859185,
       -0.0063261 , -0.00513364, -0.00462769, -0.00494897, -0.00621551,
       -0.00675679, -0.00617553, -0.00608003, -0.00681106, -0.00668059,
       -0.00610007, -0.00561079, -0.00495912, -0.00382783, -0.00303067,
       -0.00230718, -0.00180984, -0.00166353, -0.00150646, -0.00119613,
       -0.00106469, -0.00105663, -0.00087731,         nan,         nan,
               nan,         nan,         nan,         nan,         nan,
               nan,         nan,         nan,         nan,         nan,
               nan,         nan,         nan,         nan,         nan,
               nan,         nan,         nan,         nan], dtype=float32)

"""

# (32.777912 - 32.877617) / (15-2.5)

# #%% Plot test

# fig,axs= plt.subplots(2,1)

# ax = axs[0]
# ax.plot(testdzfake.isel(time=0,lat=22,lon=22).temp.data,label="Fake")
# ax.plot(testdz.isel(time=0,lat=22,lon=22).temp.data,label="differentiate")
# ax.legend()

# ax = axs[1]
# ax.plot(z.data)

# plt.show()





