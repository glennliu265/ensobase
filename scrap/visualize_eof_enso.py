#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Visualize the EOF Patterns computed via cdo

Created on Wed Nov 19 15:05:14 2025

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

#%% User Edits

datpath         = "/home/niu4/gliu8/projects/ccfs/enso_eof/"
expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C","ERA5_1979_2024"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950","ERA5 (1979 to 2024)"]

searchstr       = datpath + "%s_ENSO_EOF.nc"


figpath         = "/home/niu4/gliu8/figures/bydate/2025-11-25/"
proc.makedir(figpath)

#%%
nexps = len(expnames)
ds_all = []
for ex in range(nexps):
    
    ncname = searchstr % expnames[ex]
    ds = xr.open_dataset(ncname).load() 
    
    ds_all.append(ds)

#%% Plot Variance Explained

xtks = np.arange(0,7)
ytks = np.arange(0,0.75,0.05)

fig,ax = plt.subplots(1,1,figsize=(6,4.5))

for ex in range(nexps):
    plotvar = ds_all[ex].varexp
    ax.plot(plotvar.pc,plotvar,label=expnames_long[ex],marker="o")

ax.set_ylabel("% Variance Explained")
ax.set_xlabel("Principle Component")

ax.set_xticks(xtks)
ax.set_yticks(ytks)
ax.set_ylim([ytks[0],ytks[-1]])
ax.set_xlim([xtks[0],xtks[-1]])
ax.legend()


plt.show()


#%% Plot Variance Explained by first 2 modes




#%% Plot ENSO Pattern associated with each

proj        = ccrs.PlateCarree()
#bbox_narrow = [130,360-80,-3,3]

ex          = -1
N_mode      = 0
vmax        = 1.5

for ex in range(nexps):
    for N_mode in range(2):
        fig,ax  = ut.init_tp_map(1,1)
        
        plotvar = ds_all[ex].eofs.isel(pc=N_mode)
        pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,
                                cmap='cmo.balance',vmin=-vmax,vmax=vmax)
        #ax.set_extent(bbox_narrow)
        
        ax.set_ylim([-10,10])
        title   = "%s, Mode %i (Variance Explained = %.2f%%)" % (expnames_long[ex],N_mode+1,ds_all[ex].varexp.isel(pc=N_mode)*100)
        ax.set_title(title)
        
        
        figname = "%s%s_EOF_PC%02i.png" % (figpath,expnames[ex],N_mode+1)
        plt.savefig(figname,dpi=150,bbox_inches='tight')
        #plt.show()

#%% Try to compute rotated EOFs with PCs

rotated_pc = []
for ex in range(nexps):
    
    pc1 = ds_all[ex].pcs.sel(pc=1)
    pc2 = ds_all[ex].pcs.sel(pc=2)
    
    # Compute EP
    ep  = (pc1 - pc2) / np.sqrt(2)
    
    # Compute CP
    cp  = (pc1 + pc2) / np.sqrt(2)
    
    out = xr.merge([ep.rename('ep'),cp.rename('cp')])
    rotated_pc.append(out)

#%% Save the output
datpath = "/home/niu4/gliu8/projects/scrap/nino34/"
for ex in range(nexps):
    
    ncname = "%s%s_enso_eof_rotated.nc" % (datpath,expnames[ex])
    rotated_pc[ex].to_netcdf(ncname)
    print(ncname)
    

    
#%% Sample Plot of Temperature around Hawaii

eofsin  = [ds_all[ex].eofs.isel(pc=N_mode) for ex in range(nexps)]

bbsel   = [360-160,360-154,18.5,23]
dsreg   = [proc.sel_region_xr(ds.squeeze(),bbsel) for ds in eofsin]


msktest = [ut.load_land_mask_awi(expname) for expname in expnames]
maskreg = [proc.sel_region_xr(ds.squeeze(),bbsel) for ds in msktest]
proj    = ccrs.PlateCarree()


#%%
fig,axs  = plt.subplots(1,3,constrained_layout=True,subplot_kw={'projection':proj},figsize=(16,4))

for ii in range(3):
    
    expids = [0,1,3]
    ex = expids[ii]
    # Get Axis
    ax = axs.flatten()[ii]
    ax = viz.add_coast_grid(ax,bbsel)
    
    plotvar = dsreg[ex]
    pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,vmin=-1,vmax=1,cmap='cmo.balance')


    # Plot Mask
    plotvar = xr.where(np.isnan(maskreg[ex]),0,1).squeeze()
    lon = plotvar.lon
    lat = plotvar.lat
    viz.plot_mask(lon,lat,plotvar.T,marker="o",geoaxes=True,markersize=0.1,ax=ax,proj=proj,reverse=True) 


plt.show()
