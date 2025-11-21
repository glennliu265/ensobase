#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Investigate the change in variance between simulations
(I needed a break)

Created on Thu Nov 13 13:57:48 2025

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

#%% User Edits


figpath             = "/home/niu4/gliu8/figures/bydate/2025-11-18/"

expnames            = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
expnames_long       = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"]
datpath             = "/home/niu4/gliu8/projects/scrap/summary_stats/"
nexps               = len(expnames)

ds_all = []
for ex in tqdm.tqdm(range(nexps)):
    ncsearch = "%s%s_sst_*.nc" % (datpath,expnames[ex])
    nclist   = glob.glob(ncsearch)[0]
    ds = xr.open_dataset(nclist).load()
    ds_all.append(ds)
    

#%% Check the monthly variance for a simnulation

ex      = 0

proj    = ccrs.PlateCarree()
plotim  = np.roll(np.arange(12),1)
mons3   = proc.get_monstr()



vmin    = 0
vmax    = 3


for ex in tqdm.tqdm(range(nexps)):
    
    fig,axs = ut.init_globalmap(4,3,figsize=(18,18))
    
    for a in range(12):
        im = plotim[a]
        ax = axs[a]
        ax.set_title(mons3[im],fontsize=22)
        
        # # ---
        
        plotvar = ds_all[ex].monvar.isel(month=im)
        pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,
                                vmin=vmin,vmax=vmax,cmap='cmo.thermal'
                                )
        # # ---
    
    cb = viz.hcbar(pcm,ax=axs.flatten(),fig=fig)
    cb.set_label(r"SST Variance ($\degree C^2$)",fontsize=22)
    
    
    figname = "%s%s_sst_monthly_variance.png" % (figpath,expnames[ex],)
    plt.savefig(figname,dpi=150,bbox_inches='tight')
    #plt.show()

#%% Calculate some Differences

diff31km = ds_all[1].monvar - ds_all[0].monvar
diff9km  = ds_all[3].monvar - ds_all[2].monvar

#%%
plotdiffs = [diff31km,diff9km]
diffnames = ["31km","9km"]
vmin      = -2
vmax      = 2

for ii in range(2):
    
    fig,axs = ut.init_globalmap(4,3,figsize=(18,18))
    
    for a in range(12):
        im = plotim[a]
        ax = axs[a]
        ax.set_title(mons3[im],fontsize=22)
        
        # # ---
        
        plotvar = plotdiffs[ii].isel(month=im)
        pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,
                                vmin=vmin,vmax=vmax,cmap='cmo.balance'
                                )
        
        # # ---
    
    cb = viz.hcbar(pcm,ax=axs.flatten(),fig=fig)
    cb.set_label(r"SST Variance ($\degree C^2$)",fontsize=22)
    
    figname = "%s%s_sst_monthly_variance_difference.png" % (figpath,diffnames[ii],)
    plt.savefig(figname,dpi=150,bbox_inches='tight')
    
    




     
    
    

