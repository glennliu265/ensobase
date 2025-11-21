#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Check global lag regression output

Created on Mon Oct 27 09:55:44 2025

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

expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"]

ensoid_name     = "nino34"
standardize     = False

timecrops       = [[1993,2024]]

datpath         = "/home/niu4/gliu8/projects/scrap/global_anom_detrend2/"
outpath         = "/home/niu4/gliu8/projects/scrap/global_lag_regressions/"
figpath         = "/home/niu4/gliu8/figures/bydate/2025-10-28/"
proc.makedir(figpath)

vnames          = ['tsr','ttr','ttrc','tsrc',]#'ttcre','tscre','tsrc']#str','ssr','skt','ssh','lcc','tcc','ttr','ttrc','tsr','tsrc'] # 'sst',#

#%% Load ENSO ID

ensoids         = [ut.load_ensoid(expname,ensoid_name,standardize=standardize) for expname in expnames]

#%% Looping for each experiment

# Indicate Lead Lags
leadlags        = np.arange(-6,7,1)
sep_mon         = False

ex              = 0
nexps           = len(expnames)
expname         = expnames[ex]
vname           = vnames[0]


datpath_temp    = "/home/niu4/gliu8/projects/scrap/global_lag_regressions/mystery/"

#%%

ds_all = []

for vname in tqdm.tqdm(vnames):
    
    ncname = "%sLagRegression_AllMonths_%s_%s_%s_standardize%i_lag%ito%i.nc" % (datpath_temp,
                                                                                     expnames[ex],vname,ensoid_name,
                                                                                     standardize,leadlags[0],leadlags[-1])
    ds_all.append(xr.open_dataset(ncname).load())

#%% Try Plotting 1 Variable


def init_globalmap(nrow=1,ncol=1,figsize=(12,8)):
    proj            = ccrs.Robinson(central_longitude=-180)
    bbox            = [-180,180,-90,90]
    fig,ax          = plt.subplots(nrow,ncol,subplot_kw={'projection':proj},figsize=figsize,constrained_layout=True)
    ax.coastlines(zorder=10,lw=0.75,transform=proj)
    ax.gridlines(ls ='dotted',draw_labels=True)
    return fig,ax




#test = 


# nrow = 1
# ncol = 1
# figsize         = (12,8)
# projr           = ccrs.Robinson(central_longitude=-180)
# bbox            = [-180,180,-90,90]
# fig,ax          = plt.subplots(nrow,ncol,subplot_kw={'projection':projr},figsize=figsize,constrained_layout=True)


proj    = ccrs.PlateCarree()#central_longitude=-180) # Just set one central longitude to 180... (the one int he fucntion)
v       = 0
vmax    = -20

for v in np.arange(0,len(vnames)):
    
    for it in range(len(leadlags)):
        lag     = leadlags[it]
        
        fig,ax  =  init_globalmap()
        
        plotvar = ds_all[v][vnames[v]].sel(lag=lag)
        plotvar = ut.varcheck(plotvar,vnames[v],expname)
        plotvar = ut.standardize_names(plotvar)
        pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,vmin=-vmax,vmax=vmax,cmap='cmo.balance',zorder=3,transform=proj)
        
        # ax.coastlines(zorder=10,transform=projr)
        # ax.gridlines(ls ='dotted',draw_labels=True)
        
        viz.hcbar(pcm,ax=ax)
        ax.set_title("%s-Nino3.4 Regression, Lag %02i" % (vnames[v],lag))
        
        figname = "%sTest_Lag_regression_%s_it%02i_Lag%02i" % (figpath,vnames[v],it,lag)
        plt.savefig(figname,dpi=150,bbox_inches='tight')
        #plt.show()


#%% Check if I have actually removed the seasonal cycle :(...


datpath     = "/home/niu4/gliu8/projects/scrap/global_anom_detrend2/"
ncname      = datpath+"TCo319_ctl1950d_sst_anom.nc"
ds          = xr.open_dataset(ncname)['sst']

lonf        = 330
latf        = 50

dspt        = proc.selpt_ds(ds,lonf,latf).load()
dsscycle    = dspt.groupby('time.month').mean('time')
dspta       = proc.xrdeseason(dspt)

dsscycle.plot(),plt.show()






