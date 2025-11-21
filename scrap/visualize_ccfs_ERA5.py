#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Visualize Mean ERA5 CCFs

Created on Tue Nov 18 11:24:55 2025

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

#%% Other Functions


def mlr(X,y):
    # MLR fit using scipy 
    
    # Initialize Model and Fit
    model             = LinearRegression()
    model.fit(X,y)
    pred              = model.predict(X)
    # Calculate Error and other variables
    mlr_out           = {}
    mlr_out['pred']   = pred
    mlr_out['err']    = y - pred # Model Error # []
    mlr_out['coeffs'] = model.coef_
    mlr_out['r2']     = sklearn.metrics.r2_score(y,pred)
    
    return mlr_out

def mlr_ccfs(ccfs,flx,standardize=True,fill_value=0,verbose=False):
    # Perform MLR
    #    ccfs: LIST of DataArrays [variable x time]
    #    flx:  DataArray [time x 1]
    
    # Set up Predictors and Target (convert to Numpy Arrays)
    predictors = np.array([ds for ds in ccfs]) # [variable x time]
    if standardize:
        if verbose:
            print("Standardizing each variable")
        predictors = np.array([ds/np.nanstd(ds) for ds in list(predictors)])
    X = predictors.T
    y = flxpt.data
    
    # Replace NaN Values in Predictors
    if verbose:
        if np.any(np.isnan(X)):
            print("NaN values detected! Replace with %f" % fill_value)
    X = np.where(np.isnan(X),fill_value,X) # Set NaN to zero
    
    # Use sklearn for now (can try LSE manual later...)
    mlr_out = mlr(X,y)
    return mlr_out

def reduce_time(ds,dsst):
    dsnew,dsst = proc.match_time_month(ds,dsst)
    return dsnew

# =================
#%% User Selections
# =================

regrid1x1   = False
expname     = "ERA5_1979_2024"
datpath     = "/home/niu4/gliu8/projects/common_data/ERA5/anom_detrend1/"
outpath     = "/home/niu4/gliu8/projects/ccfs/"

meanpath    = "/home/niu4/gliu8/projects/ccfs/time_mean/"

figpath     = "/home/niu4/gliu8/figures/bydate/2025-11-18/"

if regrid1x1:
    datpath = "/home/niu4/gliu8/projects/common_data/ERA5/regrid_1x1/anom_detrend1/"
    outpath = "/home/niu4/gliu8/projects/ccfs/regrid_1x1/"

flxnames    = ['cre',]#'allsky','clearsky']# ['cre',]
ccf_vars    = ["sst","eis","Tadv","r700","w700","ws10",]#"ucc"] 




tstart   = '1979-01-01'
tend     = '2024-12-31'
timename = 'valid_time'
latname  = 'latitude'
lonname  = 'longitude'

#% Load Land Mask
landnc   = datpath + "mask_1979_2024.nc"
landmask = xr.open_dataset(landnc).mask


#%% Load time means

ds_all = []
for ccf in ccf_vars:

    ncstr       = meanpath + "%s_1979_2024.nc" % ccf
    ds = xr.open_dataset(ncstr)[ccf].load()
    ds_all.append(ds)

ds_all = [ut.standardize_names(ds.squeeze()) for ds in ds_all]

#%% Make the Scott et al PLot

import cartopy

def init_globalmap(nrow=1,ncol=1,figsize=(12,8)):
    proj            = ccrs.Robinson(central_longitude=-180)
    bbox            = [-180,180,-90,90]
    fig,ax          = plt.subplots(nrow,ncol,subplot_kw={'projection':proj},figsize=figsize,constrained_layout=True)
    
    multiax = True
    if (type(ax) == mpl.axes._axes.Axes) or (type(ax) == cartopy.mpl.geoaxes.GeoAxes):
        ax = [ax,]
        multiax = False
    
    for a in ax:
        a.coastlines(zorder=10,lw=0.75,transform=proj)
        a.gridlines(ls ='dotted',draw_labels=True)
        
    if multiax is False:
        ax = ax[0]
    return fig,ax

proj   = ccrs.PlateCarree()
varsin = ds_all

vlimsdict = dict(
    sst = [273,303],
    eis = [-8,8],
    Tadv = [-2.5,2.5],
    r700 = [20,70],
    ws10   = [3,13],
    w700 = [-.1,.1],#w700 = [-65,65]   
    )
dtday = 3600*24

for vv in tqdm.tqdm(range(len(varsin))):
    
    fig,ax      = init_globalmap(figsize=(8,3.5))
    
    plotvar     = varsin[vv]#.mean('time')
    if plotvar.name == "Tadv":
        plotvar = plotvar * dtday
        
    # if plotvar.name == "sst" and :
    #     plotvar = plotvar + 273.15
    
    vlims       = vlimsdict[plotvar.name]
    pcm         = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,cmap='cmo.balance',vmin=vlims[0],vmax=vlims[1])
    cb          = fig.colorbar(pcm,ax=ax,pad=0.01,fraction=0.025)
    

    ax.set_title(plotvar.name)
    savename    = "%s%s_%s_CCFs_TimeMean_check.png" % (figpath,"ERA5",plotvar.name)
    plt.savefig(savename,dpi=150,bbox_inches='tight')
    #plt.show()
