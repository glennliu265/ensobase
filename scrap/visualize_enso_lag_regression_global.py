#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Visualize output from [lag_regression_enso_global]

Created on Mon Oct 27 15:39:25 2025

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

#%% Additional Functions

def init_globalmap(nrow=1,ncol=1,figsize=(12,8)):
    proj            = ccrs.Robinson(central_longitude=-180)
    bbox            = [-180,180,-90,90]
    fig,ax          = plt.subplots(nrow,ncol,subplot_kw={'projection':proj},figsize=figsize,constrained_layout=True)
    
    multiax = True
    if type (ax) == mpl.axes._axes.Axes:
        ax = [ax,]
        multiax = False
    
    for a in ax:
        a.coastlines(zorder=10,lw=0.75,transform=proj)
        a.gridlines(ls ='dotted',draw_labels=True)
        
    if multiax is False:
        ax = ax[0]
    return fig,ax

#%% User Edits

expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"]
expcols = ["cornflowerblue",'lightcoral',
           "slateblue","firebrick",
           "midnightblue","k"] # Includes Glorys and different shade based on resolution
 

ensoid_name     = "nino34"
standardize     = False

timecrops       = [[1993,2024]]

datpath         = "/home/niu4/gliu8/projects/scrap/global_anom_detrend2/"
outpath         = "/home/niu4/gliu8/projects/scrap/global_lag_regressions/"
figpath         = "/home/niu4/gliu8/figures/bydate/2025-10-28/"
proc.makedir(figpath)


vnames          = ['tsr','ttr','ttrc','tsrc','ttcre','tscre','tsrc','allsky','clearsky','cre']#str','ssr','skt','ssh','lcc','tcc','ttr','ttrc','tsr','tsrc'] # 'sst',#

regrid_1x1      = True
leadlags        = np.arange(-18,19,1)
sep_mon         = False
standardize     = False
ensoid_name     = "nino34"

if regrid_1x1:
    datpath     = "/home/niu4/gliu8/projects/scrap/regrid_1x1/global_lag_regressions/"
else:
    datpath     = "/home/niu4/gliu8/projects/scrap/global_lag_regressions/"

figpath         = "/home/niu4/gliu8/figures/bydate/2025-10-28/"
proc.makedir(figpath)


#%% Load all variables

expids      = [0,2,4]
ds_byvar    = []

    
ds_byexp = []
for ii in range(3):
    ex      = expids[ii]
    
    ds_all = []
    for vname in vnames:
        ncname = "%sLagRegression_AllMonths_%s_%s_%s_standardize%i_lag%ito%i.nc" % (datpath,
                                                                                    expnames[ex],vname,ensoid_name,
                                                                                    standardize,leadlags[0],leadlags[-1])
        
        ds = xr.open_dataset(ncname)[vname].load()
        
        ds = ut.standardize_names(ds)
        ds = ut.varcheck(ds,vname,expnames[ex])
        
        ds_all.append(ds)
    
    ds_all = xr.merge(ds_all)
    ds_byexp.append(ds_all)
       
        

#%% Old Loop by Variable


for vname in ['tsr','ttr','ttrc','tsrc','tscre','tsrc','ttcre',]:

    
    ds_all = []
    for ii in range(3):
        ex      = expids[ii]
        ncname = "%sLagRegression_AllMonths_%s_%s_%s_standardize%i_lag%ito%i.nc" % (datpath,
                                                                                    expnames[ex],vname,ensoid_name,
                                                                                    standardize,leadlags[0],leadlags[-1])
        
        ds = xr.open_dataset(ncname)[vname].load()
        ds_all.append(ds)
        
    
    #%% Quick Comparison Plot (All 3 Simulations)
    

    
    proj    = ccrs.PlateCarree()
    vmax    = 10
    lag     = 0
    nlags   = len(leadlags)
    
    for it in range(nlags):
        lag = leadlags[it]
        
        fig,axs = init_globalmap(3,1,figsize=(12,12))
        
        for ii in range(3):
            ax      = axs[ii]
            ex      = expids[ii]
            
            plotvar = ds_all[ii].sel(lag=lag)
            
            plotvar = ut.varcheck(plotvar,vname,expnames[ex])
            plotvar = ut.standardize_names(plotvar)
            pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,vmin=-vmax,vmax=vmax,cmap='cmo.balance',zorder=3,transform=proj)
            
        #plt.suptitle("%s, Lag %02i" % (vname,lag))
        cb = viz.hcbar(pcm,ax=axs.flatten(),fraction=0.025)
        cb.set_label("%s [W/m2/degC], Lag %02i" % (vname,lag))
        #plt.show()
        
        savename = "%sLagRegression_AllMonths_regrid_%s_it%02i_lag%02i.png" % (figpath,vname,it,lag)
        
        
        plt.savefig(savename,dpi=150,bbox_inches='tight')
        
        #plt.show()
        plt.close()
        
    
    
    #%% Look at difference with control
    ref_ii = 0 
    
    vmax = 20
    lag  = 0
    # for it in range(nlags):
    #     lag = leadlags[it]
        
    fig,axs = init_globalmap(2,1,figsize=(12,12))
    
    for ii in range(2):
        ax      = axs[ii]
        ex      = expids[ii]
        ii_sel  = ii+1 # Offset to skip first
        
        plotvar = ds_all[ii_sel].sel(lag=lag) - ds_all[ref_ii].sel(lag=lag)
        
        plotvar = ut.varcheck(plotvar,vname,expnames[ex])
        plotvar = ut.standardize_names(plotvar)
        pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,vmin=-vmax,vmax=vmax,cmap='cmo.balance',zorder=3,transform=proj)
        
    #plt.suptitle("%s, Lag %02i" % (vname,lag))
    cb = viz.hcbar(pcm,ax=axs.flatten(),fraction=0.025)
    cb.set_label("%s [W/m2/degC], Lag %02i" % (vname,lag))
    #plt.show()
    
    savename = "%sLagRegression_AllMonths_regrid_%s_it%02i_lag%02i.png" % (figpath,vname,it,lag)
    plt.show()
    
    #plt.savefig(savename,dpi=150,bbox_inches='tight')
        
    #plt.show()
    plt.close()
# ======================================================================
#%% Make Individual Difference Plots

v           = 1

vname       = vnames[v]
for vname in vnames:
    vmax        = 10
    proj        = ccrs.PlateCarree()
    lag         = 0
    nlags = len(leadlags)
    
    # Set Reference
    for it in tqdm.tqdm(range(nlags)):
        lag         = leadlags[it]
        
        ref_ii      = 0
        ref_ex      = expids[ref_ii]
        refvar      = ds_byexp[ref_ii][vname].sel(lag=lag)
        #refvar      = ut.varcheck(,vname,expnames[ref_ex])
        #refvar      = ut.standardize_names(refvar)
        
        
        fig,axs     = init_globalmap(3,1,figsize=(12,12))
        
        for ii in range(3):
            
            ax        = axs[ii]
            ii_target = ii
            ex        = expids[ii_target]
            
            
            if ii == 0: # Just Plot the Control Values
                plotvar     = refvar
                lab         = expnames[ref_ex]
            else:
                targetvar   = ds_byexp[ii_target][vname].sel(lag=lag)
                #targetvar   = ut.varcheck(,vname,expnames[ex])
                #targetvar   = ut.standardize_names(targetvar)
                plotvar     = refvar - targetvar
                lab         = "%s - %s" % (expnames[ref_ex], expnames[ex])
                
                
            
            pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,vmin=-vmax,vmax=vmax,cmap='cmo.balance',zorder=3,transform=proj)
            ax.set_title(lab)
             
        cb = viz.hcbar(pcm,ax=axs.flatten(),fraction=0.025)
        cb.set_label("%s [W/m2/degC], Lag %02i" % (vname,lag))
        #plt.show()
        
        savename = "%sLagRegression_AllMonths_regrid_Differencesto31km_%s_it%02i_lag%02i.png" % (figpath,vname,it,lag)
        plt.savefig(savename,dpi=150,bbox_inches='tight')
        plt.close()



#%% Check if global mean is equal to the lag regression relationships

gmean_byexp = [proc.area_avg_cosweight(ds) for ds in ds_byexp]


lags = leadlags
#% Plot Lag Regression of global mean
plot_vnames = ['allsky','clearsky','cre']#['ttcre','ttrc','tscre','tsrc']

for vname in plot_vnames:
    
    fig,ax = plt.subplots(1,1,constrained_layout=True)
    
    for ii in range(3):
        
        ex      = expids[ii]
        plotvar = gmean_byexp[ii][vname]
        
        ax.plot(plotvar.lag,plotvar,label=expnames_long[ex],c=expcols[ex],lw=2.5)
    
    
    ax.set_xticks(lags[::2])
    ax.set_xlim([lags[0],lags[-1]])
    
    ax.axhline([0],lw=0.55,c='k')
    ax.axvline([0],lw=0.55,c='k')
    ax.grid(True,ls='dotted',c='gray',lw=0.55)
    ax.set_ylim([-.6,.6])
    ax.tick_params(labelsize=10)
    
    ax.set_xlabel("<-- Flux Leads | ENSO Leads -->")
    
    
    ax.set_title(vname)
    ax.legend()
    figname = "%sGlobal_Mean_regrid1x1_LeadLag_%s.png" % (figpath,vname)
    plt.savefig(figname,dpi=150,bbox_inches='tight')

    plt.show()



    
        
    
    
    
    
    
    
    
    



    
    
    



