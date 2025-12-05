#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Global Mean Lead Lag Analysis

Uses global means computed by 
    - ERA5 and CERES: [global_mean_ERA5.sh]
    

Data has been anomalized (deseasoned) and linearly detrended thru several scripts
    - ERA5: [anomalize_detrend_era5.sh]
    - AWI-CM3: [anom_detrend1_awiloop.sh]
    - CERES: [setup_input_data_CERES.sh]

TOA Fluxes computed in
    - [calc_toa_fluxes_awi_manual]

ENSO was computed using [calc_nino34.py]

Created on Thu Dec  4 14:29:20 2025

@author: gliu

"""


import sys
import time
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import xarray as xr
import sys
import glob 
import scipy as sp
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
from scipy.io import loadmat
import matplotlib as mpl
import importlib

import importlib
from tqdm import tqdm

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% Indicate Paths

gmeanpath   = '/home/niu4/gliu8/projects/ccfs/global_mean/'
ninopath    = "/home/niu4/gliu8/projects/scrap/nino34/"
figpath     = "/home/niu4/gliu8/figures/bydate/2025-12-AWI-Hackathon/"
proc.makedir(figpath)

flxnames    = ["allsky","clearsky","cre","ttcre","tscre",]
flxnames_long = ["All Sky", "Clear Sky", "Cloud-Radiative Effect", "Longwave CRE", "Shortwave CRE"]
ninoname    = "nino34"
standardize = False

# Simulation Names -----
expnames = ["TCo319_ctl1950d", "TCo319_ssp585",
            "TCo1279-DART-1950", "TCo1279-DART-2090", "TCo2559-DART-1950C","ERA5_1979_2024","CERES_EBAF"]
expnames_long = ["31km Control", "31km SSP585",
                 "9km 1950", "9km 2090", "5km 1950","ERA5 (1979-2024)","CERES-EBAF and ERA5 (2000-2024)"]
expcols = ["cornflowerblue", 'lightcoral',
           "slateblue", "firebrick",
           "midnightblue", "k","gray"]  # Includes Glorys and different shade based on resolution

#%% Load ENSO indices

ensoids = []
for expname in expnames:
    
    if expname == "CERES_EBAF":
        ensoid = ut.load_ensoid("ERA5_1979_2024", ninoname,standardize=standardize)
    else:
        ensoid = ut.load_ensoid(expname, ninoname,standardize=standardize)
    ensoids.append(ensoid)
    

#%% Load The Anomalized Fluxes...

dtday = 3600*24

def get_time_range(ds,timename="time"):
    return [ds[timename].data[0],ds[timename].data[-1]]

nvars      = len(flxnames)
flxa_byexp = []
for expname in expnames:
    
    flxa_bytype = []
    for ff in range(nvars):
        flxname = flxnames[ff]
        ncname  = "%s%s_%s_global_mean.nc" % (gmeanpath,expname,flxname,)
        try:
            ds      = xr.open_dataset(ncname)[flxname]
            ds      = ut.standardize_names(ds)
            
            if "CERES" in expname:
                ds = ut.shift_time_monthstart(ds)
            if "ERA5" in expname:
                ds = ds / dtday
            if "TCo" in expname:
                ds = ut.varcheck(ds,flxname,expname)
            
            flxa_bytype.append(ds)
        except:
            print("Warning! Could not find %s for %s" % (flxname,expname))
            print(ncname)
            break
        

            
            
            
        
        print("%s, %s:" % (expname,flxname))
        print("\t" +  str(get_time_range(ds)))
            
    flxa_byexp.append(flxa_bytype)

#%% Check Monthly Variance of each flux

monvar_byexp = []
for ex in range(len(expnames)):
    monvar_bytype = []
    for ff in range(nvars):
        
        flxin = flxa_byexp[ex][ff].groupby('time.month').std('time')
        monvar_bytype.append(flxin)
    monvar_byexp.append(monvar_bytype)

#%% Plot values
mons3   = proc.get_monstr()
fig,axs = viz.init_monplot(1,5,figsize=(20,3),constrained_layout=True)

for ii in range(5):
    ax = axs[ii]
    ax.set_title(flxnames[ii])
    ax.grid(True,ls='dotted',c='gray',lw=0.55)
    #ax.set_ylim([-.75,.75])
    ax.tick_params(labelsize=10)
    
    for ex in range(len(expnames)):
        plotvar = monvar_byexp[ex][ii].squeeze()
        ax.plot(mons3,plotvar,label=expnames_long[ex],c=expcols[ex],lw=2)
        
    
    ax.set_ylim([0,2])
    if ii == 3:
        ax.legend(ncol=1,fontsize=8)

figname = "%sGlobal_Mean_TOA_Fluxes_Monthly_Stdev.png" % (figpath)
plt.savefig(figname,dpi=150)
    
plt.show()
    

#%% Compute Lag Regression
lags          = np.arange(-24,25,)
nlags         = len(lags)
nexps         = len(expnames)
preproc       = False
betas_all     = np.zeros((nvars,nexps,nlags)) * np.nan
#betas_all_mon = np.zeros((nvars,nexps,12,nlags)) * np.nan # Skip for now because CERES is not monthly

# Do additional subsampling (currently only running for historical)
mciter        = 1000
minlen        = 8*12
mcbetas       = np.zeros((nvars,mciter,nlags)) * np.nan

for v in range(nvars):
    
    for ex in range(nexps):
        
        
        sst_in            = ensoids[ex].squeeze()
        flx_in            = flxa_byexp[ex][v].squeeze()
        
        sst_in            = ut.remove_duplicate_times(sst_in)
        flx_in            = ut.remove_duplicate_times(flx_in)
        
        sst_in,flx_in     = proc.match_time_month(sst_in,flx_in)
        
        flxname           = flxnames[v]
        
        
        dsout             = ut.calc_leadlag_regression_2d(sst_in,flx_in,lags,sep_mon=False)
        betas_all[v,ex,:] = dsout[flxname].data
        
        #dsoutmon          = ut.calc_leadlag_regression_2d(sst_in,flx_in,lags,sep_mon=True)
        #betas_all_mon[v,ex,:,:] = dsoutmon[flxname].data.T
        
        
        
        # Do some additional sampling
        if ex == 0:
            
            outdict_31km = ut.mcsampler(sst_in,minlen,mciter,target_timeseries=[flx_in,],preserve_month=True,scramble_year=False)
            ts1          = outdict_31km['samples']
            ts2          = outdict_31km['other_sampled_timeseries'][0]
            for mc in tqdm(range(mciter)):
                bb              = ut.calc_lag_regression_1d(ts2[mc,:],ts1[mc,:],lags)
                mcbetas[v,mc,:] = bb

#%% Plot the Lag Regression Output

# Part (1)
fig,axs = plt.subplots(1,5,figsize=(22,3),constrained_layout=True)

for ii in range(5):
    ax = axs[ii]
    
    for ex in range(nexps):
        
        plotvar = betas_all[ii,ex,:]
        ax.plot(lags,plotvar,label=expnames_long[ex],c=expcols[ex],lw=2)
        
        if ex == 0: # Plot uncertainty for the first experiment
            mu    = mcbetas[ii,:,:].mean(0)
            sigma = mcbetas[ii,:,:].std(0)
            ax.plot(lags,mu,color='dimgray',ls='dashed',lw=0.75)
            ax.fill_between(lags,mu-sigma,mu+sigma,color='dimgray',alpha=0.15,label="")
    
    if ii == 4:
        
        ax.legend(ncol=2,fontsize=8)
    
    ax.set_xticks(lags[::2])
    ax.set_xlim([lags[0],lags[-1]])
    
    ax.axhline([0],lw=0.55,c='k')
    ax.axvline([0],lw=0.55,c='k')
    
    ax.set_xlabel("<-- Flux Leads | ENSO Leads -->")
    
    ax.set_title(flxnames[ii])
    ax.grid(True,ls='dotted',c='gray',lw=0.55)
    ax.set_ylim([-.75,.75])
    #ax.set_aspect(1)
    ax.tick_params(labelsize=8,rotation=45)

#figname = "%sFluxes_Lag_Regression_%s.png" % (figpath,ensoid_name)
#plt.savefig(figname, dpi=150, bbox_inches='tight')
figname = "%sGlobal_Mean_TOA_Fluxes_Nino34_LagRegression.png" % (figpath)
plt.savefig(figname,dpi=150)
    
plt.show()
    

#%% Need to Split this plot into different sections


els     = ["solid","dashed","solid","dashed","solid",'dotted',"dotted"]
emks    = ["v","v","o","o","^","d","d"]

ff     = 0
for ff in range(nvars):
    
    fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(8,6.5))
    
    for ex in range(nexps):
        
        plotvar = betas_all[ff,ex,:]
        ax.plot(lags,plotvar,label=expnames_long[ex],c=expcols[ex],lw=2,marker=emks[ex],ls='solid')#ls=els[ex])
        
        if ex == 0: # Plot uncertainty for the first experiment
            mu    = mcbetas[ff,:,:].mean(0)
            sigma = mcbetas[ff,:,:].std(0)
            ax.plot(lags,mu,color='dimgray',ls='dashed',lw=0.75)
            ax.fill_between(lags,mu-sigma,mu+sigma,color='dimgray',alpha=0.15,label="")
            
            
    ax.set_xticks(lags[::2])
    ax.set_xlim([lags[0],lags[-1]])
    
    ax.axhline([0],lw=0.55,c='k')
    ax.axvline([0],lw=0.55,c='k')
    ax.legend(fontsize=12)
    ax.set_title(flxnames_long[ff])
    ax.set_xlabel("<-- Flux Leads | ENSO Leads -->")
    ax.set_ylim([-.75,.75])
    
    figname = "%sGlobal_Mean_TOA_%s_Nino34_LagRegression.png" % (figpath,flxnames[ff])
    plt.savefig(figname,dpi=150,bbox_inches='tight')

#plt.show()

#nexps = len
#for ex in range(nexps):
    


    





