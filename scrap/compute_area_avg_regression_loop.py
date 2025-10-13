#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 12:01:46 2025

@author: gliu
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copied from 
Area Average Regressions


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

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% Commons Variables


#%% Shared variable names and experiment
# Copied from visualize composites

#datpath         = "/Users/gliu/Downloads/02_Research/01_Projects/07_ENSO/01_Data/TP_Crop/composites/"


# Simulation Names -----
expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"]
expcols         = ["cornflowerblue","red","cornflowerblue","red","cornflowerblue"]
expls           = ["solid","solid",'dashed','dashed','dotted','dotted']
emks            = ["v","v","o","o","^","^"]

# Initial Variable Analysis -----
vnames          = ["sst","ssr","str","tx_sur","D20","Dmaxgrad"]
vunits          = [r"$\degree C$",r"$\frac{W}{m^2}$",r"$\frac{W}{m^2}$",r"$\frac{m}{s^2}$","m","m"]
vnames_long     = ["SST","Surface Shortwave","Surface Longwave","Zonal Wind Stress","Thermocline (20$\degree$ Isotherm)","Thermocline (Max Vertical Gradient)"]
vmaxes          = [2,40,20,0.02,20,20]

# ENSO Names -----
ninoname        = [r"$El$ $Ni\tilde{n}o$",r"$La$ $Ni\tilde{n}a$"]
ninocol         = ["cornflowerblue","firebrick"]
ninoshort       = ['nino','nina']

# Conversion for STR and SSR considering 3h Accumulation -----
conversion      = 1/(3*3600) # 3 h accumulation time...? #1/(24*30*3600)
# https://forum.ecmwf.int/t/surface-radiation-parameters-joule-m-2-to-watt-m-2/1588

# Bounding Boxes from Jin et al. 2020 Eqn. 6.6  -----
bbox_cep        = [150      , -130+360 , -5, 5]   # Central Equatorial Pacific, for [tau_x], 
bbox_nino3      = [-150+360 , -90+360  , -5, 5]  # Nino 3 Box: For SST, <tau_x>
bbox_nino34     = [-170+360 , -120+360 , -5, 5]
bbox_epac       = [-155+360 , -80+360  , -5, 5]  # Eastern Pacific (for h_e calculation)
bbox_wpac       = [120      , -155+360 , -5, 5]  # Western Pacific (for h_w calculation)

bboxes          = [bbox_cep,bbox_nino3,bbox_nino34,bbox_epac,bbox_wpac]
bbnames_long    = ["Central Equatorial Pacific",r"$Ni\tilde{n}o3$",r"$Ni\tilde{n}o3.4$","Tropical Eastern Pacific","Tropical Western Pacific"]
bbnames         = ["CEP","nino3","nino34","EPac","WPac"]

#%% Indicate paths

figpath         = "/home/niu4/gliu8/figures/bydate/2025-10-14/"
proc.makedir(figpath)

datpath         = "/home/niu4/gliu8/projects/scrap/TP_crop/"
datpath_anom    = datpath + "anom_detrend2/"



ninopath        = "/home/niu4/gliu8/projects/scrap/nino34/"

nexps           = len(expnames)


#%% Load ENSO Indices

ninoid_name   = "nino34"
unstandardize = True
ds_enso       = []
ensoids       = []
for ex in range(nexps):
    
    # ninonc      = "%s%s_%s.nc" % (ninopath,expnames[ex],ninoid_name)
    # ds          = xr.open_dataset(ninonc).load()
    
    ds = ut.load_ensoid(expnames[ex],ninoid_name,standardize=True)
    
    ds_enso.append(ds)
    ensoids.append(ds * ds['std'].data.item())
    
if unstandardize:
    ensoids = [ds * ds['std'].item() for ds in ensoids]

#%% Analyze relationship between loaded Nino Index above and a selected variable
# with area-averages across several regions...


def remove_duplicate_times(ds,verbose=True):
    # From : https://stackoverflow.com/questions/51058379/drop-duplicate-times-in-xarray
    _, index = np.unique(ds['time'], return_index=True)
    print("Found %i duplicate times. Taking first entry." % (len(ds.time) - len(index)))
    return ds.isel(time=index)
    
# Select Variable and Regions

vnames = ["sst","str","ssr","tx_sur","ttr","tsr","ttrc","tsrc","eis"]
for vname in vnames:
    
    if vname == "eis":
        vunit   = "$EIS$"
    elif vname == "tx_sur":
        vunit   = "$m/s2"
    elif vname == "sst":
        vunit   = "degC"
    else:
        vunit   = "J/m2"
    
    
    
    regions = bboxes
    nreg    = len(regions)
    
    
    # Select and take area-average for each region
    ds_byreg    = []
    aavgs_byreg = []
    for r in range(nreg):
        
        aavgs  = []
        ds_var = []
        for ex in tqdm.tqdm(range(nexps)):
            
            # Load Dataset
            ncname = "%s%s_%s_anom.nc" % (datpath_anom,expnames[ex],vname)
            
            try:
                ds = xr.open_dataset(ncname)
            except:
                ds_var.append(None)
                aavgs.append(None)
                print("Could not find %s %s" % (vname,expnames[ex]))
                continue
            
            # Do renaming
            if vname.upper() in list(ds.keys()):
                print("Renaming %s to %s" % (vname.upper(),vname))
                ds = ds.rename({vname.upper():vname})
            
            ds = ut.standardize_names(ds)[vname]
            ds = ut.varcheck(ds,vname,expnames[ex])
            
            # Subset Regions
            bbox = regions[r]
            ds  = proc.sel_region_xr(ds,bbox).load()
            ds_var.append(ds)
            
            # Take Area Average
            dsaavg = proc.area_avg_cosweight(ds)
            aavgs.append(dsaavg)
    
        ds_byreg.append(ds_var)
        aavgs_byreg.append(aavgs)
    
    "#%% Compute Regressions for each case (nino index to variable over each region)
    
    lags       = np.arange(-25,26,1)
    
    beta_byreg = []
    for r in range(nreg):
        beta_byexp = []
        for ex in range(nexps):
            
            sst_in = ensoids[ex]
            tau_in = aavgs_byreg[r][ex]
            
            if tau_in is None:
                beta_byexp.append(None)
                print("No %s for %s" % (vname,expnames[ex]))
                continue
            
            sst_in,tau_in = proc.match_time_month(sst_in,tau_in)
            
            if len(sst_in) != len(tau_in):
                print("Warning, length is still not matching. Checking for duplicate times.")
                sst_in = remove_duplicate_times(sst_in)
                tau_in = remove_duplicate_times(tau_in)
                
            
            betas  = ut.calc_lag_regression_1d(tau_in,sst_in,lags)
            
            beta_byexp.append(betas)
        beta_byreg.append(beta_byexp)
        
    #%% Resolution Effects
    
    # SST Nino3 vs Heat Flux Nino3 (copied from Taux example above)
    
    # Region and Scenario Options
    comparison      = "control"
    regid           = 1
    for regid in range(nreg):
        regname         = bbnames_long[regid]
        regname_short   = bbnames[regid]
        
        # Indicate Variable
        invar_aavgs     = aavgs_byreg
        
        # Other Options
        outcols  = ['red','blue','yellow']
        
        # Test significance
        expids_in       = [0,2,4]
        tlens           = [len(ensoids[e]) for e in expids_in]
        minlen          = np.min(tlens) # Get Minimum Length
        # Test for 31 km
        mciter            = 1000
        ts_full           = ensoids[0]
        tau_in            = invar_aavgs[regid][0]
        if tau_in is not None:
            ts_full,tau_in    = proc.match_time_month(ts_full,tau_in)
            target_timeseries = [tau_in.data,]
            ts_full           = ts_full.data
            outdict_31km      = ut.mcsampler(ts_full,minlen,mciter,target_timeseries=target_timeseries,preserve_month=True,scramble_year=False)
            indict            = outdict_31km
            ts1               = indict['samples'] # iter x time
            ts2               = indict['other_sampled_timeseries'][0] # iter x time
            outval = []
            for mc in tqdm.tqdm(range(mciter)):
                betas  = ut.calc_lag_regression_1d(ts2[mc,:],ts1[mc,:],lags)
                outval.append(betas)
            out_result = np.array(outval)
        else:
            out_result = None
        
        # Make the Plot -------------------------------------------------------------
        fig,ax        = plt.subplots(1,1,figsize=(8,4.5),constrained_layout=True)
        # Plot the wind stress-sst relationship
        for ex in expids_in:
            plotvar = beta_byreg[regid][ex]
            if plotvar is None:
                continue
            ax.plot(lags,plotvar,label=expnames_long[ex],c=expcols[ex],ls=expls[ex],marker=emks[ex])
        ax.legend()
        
        # Plot Significance
        plotvar = out_result
        if plotvar is not None: # PLot only if the data exists..
            mu      = np.nanmean(plotvar,0)
            sigma   = np.nanstd(plotvar,0)
            
            #clow,chi = np.percentile(plotvar,[0.025,0.975],axis=0)
            #ax.fill_between(lags,clow,chi,color=outcols[ii],alpha=0.15)
            ax.plot(lags,mu,color='dimgray',ls='dashed',lw=0.75)
            ax.fill_between(lags,mu-sigma,mu+sigma,color='dimgray',alpha=0.15,label="")
            
        # Set Limits
        ax.set_xlim([-24,24])
        ax.set_xticks(np.arange(-24,25,2))
        #ax.set_ylim([-.015,.015])
        
        ax.axhline([0],ls='solid',lw=0.75,c="k")
        ax.axvline([0],ls='solid',lw=0.75,c="k")
        
        ax.grid(True,ls='dotted',c='gray',lw=0.55)
        ax.set_title(r"%s (%s) vs. sst (%s) Lag Regression" % (vname,regname_short,ninoid_name))
        
        ax.set_xlabel("<-- %s Leads | SST Leads -->" % vname)
        ax.set_ylabel(r"Regression Coefficient [%s per $\degree C]$" % vunit)
        
        figname = "%s%s_%s_%s_lag_regression_relationship.png" % (figpath,vname,ninoid_name,regname_short)
        plt.savefig(figname,dpi=150,bbox_inches='tight')
        #plt.show()
        plt.close()


#%% 10.10.2025 General Formulation: correlate (nino index) with [vname]




