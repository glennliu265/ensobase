#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Area Average Regressions

take lag regressions between ENSO indices and certain variables, based on the Jin et al. 2020 Paper

Created on Thu Oct  2 11:28:16 2025

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
#expcols         = ["cornflowerblue","red","cornflowerblue","red","cornflowerblue"]
expls           = ["solid","solid",'dashed','dashed','dotted','dotted']
emks            = ["v","v","o","o","^","^"]

expcols = ["cornflowerblue",'lightcoral',
           "slateblue","firebrick",
           "midnightblue","k"] # Includes Glorys and different shade based on resolution
 


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

#%% More Functions
# def preprocess_enso(ds):
#     # Remove Mean Seasonal Cycle and the Quadratic Trend
#     dsds   = proc.xrdeseason(ds)
    
#     # Compute Spectra
#     def detrend_quadratic(ds):
#         x = np.arange(len(ds))
#         y = ds.data
#         if np.any(np.isnan(y)):
#             return np.ones(y.shape)*np.nan
#         ydetrended,model=proc.detrend_poly(x,y,2)
#         return ydetrended
        
#     st = time.time()
#     dsanom = xr.apply_ufunc(
#         detrend_quadratic,  # Pass the function
#         dsds,  # The inputs in order that is expected
#         # Which dimensions to operate over for each argument...
#         input_core_dims=[['time'],],
#         output_core_dims=[['time'],],  # Output Dimension
#         vectorize=True,  # True to loop over non-core dims
#         )
#     print("Detrended in %.2fs" % (time.time()-st))
#     #dsanom = proc.xrdetrend_1d(ds,2) 
#     return dsanom

# def swap_rename(ds,chkvar,newvar):
#     if chkvar in list(ds.coords):
#         print("Renaming [%s] to [%s]" % (chkvar,newvar))
#         ds = ds.rename({chkvar:newvar})
#     return ds

# def standardize_names(ds):
    
#     ds = swap_rename(ds,'time_counter','time')
#     ds = swap_rename(ds,"TIME_COUNTER",'time')
#     ds = swap_rename(ds,"LON","lon")
#     ds = swap_rename(ds,"LAT","lat")
#     ds = swap_rename (ds,"LAT232_409","lat")
#     return ds

# def mcsampler(ts_full,sample_len,mciter,preserve_month=True,scramble_year=False,target_timeseries=None):
#     # Given a monthly timeseries [time] and sample length (int), take [mciter] # of random samples.
#     # if preserve_month = True, preserve the 12-month sequence as a chunk
#     # if scramble_year = True, randomize years that you are selecting from (do not preserve year order)
#     # if target_timeseries is not None: also select random samples from list of timeseries (must be same length as ts_full)
    
#     # Function Start
#     ntime_full        = len(ts_full)
    
#     # 1 -- month agnostic (subsample sample length, who cares when)
#     if not preserve_month:
        
#         print("Month with not be preserved.")
#         istarts    = np.arange(ntime_full-sample_len)
        
#         sample_ids = []
#         samples    = []
#         for mc in range(mciter):
#             # ts_full[istarts[-1]:(istarts[-1]+sample_len)] Test last possible 
#             iistart = np.random.choice(istarts)
#             idsel   = np.arange(iistart,iistart+sample_len) 
#             msample = ts_full[idsel]
            
            
            
#             sample_ids.append(idsel)
#             samples.append(msample) # [iter][sample]
            
#         samples = np.array(samples) # [iter x sample]
#         # Returns 
            
#     elif preserve_month:
#         # 2 -- month aware (must select starting points of January + maintain the chunk, preserving the month + year to year autocorrelation)
#         if not scramble_year:
            
#             # Only start on the year  (to preserve month sequence)
#             istarts    = np.arange(0,ntime_full-sample_len,12)
            
#             # -------------------- Same as Above
#             sample_ids = []
#             samples    = []
#             for mc in range(mciter):
#                 # ts_full[istarts[-1]:(istarts[-1]+sample_len)] Test last possible 
#                 iistart = np.random.choice(istarts)
#                 idsel   = np.arange(iistart,iistart+sample_len) 
#                 msample = ts_full[idsel]
                
#                 sample_ids.append(idsel)
#                 samples.append(msample) # [var][iter][sample]
#             samples = np.array(samples) # [var x iter x sample]
#             # -------------------- 
            
#         # 3 -- month aware, year scramble (randomly select the year of each month, but preserve each month)
#         elif scramble_year: # Scrample Year and Month
            
#             # Reshape to the year and month
#             nyr_full        = int(ntime_full/12)
#             ts_yrmon        = ts_full.reshape(nyr_full,12)
#             ids_ori         = np.arange(ntime_full)
#             ids_ori_yrmon   = ids_ori.reshape(ts_yrmon.shape)
            
#             nyr_sample      = int(sample_len/12)
#             sample_ids      = []
#             samples         = []
#             for mc in range(mciter): # For each loop
                
#                 # Get start years
#                 startyears = np.random.choice(np.arange(nyr_full),nyr_sample)
#                 # Select random years equal to the sample length and combine
#                 idsel      = ids_ori_yrmon[startyears,:].flatten() 
#                 # ------
#                 msample    = ts_full[idsel]
#                 sample_ids.append(idsel)
#                 samples.append(msample) # [var][iter][sample]
#             samples = np.array(samples) # [var x iter x sample]
#             # -----
    
#     outdict = dict(sample_ids = sample_ids, samples=samples)    
#     if target_timeseries is not None:
        
#         sampled_timeseries = []
#         for ts in target_timeseries:
#             if len(ts) != len(ts_full):
#                 print("Warning... timeseries do not have the same length")
#             randsamp = [ts[sample_ids[mc]] for mc in range(mciter)]
#             randsamp = np.array(randsamp)
#             sampled_timeseries.append(randsamp) # [var][iter x time]
#     outdict['other_sampled_timeseries'] = sampled_timeseries
    
#     return outdict


            
# def calc_lag_regression_1d(var1,var2,lags): # CAn make 2d by mirroring calc_lag_covar_annn
#     # Calculate the regression where
#     # (+) lags indicate var1 lags  var2 (var 2 leads)
#     # (-) lags indicate var1 leads var2 (var 1 leads)
    
#     ntime = len(var1)
#     betalag = []
#     poslags = lags[lags >= 0]
#     for l,lag in enumerate(poslags):
#         varlag   = var1[lag:]
#         varbase  = var2[:(ntime-lag)]
        
#         # Calculate correlation
#         beta = np.polyfit(varbase,varlag,1)[0]   
#         betalag.append(beta.item())
    
    
#     neglags = lags[lags < 0]
#     neglags_sort = np.sort(np.abs(neglags)) # Sort from least to greatest #.sort
#     betalead = []
    
#     for l,lag in enumerate(neglags_sort):
#         varlag   = var2[lag:] # Now Varlag is the base...
#         varbase  = var1[:(ntime-lag)]
#         # Calculate correlation
#         beta = np.polyfit(varlag,varbase,1)[0]   
#         betalead.append(beta.item())
        
        
#     # Append Together
#     return np.concat([np.flip(np.array(betalead)),np.array(betalag)])


#%% Load ENSO Indices

ninoid_name   = "nino3"
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
    
# ================================
# Part (1): Zonal Wind Stress
# ================================
#%% Load Variable of Choice for Analysis. Let's start with tx_sur

vname = "tx_sur"
regnames     = ["Central Equatorial Pacific",r"$Ni\tilde{n}o 3$"]
ds_var_cep   = []
ds_var_nino3 = []

for ex in tqdm.tqdm(range(nexps)):
    
    # Load Dataset
    ncname = "%s%s_%s_anom.nc" % (datpath_anom,expnames[ex],vname)
    ds = xr.open_dataset(ncname)
    
    # Do renaming
    if vname.upper() in list(ds.keys()):
        print("Renaming %s to %s" % (vname.upper(),vname))
        ds = ds.rename({vname.upper():vname})
    
    # Subset Regions
    ds1 = proc.sel_region_xr(ds,bbox_cep).load()
    print(ds1[vname].shape)
    ds_var_cep.append(ds1)
    
    ds2 = proc.sel_region_xr(ds,bbox_nino3).load()
    ds_var_nino3.append(ds2)
    
# Preprocess
tau_anoms = []
for ds_var in [ds_var_cep,ds_var_nino3]:
    
    ds_var      = [ut.standardize_names(ds)[vname] for ds in ds_var]
    
    # if vname in ["D20","Dmaxgrad"]:
    #     ds_anoms = [ut.preprocess_enso(ds['nz1']) for ds in ds_var]
    # else:
    #     ds_anoms = [ut.preprocess_enso(ds[vname]) for ds in ds_var]
    
    tau_anoms.append(ds_var)

# Take Area Averages
tau_aavgs = []
for dsreg in tau_anoms:
    aavgs = [proc.area_avg_cosweight(ds) for ds in dsreg]
    tau_aavgs.append(aavgs)

#%% Exploratory Plot: Monthly Variance


fig,axs = viz.init_monplot(3,1,figsize=(6,12))
mons3 = proc.get_monstr()


# Plot the Monthly Variance of Wind Stress
for rr in range(2):
    ax = axs[rr]
    ax.set_title("Wind Stress Monthly Stdev. (%s)" % regnames[rr])
    for ex in range(nexps):
        plotvar = tau_aavgs[rr][ex].groupby('time.month').std('time')
        
        ax.plot(mons3,plotvar,
                label=expnames_long[ex],c=expcols[ex],
                ls=expls[ex],marker=emks[ex])

    
    if rr == 1:
        ax.legend(ncol=3)
    ax.set_ylim([0,0.02])
    ax.set_ylabel(r"Zonal Wind Stress ($\frac{m}{s^2}$)")
    
    

# Plot SST Monthly Variance
ax = axs[2]
ax.set_title("SST Monthly Stdev. (%s)" % ninoid_name) 
for ex in range(nexps):
    plotvar = ensoids[ex].groupby('time.month').std('time')
    
    ax.plot(mons3,plotvar,
            label=expnames_long[ex],c=expcols[ex],
            ls=expls[ex],marker=emks[ex])
    

ax.set_ylim([0,2.5])
ax.set_ylabel(r"SST ($\degree C$)")

figname = "%sMonvar_TauRegs_SST_%s.png" % (figpath,ninoid_name)
plt.savefig(figname,dpi=150,bbox_inches='tight')

plt.show()


#%% Calculate Lag Regression

lags   = np.arange(-24,25,1)
#out = calc_lag_regression_1d(var1,var2,lags)

#%% Calculate Relationship for each one

regid = 0
beta_byreg = []
for regid in range(2):
    beta_byexp = []
    for ex in range(nexps):
        
        sst_in = ensoids[ex]
        tau_in = tau_aavgs[regid][ex]
        
        sst_in,tau_in = proc.match_time_month(sst_in,tau_in)
        
        betas         = ut.calc_lag_regression_1d(tau_in,sst_in,lags)
        
        beta_byexp.append(betas)
    beta_byreg.append(beta_byexp)
    
#%% Also compute lag 0 regression for each month separately

regid = 0
betamon_byreg = []
for regid in range(2):
    betamon_byexp = []
    for ex in range(nexps):
        
        sst_in        = ensoids[ex]
        tau_in        = tau_aavgs[regid][ex]
        
        sst_in,tau_in = proc.match_time_month(sst_in,tau_in)
        
        ntime = sst_in.shape[0]
        nyr   = int(ntime/12)
        
        sst_monyr = sst_in.data.reshape((nyr,12))
        tau_monyr = tau_in.data.reshape((nyr,12))
        
        lagreg_sepmon = np.ones(12) * np.nan
        for im in range(12):
            
            betas         = ut.calc_lag_regression_1d(tau_monyr[:,im],sst_monyr[:,im],np.array([0,]))
            lagreg_sepmon[im] = betas[0]
        
        betamon_byexp.append(lagreg_sepmon)
    betamon_byreg.append(betamon_byexp)
    
#%% Do a quick plot of this



fig,axs = viz.init_monplot(1,2)
for regid in range(2):
    
    ax = axs[regid]
    
    for ex in range(nexps):
        
        ax.plot(mons3,betamon_byreg[regid][ex],label=expnames_long[ex],c=expcols[ex],ls=expls[ex],marker=emks[ex])
        
    ax.legend()
    ax.set_title(regnames[regid])

plt.show()


#%% Testing Section Below =====================================================
#%% Subsample random chunk and redo calculations (to see if this is due to resolution Change)
#=====================================================
ntimes_all = [len(ts) for ts in ensoids]

# Need to test 3 versions

comparison = "ssps85"
regid      = 1


# Function Arguments
if comparison == "Control":
    sample_len        = 240 # months
    ts_full           = ensoids[0].data
    target_timeseries = [tau_aavgs[regid][0].data,] # Additional Timeseries to Sample, must be same time as ts_full#[]
elif comparison == "ssps85":
    sample_len        = 120 # months
    ts_full           = ensoids[1]
    ts_full,tauinn    = proc.match_time_month(ts_full,tau_aavgs[regid][1])
    ts_full           = ts_full.data
    target_timeseries = [tauinn.data,] # Additional Timeseries to Sample, must be same time as ts_full#[]
preserve_month    = False
scramble_year     = False
mciter            = 1000


# Test the three methods
mciter   = 1000
outdict1 = ut.mcsampler(ts_full,sample_len,mciter,target_timeseries=target_timeseries,preserve_month=False)
outdict2 = ut.mcsampler(ts_full,sample_len,mciter,target_timeseries=target_timeseries,preserve_month=True,scramble_year=False)
outdict3 = ut.mcsampler(ts_full,sample_len,mciter,target_timeseries=target_timeseries,preserve_month=True,scramble_year=True)
outdicts = [outdict1,outdict2,outdict3]

out_result = []
for ii in range(3):
    print(ii)
    indict = outdicts[ii]
    
    ts1    = indict['samples'] # iter x time
    ts2    = indict['other_sampled_timeseries'][0] # iter x time
    
    outval = []
    for mc in tqdm.tqdm(range(mciter)):
        betas  = ut.calc_lag_regression_1d(ts2[mc,:],ts1[mc,:],lags)
        outval.append(betas)
        
    outval = np.array(outval)
    out_result.append(outval)
        
#%% Plot Sensitivity to Random Sampling Method

fig,ax   = plt.subplots(1,1,figsize=(10,4.5),constrained_layout=True)

outcols  = ['red','blue','yellow']
outnames = ["Completely Random","Preserve Month","Preserve Month + Scramble Years"]

for ii in range(3):
    plotvar = out_result[ii]
    
    mu    = np.nanmean(plotvar,0)
    sigma = np.nanstd(plotvar,0)
    
    #clow,chi = np.percentile(plotvar,[0.025,0.975],axis=0)
    #ax.fill_between(lags,clow,chi,color=outcols[ii],alpha=0.15)
    ax.plot(lags,mu,color=outcols[ii],ls='dashed',lw=0.75)
    ax.fill_between(lags,mu-sigma,mu+sigma,color=outcols[ii],alpha=0.15,label=outnames[ii])
    
if comparison == "Control":

    ex    = 0
elif comparison == "ssps85":

    ex    = 1
    
plotvar = beta_byreg[regid][ex]
ax.plot(lags,plotvar,label=expnames_long[ex],c=expcols[ex],ls=expls[ex],marker=emks[ex])

if comparison == "Control":
    ex    = 2
elif comparison == "ssps85":
    ex    = 3

plotvar = beta_byreg[regid][ex]
ax.plot(lags,plotvar,label=expnames_long[ex],c=expcols[ex],ls=expls[ex],marker=emks[ex])


ax.legend()


ax.set_xlabel("<-- Taux Leads | SST Leads -->")
ax.set_ylabel(r"Regression Coefficient [$\frac{m}{s^2}$ per $\degree C]$")


ax.set_xlim([-24,24])
ax.set_xticks(np.arange(-24,25,2))
ax.set_ylim([-.015,.015])

ax.axhline([0],ls='solid',lw=0.75,c="k")
ax.axvline([0],ls='solid',lw=0.75,c="k")

ax.grid(True,ls='dotted',c='gray',lw=0.55)


figname = "%sResolution_Significance_Tau_v_Nino3_%s_v_1950_reg%s.png" % (figpath,comparison,regnames_short[regid])
plt.savefig(figname)

plt.show()

#%% Update 2025.10.25: Do comparison on resolution differences for the shortest period for [Control]

comparison      = "control"
regid           = 1
regnames        = ["Central Equatorial Pacific",r"$Ni\tilde{n}o 3$"]
regnames_short  = ["CEP","Nino3"]

# Test significance
expids_in       = [0,2,4]
tlens           = [len(ensoids[e]) for e in expids_in]
minlen          = np.min(tlens) # Get Minimum Length
# Test for 31 km
mciter            = 1000
ts_full           = ensoids[0]
tau_in            = tau_aavgs[regid][0]
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

fig,ax        = plt.subplots(1,1,figsize=(8,4.5),constrained_layout=True)
# Plot the wind stress-sst relationship
for ex in expids_in:
    plotvar = beta_byreg[regid][ex]
    ax.plot(lags,plotvar,label=expnames_long[ex],c=expcols[ex],ls=expls[ex],marker=emks[ex])
ax.legend()

# Plot Significance
plotvar = out_result
mu      = np.nanmean(plotvar,0)
sigma   = np.nanstd(plotvar,0)

#clow,chi = np.percentile(plotvar,[0.025,0.975],axis=0)
#ax.fill_between(lags,clow,chi,color=outcols[ii],alpha=0.15)
ax.plot(lags,mu,color=outcols[ii],ls='dashed',lw=0.75)
ax.fill_between(lags,mu-sigma,mu+sigma,color=outcols[ii],alpha=0.15,label=outnames[ii])

# Set Limits
ax.set_xlim([-24,24])
ax.set_xticks(np.arange(-24,25,2))
ax.set_ylim([-.015,.015])

ax.axhline([0],ls='solid',lw=0.75,c="k")
ax.axvline([0],ls='solid',lw=0.75,c="k")

ax.grid(True,ls='dotted',c='gray',lw=0.55)
ax.set_title(regnames[regid])

ax.set_xlabel("<-- Taux Leads | SST Leads -->")
ax.set_ylabel(r"Regression Coefficient [$\frac{m}{s^2}$ per $\degree C]$")


figname = "%s%s_txsur_%s_lag_regression_relationship.png" % (figpath,ninoid_name,regnames_short[regid])
plt.savefig(figname,dpi=150,bbox_inches='tight')
plt.show()





#outdict3   = ut.mcsampler(ts_full,sample_len,mciter,target_timeseries=target_timeseries,preserve_month=True,scramble_year=True)

#%% Plot Output

regnames        = ["Central Equatorial Pacific",r"$Ni\tilde{n}o 3$"]
regnames_short = ["CEP","Nino3"]
fig,axs        = plt.subplots(1,2,figsize=(16,4.5),constrained_layout=True)

for regid in range(2):
    ax = axs[regid]
    for ex in range(nexps):
        
        
        plotvar = beta_byreg[regid][ex]
        ax.plot(lags,plotvar,label=expnames_long[ex],c=expcols[ex],ls=expls[ex],marker=emks[ex])

    ax.legend()

    ax.set_xlim([-24,24])
    ax.set_xticks(np.arange(-24,25,2))
    ax.set_ylim([-.015,.015])
    
    ax.axhline([0],ls='solid',lw=0.75,c="k")
    ax.axvline([0],ls='solid',lw=0.75,c="k")

    ax.grid(True,ls='dotted',c='gray',lw=0.55)
    ax.set_title(regnames[regid])
    
    ax.set_xlabel("<-- Taux Leads | SST Leads -->")

axs[0].set_ylabel(r"Regression Coefficient [$\frac{m}{s^2}$ per $\degree C]$")
plt.show()

# ================================
#%% Part (2), Thermocline Gradient
# ================================
    
vname  = "D20"

dsvarw = []
dsvare = []
for ex in tqdm.tqdm(range(nexps)):
    
    
    if vname == "tx_sur":
        vname_file = "tx_surf_1m"
        ncname = "%s%s_%s.nc" % (datpath,expnames[ex],vname_file)
    else:
        ncname = "%s%s_%s.nc" % (datpath,expnames[ex],vname)
    
    ds = xr.open_dataset(ncname)
    
    if vname.upper() in list(ds.keys()):
        print("Renaming %s to %s" % (vname.upper(),vname))
        ds = ds.rename({vname.upper():vname})
    
    ds1 = proc.sel_region_xr(ds,bbox_wpac).load()
    
    dsvarw.append(ds1)
    
    ds2 = proc.sel_region_xr(ds,bbox_epac).load()
    
    dsvare.append(ds2)
    
# Preprocess
thermocline_anoms = []
for ds_var in [dsvarw,dsvare]:
    
    ds_var      = [standardize_names(ds) for ds in ds_var]
    
    if vname in ["D20","Dmaxgrad"]:
        ds_anoms = [preprocess_enso(ds['nz1']) for ds in ds_var]
    else:
        ds_anoms = [preprocess_enso(ds[vname]) for ds in ds_var]
    
    thermocline_anoms.append(ds_anoms)

# Take Area Averages
thermocline_aavgs = []

for dsreg in thermocline_anoms:
    aavgs = [proc.area_avg_cosweight(ds) for ds in dsreg]
    thermocline_aavgs.append(aavgs)

# Do he - hw
thermocline_diff  = []
for ex in range(nexps):
    print(ex)
    tdiff = thermocline_aavgs[1][ex] - thermocline_aavgs[0][ex]
    thermocline_diff.append(tdiff)
    
    
    

#%% Look at the difference


regid = 0
beta_byreg_thermocline = []
for ex in range(nexps):
    
    sst_in = ensoids[ex]
    h_in   = thermocline_diff[ex]
    sst_in,h_in = proc.match_time_month(sst_in,h_in)
    betas       = calc_lag_regression_1d(h_in,sst_in,lags)
    beta_byreg_thermocline.append(betas)

#%% Plot Result

fig,ax  = plt.subplots(1,1,figsize=(16,4.5),constrained_layout=True)

for ex in range(nexps):
    
    
    plotvar = beta_byreg_thermocline[ex]
    ax.plot(lags,plotvar,label=expnames_long[ex],c=expcols[ex],ls=expls[ex],marker=emks[ex])

ax.legend()

ax.set_xlim([-24,24])
ax.set_xticks(np.arange(-24,25,2))
#ax.set_ylim([-.015,.015])

ax.axhline([0],ls='solid',lw=0.75,c="k")
ax.axvline([0],ls='solid',lw=0.75,c="k")

ax.grid(True,ls='dotted',c='gray',lw=0.55)
#ax.set_title(regnames[regid])

ax.set_xlabel("<-- [$h_e - h_w$] Leads | SST Leads -->")

ax.set_ylabel(r"Regression Coefficient [$m$ per $\degree C]$")
plt.show()



#%% Redo same significance test as above


comparison = "ssps85"#"ssps85"

# Function Arguments
if comparison == "Control":
    sample_len        = 240 # months
    ts_full           = ensoids[0].data
    
    target_timeseries = [thermocline_diff[0].data,] # Additional Timeseries to Sample, must be same time as ts_full#[]
elif comparison == "ssps85":
    sample_len        = 120 # months
    #ts_full           = ensoids[1].data
    ts_full,thermo_in = proc.match_time_month( ensoids[1],thermocline_diff[1])
    ts_full = ts_full.data
    
    target_timeseries = [thermo_in.data,] # Additional Timeseries to Sample, must be same time as ts_full#[]
preserve_month    = False
scramble_year     = False
mciter            = 1000
            
# Test the three methods
mciter   = 1000
outdict1 = mcsampler(ts_full,sample_len,mciter,target_timeseries=target_timeseries,preserve_month=False)
outdict2 = mcsampler(ts_full,sample_len,mciter,target_timeseries=target_timeseries,preserve_month=True,scramble_year=False)
outdict3 = mcsampler(ts_full,sample_len,mciter,target_timeseries=target_timeseries,preserve_month=True,scramble_year=True)
outdicts = [outdict1,outdict2,outdict3]

out_result = []
for ii in range(3):
    
    indict = outdicts[ii]
    
    ts1    = indict['samples'] # iter x time
    ts2    = indict['other_sampled_timeseries'][0] # iter x time
    
    outval = []
    for mc in tqdm.tqdm(range(mciter)):
        betas  = calc_lag_regression_1d(ts2[mc,:],ts1[mc,:],lags)
        outval.append(betas)
        
    outval = np.array(outval)
    out_result.append(outval)
        
#%% Plot Sensitivity to Random Sampling Method


fig,ax  = plt.subplots(1,1,figsize=(10,4.5),constrained_layout=True)


outcols  = ['red','blue','yellow']
outnames = ["Completely Random","Preserve Month","Preserve Month + Scramble Years"]

for ii in range(3):
    plotvar = out_result[ii]
    
    mu    = np.nanmean(plotvar,0)
    sigma = np.nanstd(plotvar,0)
    
    #clow,chi = np.percentile(plotvar,[0.025,0.975],axis=0)
    #ax.fill_between(lags,clow,chi,color=outcols[ii],alpha=0.15)
    ax.plot(lags,mu,color=outcols[ii],ls='dashed',lw=0.75)
    ax.fill_between(lags,mu-sigma,mu+sigma,color=outcols[ii],alpha=0.15,label=outnames[ii])
    
if comparison == "Control":

    ex    = 0
elif comparison == "ssps85":

    ex    = 1
    
plotvar = beta_byreg_thermocline[ex]
ax.plot(lags,plotvar,label=expnames_long[ex],c=expcols[ex],ls=expls[ex],marker=emks[ex])

if comparison == "Control":

    ex    = 2
elif comparison == "ssps85":

    ex    = 3
    

plotvar = beta_byreg_thermocline[ex]
ax.plot(lags,plotvar,label=expnames_long[ex],c=expcols[ex],ls=expls[ex],marker=emks[ex])


ax.legend()


ax.set_xlabel("<-- [$h_e - h_w$] Leads | SST Leads -->")
ax.set_ylabel(r"Regression Coefficient [$m$ per $\degree C]$")


ax.set_xlim([-24,24])
ax.set_xticks(np.arange(-24,25,2))
#ax.set_ylim([-.015,.015])

ax.axhline([0],ls='solid',lw=0.75,c="k")
ax.axvline([0],ls='solid',lw=0.75,c="k")

ax.grid(True,ls='dotted',c='gray',lw=0.55)


figname = "%sResolution_Significance_ThermoclineDiff_v_Nino3_%s_v_1950_reg%s.png" % (figpath,comparison,regnames_short[regid])
plt.savefig(figname)

plt.show()


# ================================
#%% Part (3), Heat Fluxes
# ================================



vnames = ["str","ssr"]#E"tx_sur"

ds_var_cep   = []
# = []

dsvarflx = []
for vname in vnames:
    ds_byvar = []
    for ex in tqdm.tqdm(range(nexps)):
        
        if vname == "tx_sur":
            vname_file = "tx_surf_1m"
            ncname = "%s%s_%s.nc" % (datpath,expnames[ex],vname_file)
        else:
            ncname = "%s%s_%s.nc" % (datpath,expnames[ex],vname)
        
        ds = xr.open_dataset(ncname)
       
        
        if vname.upper() in list(ds.keys()):
            print("Renaming %s to %s" % (vname.upper(),vname))
            ds = ds.rename({vname.upper(): vname})
        ds = standardize_names(ds)
         
        ds1 = proc.sel_region_xr(ds,bbox_nino3).load()
        ds_byvar.append(ds1[vname] * conversion)
    
    dsvarflx.append(ds_byvar)
    
    # ds_var_cep.append(ds1)
    
    # ds2 = proc.sel_region_xr(ds,bbox_nino3).load()
    
    # ds_var_nino3.append(ds2)

# Preprocess
flx_anoms = []
for ds_var in [dsvarflx[0],dsvarflx[1]]:
    
    #ds_var      = [standardize_names(ds) for ds in ds_var]
    
    if vname in ["D20","Dmaxgrad"]:
        ds_anoms = [preprocess_enso(ds['nz1']) for ds in ds_var]
    else:
        ds_anoms = [preprocess_enso(ds) for ds in ds_var]
    
    flx_anoms.append(ds_anoms)

# Take Area Averages
flx_aavgs = []
for dsreg in flx_anoms:
    aavgs = [proc.area_avg_cosweight(ds) for ds in dsreg]
    flx_aavgs.append(aavgs)
    
#%% Compute the Regressions

#lags = 
beta_byvar = []
for regid in range(2):
    beta_byexp = []
    for ex in range(nexps):
        
        sst_in = ensoids[ex]
        tau_in = flx_aavgs[regid][ex]
        
        sst_in,tau_in = proc.match_time_month(sst_in,tau_in)
        
        betas  = calc_lag_regression_1d(tau_in,sst_in,lags)
        
        beta_byexp.append(betas)
    beta_byvar.append(beta_byexp)
    
    
#%%




#%% For each case, do some calculations

for v in range(2):
    
    fig,ax  = plt.subplots(1,1,figsize=(16,4.5),constrained_layout=True)
    
    for ex in range(nexps):
        
        
        plotvar = beta_byvar[v][ex]
        ax.plot(lags,plotvar,label=expnames_long[ex],c=expcols[ex],ls=expls[ex],marker=emks[ex])
    
    ax.legend()
    
    ax.set_xlim([-24,24])
    ax.set_xticks(np.arange(-24,25,2))
    #ax.set_ylim([-.015,.015])
    
    ax.axhline([0],ls='solid',lw=0.75,c="k")
    ax.axvline([0],ls='solid',lw=0.75,c="k")
    
    if v == 0:
        ax.set_ylim([-10,10])
    else:
        ax.set_ylim([-18,18])
    
    ax.grid(True,ls='dotted',c='gray',lw=0.55)
    #ax.set_title(regnames[regid])
    
    ax.set_xlabel("<-- [$%s$] Leads | SST Leads -->" % vnames[v]) 
    
    ax.set_ylabel(r"Regression Coefficient [$\frac{W}{m^2}$ per $\degree C]$")
    
    figname = "%sLeadLag_Regression_%s_v_Nino3_SST.png" % (figpath,vnames[v])
    plt.savefig(figname)
    
    
    plt.show()



#%% Evaluate significance of such changes


comparison = "ssps85"#"ssps85"
varid      = 0

for varid in range(2):
    for comparison in ['ssps85','Control']:
        # Function Arguments
        if comparison == "Control":
            sample_len        = 240 # months
            ts_full           = ensoids[0].data
            
            target_timeseries = [flx_aavgs[varid][0].data,] # Additional Timeseries to Sample, must be same time as ts_full#[]
        elif comparison == "ssps85":
            sample_len        = 120 # months
            #ts_full           = ensoids[1].data
            ts_full,thermo_in = proc.match_time_month( ensoids[1],flx_aavgs[varid][1])
            ts_full = ts_full.data
            target_timeseries = [thermo_in.data,] # Additional Timeseries to Sample, must be same time as ts_full#[]
        preserve_month    = False
        scramble_year     = False
        mciter            = 1000
                    
        # Test the three methods
        mciter   = 1000
        outdict1 = mcsampler(ts_full,sample_len,mciter,target_timeseries=target_timeseries,preserve_month=False)
        outdict2 = mcsampler(ts_full,sample_len,mciter,target_timeseries=target_timeseries,preserve_month=True,scramble_year=False)
        outdict3 = mcsampler(ts_full,sample_len,mciter,target_timeseries=target_timeseries,preserve_month=True,scramble_year=True)
        outdicts = [outdict1,outdict2,outdict3]
        
        out_result = []
        for ii in range(3):
            
            indict = outdicts[ii]
            
            ts1    = indict['samples'] # iter x time
            ts2    = indict['other_sampled_timeseries'][0] # iter x time
            
            outval = []
            for mc in tqdm.tqdm(range(mciter)):
                betas  = calc_lag_regression_1d(ts2[mc,:],ts1[mc,:],lags)
                outval.append(betas)
                
            outval = np.array(outval)
            out_result.append(outval)
        
        
        #%% Plot significance od ifferent changes
        fig,ax  = plt.subplots(1,1,figsize=(10,4.5),constrained_layout=True)
        
        
        outcols  = ['red','blue','yellow']
        outnames = ["Completely Random","Preserve Month","Preserve Month + Scramble Years"]
        
        for ii in range(3):
            plotvar = out_result[ii]
            
            mu    = np.nanmean(plotvar,0)
            sigma = np.nanstd(plotvar,0)
            
            #clow,chi = np.percentile(plotvar,[0.025,0.975],axis=0)
            #ax.fill_between(lags,clow,chi,color=outcols[ii],alpha=0.15)
            ax.plot(lags,mu,color=outcols[ii],ls='dashed',lw=0.75)
            ax.fill_between(lags,mu-sigma,mu+sigma,color=outcols[ii],alpha=0.15,label=outnames[ii])
            
        if comparison == "Control":
        
            ex    = 0
        elif comparison == "ssps85":
        
            ex    = 1
            
        plotvar = beta_byvar[varid][ex]
        ax.plot(lags,plotvar,label=expnames_long[ex],c=expcols[ex],ls=expls[ex],marker=emks[ex])
        
        if comparison == "Control":
        
            ex    = 2
        elif comparison == "ssps85":
        
            ex    = 3
            
        
        plotvar = beta_byvar[varid][ex]
        ax.plot(lags,plotvar,label=expnames_long[ex],c=expcols[ex],ls=expls[ex],marker=emks[ex])
        
        
        ax.legend()
        
        
        ax.set_xlabel("<-- [$%s$] Leads | SST Leads -->" % vnames[varid])
        ax.set_ylabel(r"Regression Coefficient [$\frac{W}{m^2}$ per $\degree C]$")
        
        
        ax.set_xlim([-24,24])
        ax.set_xticks(np.arange(-24,25,2))
        #ax.set_ylim([-.015,.015])
        
        if varid == 0:
            ax.set_ylim([-10,10])
        else:
            ax.set_ylim([-18,18])
        
        ax.axhline([0],ls='solid',lw=0.75,c="k")
        ax.axvline([0],ls='solid',lw=0.75,c="k")
        
        ax.grid(True,ls='dotted',c='gray',lw=0.55)
        
        
        figname = "%sResolution_Significance_%s_v_Nino3_%s_v_1950.png" % (figpath,vnames[varid],comparison)
        plt.savefig(figname)
        
        #plt.show()
        
        
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
    elif vname == "J/m2"
    
    
    
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

#%% Sanity Check on how polyfit works

x    = np.arange(1000)
test = 2.5*x + -1 + np.random.normal(0,1.25,len(x))

out = np.polyfit(x,test,1)
#%%

#lagout = proc.calc_lag_covar_ann()





#%% 10.10.2025 General Formulation: correlate (nino index) with [vname]




