#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Global Mean and Nino3.4 Regression for TOA Fluxes

Goal: Trying to reproduce the results from the 
Ceppi and Fueglistaler 2021 and Hanke and Proistosescu 2024 papers
but for the AWI-CM3 simulations

Created on Wed Oct 15 11:34:32 2025

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

# %% Niu Paths for Custom Modules

scmpath = "/home/niu4/gliu8/scripts/commons/stochmod/model"
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)

sys.path.append(scmpath)

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)

from amv import proc, viz
import utils as ut
import scm
# %% Other User Edits

gmeanpath = "/home/niu4/gliu8/projects/scrap/global_mean/"
ninopath = "/home/niu4/gliu8/projects/scrap/nino34/"
figpath = "/home/niu4/gliu8/figures/bydate/2025-11-04/"
proc.makedir(figpath)

vnames = ["ttr", "ttrc", "tsr", "tsrc"]

standardize = False
ensoid_name = "nino34"

# Simulation Names -----
expnames = ["TCo319_ctl1950d", "TCo319_ssp585",
            "TCo1279-DART-1950", "TCo1279-DART-2090", "TCo2559-DART-1950C"]
expnames_long = ["31km Control", "31km SSP585",
                 "9km 1950", "9km 2090", "5km 1950"]
expcols = ["cornflowerblue", 'lightcoral',
           "slateblue", "firebrick",
           "midnightblue", "k"]  # Includes Glorys and different shade based on resolution


vnames_other = ['sst','eis',]#'w','r']


# %% Load ENSO Indices


ensoids = [ut.load_ensoid(expname, ensoid_name,
                          standardize=standardize) for expname in expnames]

# %%
loadmaskeis=True
apply_conversion = True
nexps = len(expnames)
vars_byexp = []
skipexps = []
for ex in range(nexps):

    ds_byvar = []

    for v in range(4):
        vname = vnames[v]

        try:
            ncname = "%s%s_%s_global_mean.nc" % (
                gmeanpath, expnames[ex], vname)
            ds = xr.open_dataset(ncname).load()[vname]
            print("Found %s for %s" % (vname, expnames[ex]))
            # Standardize Names
            
            ds = ut.standardize_names(ds)
            ds = ut.remove_duplicate_times(ds)
            if apply_conversion:
                ds = ut.varcheck(ds, vname, expnames[ex])

            ds_byvar.append(ds.copy())
        except:
            print("Could not find %s for %s" % (vname, expnames[ex]))
            if vname == 'eis' and loadmaskeis:
                ncname = "%s%s_%s_global_mean_masked.nc" % (
                    gmeanpath, expnames[ex], vname)
            else:
                ncname = "%s%s_%s_global_mean.nc" % (
                    gmeanpath, expnames[ex], vname)
                

            
            print("\t" + ncname)
            ds = None
            #ds_byvar.append(None)
            skipexps.append(ex)
    vars_byexp.append(ds_byvar)


# Check the shape
for ex in range(nexps):
    print(expnames_long[ex])
    if ex in skipexps:
        continue
    # if None not in vars_byexp[ex]:
    [print("\t %s" % (str(ds.shape))) for ds in vars_byexp[ex]]
    print("\n")

# Merge Variables into datasets
for ex in range(nexps): 
    if ex in skipexps:
        continue
    else:
        vars_byexp[ex] = xr.merge(vars_byexp[ex])

# %% Plot the values

fig, axs = plt.subplots(4, 1, figsize=(12.5, 4.5), constrained_layout=True)

for v in range(4):

    ax = axs[v]
    vname = vnames[v]
    for ex in range(nexps):

        if ex not in skipexps:
            plotvar = vars_byexp[ex][vname].squeeze()
            ax.plot(plotvar.time, plotvar, lw=2, c=expcols[ex])
    ax.set_title(vnames[v])
    ax.set_xlim([plotvar.time.sel(time='1950-01-01', method='nearest'),
                 plotvar.time.sel(time="1960-01-01", method='nearest')])

figname = "%sGlobal_Mean_TOA_Fluxes.png" % figpath
plt.savefig(figname, dpi=150, bbox_inches='tight')

plt.show()

#%% Plot Mean Seasonal Cycle

fig,axs = viz.init_monplot(4,1,figsize=(12,10))
mons3   = proc.get_monstr()

for v in range(4):

    ax = axs[v]
    vname = vnames[v]
    for ex in range(nexps):

        if ex not in skipexps:
            plotvar = vars_byexp[ex][vname].groupby('time.month').mean('time').squeeze()
            ax.plot(mons3, plotvar, lw=2, c=expcols[ex])
    ax.set_title(vnames[v])
            
    
plt.show()

# %% Calculate CRE and Net

cre_lw    = []
cre_sw    = []
net       = []
net_clear = []
for ex in range(nexps):

    if ex in skipexps:
        cre_lw.append(None)
        cre_sw.append(None)
        net.append(None)
        net_clear.append(None)
        continue
    
    dsin = vars_byexp[ex]
    cl   = dsin.ttr - dsin.ttrc     # CRE  LW
    cs   = dsin.tsr - dsin.tsrc     # CRE  SW
    nn   = dsin.tsr + dsin.ttr      # All Sky (Net)
    nc   = dsin.tsrc + dsin.ttrc    # Clear Sky
    
    cre_lw.append(cl)
    cre_sw.append(cs)
    net.append(nn)
    net_clear.append(nc)

net_cre     = [net[ex] - net_clear[ex] if ex not in skipexps else None for ex in range(nexps)]
net_cre_add = [cre_sw[ex] + cre_lw[ex] if ex not in skipexps else None for ex in range(nexps)]


#%% For each experiment, plot the values

plotvars    = [net,net_clear,net_cre,cre_lw,cre_sw]
plotnames   = ["all sky","clear sky","cre","cre (LW)","cre (SW)"]

ccols       = ['k','cornflowerblue','gray','red','blue']
ccls        = ['solid','dashed','dotted','dotted','dotted']

fig,axs     = viz.init_monplot(1,5,figsize=(16,4.5))
mons3       = proc.get_monstr()

for ex in range(nexps):
    ax = axs[ex]
    for v in range(5):
        
        plotvar = plotvars[v][ex].groupby('time.month').mean('time').squeeze()
        ax.plot(mons3, plotvar, lw=2, c=ccols[v],ls=ccls[v],label=plotnames[v])
    ax.set_title(expnames_long[ex])
    
    if ex == 0:
        ax.legend()
    
    
figname = "%sNet_TOA_Fluxes_Scycle.png" % figpath
plt.savefig(figname, dpi=150, bbox_inches='tight')

plt.show()

# %% Perform Lag Regressions

fluxes_in = [net,net_clear,net_cre,cre_sw,cre_lw]
flxnames  = ["All Sky","Clear Sky","CRE","CRE (SW)","CRE (LW)"]

#%% Compute Lag Regressions

preproc  = True
lags     = np.arange(-18,19)
betas_all = np.zeros((len(fluxes_in),nexps,len(lags))) * np.nan
betas_all_mon = np.zeros((len(fluxes_in),nexps,12,len(lags))) * np.nan

# Do additional subsampling (currently only running for historical)
mciter    = 1000
minlen    = 8*12
mcbetas   = np.zeros((len(fluxes_in),mciter,len(lags))) * np.nan

for v in range(len(fluxes_in)):
    for ex in range(nexps):
        
        sst_in            = ensoids[ex].squeeze()
        try:
            flx_in        = fluxes_in[v][ex].squeeze()
            if preproc:
                flx_in = ut.preprocess_enso(flx_in)
            
        except:
            continue
        
        sst_in,flx_in       = proc.match_time_month(sst_in,flx_in)
        #flx_in,sst_in       = proc.match_time_month(sst_in,flx_in)
        #betas               = ut.calc_lag_regression_1d(flx_in.data,sst_in.data,lags)
        #betas_all[v,ex,:] = betas
        
        dsout             = ut.calc_leadlag_regression_2d(sst_in,flx_in.drop(('lat','lon')).squeeze(),lags,sep_mon=False)
        betas_all[v,ex,:] = dsout.regression_coefficient.data
        
        dsoutmon          = ut.calc_leadlag_regression_2d(sst_in,flx_in.drop(('lat','lon')).squeeze(),lags,sep_mon=True)
        betas_all_mon[v,ex,:,:] = dsoutmon.regression_coefficient.data.T
        
        
        
        # Do some additional sampling
        if ex == 0:
            
            outdict_31km = ut.mcsampler(sst_in,minlen,mciter,target_timeseries=[flx_in,],preserve_month=True,scramble_year=False)
            ts1 = outdict_31km['samples']
            ts2 = outdict_31km['other_sampled_timeseries'][0]
            for mc in tqdm.tqdm(range(mciter)):
                bb              = ut.calc_lag_regression_1d(ts2[mc,:],ts1[mc,:],lags)
                mcbetas[v,mc,:] = bb

            print(v)



#%%

fig,axs = plt.subplots(1,5,figsize=(20,3),constrained_layout=True)

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
    ax.set_ylim([-.6,.6])
    ax.tick_params(labelsize=10)

figname = "%sFluxes_Lag_Regression_%s.png" % (figpath,ensoid_name)
plt.savefig(figname, dpi=150, bbox_inches='tight')
plt.show()

#%% Print the Maximum Values

for ii in range(5):
    print(flxnames[ii])
    for ex in range(nexps):
        plotvar = betas_all[ii,ex,:]
        ilag    = np.argmax(np.abs(plotvar))
        maxval  = plotvar[ilag]
        print("\t%.02f @ Lag %02i (%s)" % (maxval,lags[ilag],expnames_long[ex]))
    print("\n")
    
#%% Try computing the monthly feedback (note it doesnt look good... very noiy. Need to think more carefully)


#%% Plot to see if there are any unique monthly features

v     = 0
ex    = 0
vmax  = .75
mons3 = proc.get_monstr()

fig,axs = plt.subplots(5,1,constrained_layout=True,figsize=(12.5,8))


for ex in range(5):
    ax      = axs[ex]
    plotvar = betas_all_mon[v,ex,:,:]
    pcm     = ax.pcolormesh(lags,mons3,plotvar,cmap='cmo.balance',vmin=-vmax,vmax=vmax)

cb = viz.hcbar(pcm,ax=axs.flatten())

plt.show()

#%% Try Plotting Separately for each month

v       = 0



fig,axs = plt.subplots(1,12,constrained_layout=True,figsize=(12.5,4),sharey=True)
for im in range(12):
    ax = axs[im]
    for ex in range(5):
        plotvar = betas_all_mon[v,ex,im,:]
        ax.plot(lags,plotvar,label=expnames_long[ex],c=expcols[ex],lw=2)
    
    ax.set_ylim([-.75,.75])
    
plt.show()



#%%

# Note seems to roughly work... (Note Moved over to utils folder) 
def calc_leadlag_regression_2d(ensoid,dsvar,leadlags,sep_mon=False):
    adddim=False
    if len(dsvar.shape) == 1:
        dummylatlon = dict(lat=[1,],lon=[1,])
        dsvar = dsvar.expand_dims(dim=dummylatlon,axis=(0,1))
        adddim=True
    
    if dsvar.name is None:
        dsvar = dsvar.rename("regression_coefficient")
    
    # Check to make sure the time matches
    dsvar,ensoid    = proc.match_time_month(dsvar,ensoid)
    
    # Get Dimension Lengths
    dsvar           = dsvar.transpose('lon','lat','time')
    nlon,nlat,ntime = dsvar.shape

    if sep_mon is False: # Do for all months
        # Do the Leads (variable leads)
        leads       = np.abs(leadlags[leadlags <=0])
        nleads      = len(leads)
        beta_leads  = np.zeros((nlon,nlat,nleads)) * np.nan
        sig_leads   = beta_leads.copy()
        for ll in range(nleads):
            lag                = leads[ll] 
            ints               = ensoid.data[lag:]
            invar              = dsvar.data[:,:,:(ntime-lag)]
            rout               = proc.regress_ttest(invar,ints,verbose=False)
            beta_leads[:,:,ll] = rout['regression_coeff']
            sig_leads[:,:,ll]  = rout['sigmask']
            
        # Do the lags
        lags        = leadlags[leadlags > 0]
        nlags       = len(lags)
        beta_lags   = np.zeros((nlon,nlat,nlags)) * np.nan
        sig_lags    = beta_lags.copy()
        for ll in range(nlags):
            lag   = lags[ll] 
            ints  = ensoid.data[:(ntime-lag)]
            invar = dsvar.data[:,:,lag:]
            rout  = proc.regress_ttest(invar,ints,verbose=False)
            beta_lags[:,:,ll] = rout['regression_coeff']
            sig_lags[:,:,ll]  = rout['sigmask']
        
        # Concatenate
        betas = np.concatenate([beta_leads,beta_lags],axis=2)
        sigs  = np.concatenate([sig_leads,sig_lags],axis=2)
        
        # Replace into DataArray
        coords   = dict(lon=dsvar.lon,lat=dsvar.lat,lag=leadlags)
        da_betas = xr.DataArray(betas,coords=coords,dims=coords,name=dsvar.name)
        da_sigs  = xr.DataArray(sigs,coords=coords,dims=coords,name="sig")
        da_out   = [ds.transpose('lag','lat','lon') for ds in [da_betas,da_sigs]]
        
    else: # Calculate Separately by Month
        
        # Reshape the variables
        nyr          = int(ntime/12)
        ints_yrmon   = ensoid.data.reshape(nyr,12)
        invar_yrmon  = dsvar.data.reshape(nlon,nlat,nyr,12)
        
        # Calculate Leads (a bit silly to write it this way I know... should prob make a fuction)
        leads       = np.abs(leadlags[leadlags <=0])
        nleads      = len(leads)
        beta_leads  = np.zeros((12,nlon,nlat,nleads)) * np.nan
        sig_leads   = beta_leads.copy()
        for ll in range(nleads):
            lag                = leads[ll]
            if (lag >= nyr):
                print("Cannot perform calculation since the lag (%i) exceeds the # of years %02i. Skipping" % (lag,nyr))
                continue
            for im in range(12):
                ints                    = ints_yrmon[lag:,im]
                invar                   = invar_yrmon[:,:,:(nyr-lag),im] #dsvar.data[:,:,:(ntime-lag)]
                rout                    = proc.regress_ttest(invar,ints,verbose=False)
                beta_leads[im,:,:,ll]   = rout['regression_coeff']
                sig_leads[im,:,:,ll]    = rout['sigmask']
                
        # Calculate Lags
        lags        = leadlags[leadlags > 0]
        nlags       = len(lags)
        beta_lags   = np.zeros((12,nlon,nlat,nlags)) * np.nan
        sig_lags    = beta_lags.copy()
        for ll in range(nlags):
            lag   = lags[ll] 
            if (lag >= nyr):
                print("Cannot perform calculation since the lag (%i) exceeds the # of years %02i. Skipping" % (lag,nyr))
                continue
            for im in range(12):
                ints  = ints_yrmon[:(nyr-lag),im]
                invar = invar_yrmon[:,:,lag:,im]
                rout  = proc.regress_ttest(invar,ints,verbose=False)
                beta_lags[im,:,:,ll] = rout['regression_coeff']
                sig_lags[im,:,:,ll]  = rout['sigmask']
                
        # Concatenate
        betas = np.concatenate([beta_leads,beta_lags],axis=3)
        sigs  = np.concatenate([sig_leads,sig_lags],axis=3)
        
        # Replace into DataArray
        coords   = dict(mon=np.arange(1,13,1),lon=dsvar.lon,lat=dsvar.lat,lag=leadlags)
        da_betas = xr.DataArray(betas,coords=coords,dims=coords,name=dsvar.name)
        da_sigs  = xr.DataArray(sigs,coords=coords,dims=coords,name="sig")
        da_out   = [ds.transpose('lag','mon','lat','lon') for ds in [da_betas,da_sigs]]
    
    # Merge Variables
    if adddim:
        da_out = [ds.squeeze() for ds in da_out]
    ds_out   = xr.merge(da_out)
    
    return ds_out




# %% Check the Mean of each timeseries to see conversion (debugging)

inds = [vars_byexp[ex][vname].mean(
    'time') if ex not in skipexps else None for ex in np.arange(nexps)]

inds[0]/inds[2]
inds[0]/inds[3]
inds[0]/inds[4]

# --------------------------------------------
#%% Examine Relationship with other variables
# --------------------------------------------

othervars = []
for vname in vnames_other:
    print(vname)
    if vname in ['w','r','t']:
        vname_in = vname+"700"
    else:
        vname_in = vname
        
    expvars = []
    for ex in tqdm.tqdm(range(nexps)):
        expname = expnames[ex]
        ncname  = "%s%s_%s_global_mean.nc" % (gmeanpath,expname,vname_in)
        try:
            ds      = xr.open_dataset(ncname)[vname].load()
            ds      = ut.standardize_names(ds)
            ds      = ut.varcheck(ds,vname,expname)
            expvars.append(ds)
        except:
            print("%s not found for %s, skipping" % (vname,expname))
            expvars.append(None)
    othervars.append(expvars)
        
#%% Calculate Lag Regression (copied from above)

preproc  = True

betas_all_othervars = np.zeros((len(othervars),nexps,len(lags))) * np.nan

# Do additional subsampling (currently only running for historical)
mciter    = 1000
minlen    = 8*12
mcbetas   = np.zeros((len(fluxes_in),mciter,len(lags))) * np.nan

for v in range(len(othervars)):
    for ex in range(nexps):
        print(vnames_other[v])
        vname = vnames_other[v]
        
        sst_in            = ensoids[ex].squeeze()
        try:
            flx_in        = othervars[v][ex].squeeze()
            if preproc:
                flx_in = ut.preprocess_enso(flx_in.squeeze())
            
        except:
            print("Issue on %s for %s, skipping" % (vname,expname))
            continue
        
        sst_in,flx_in        = proc.match_time_month(sst_in,flx_in)
        #flx_in,sst_in       = proc.match_time_month(sst_in,flx_in)
        #betas               = ut.calc_lag_regression_1d(flx_in.data,sst_in.data,lags)
        #betas_all[v,ex,:] = betas
        
        dsout                       = ut.calc_leadlag_regression_2d(sst_in,flx_in.drop(('lat','lon')).squeeze(),lags,sep_mon=False)
        betas_all_othervars[v,ex,:] = dsout[vname].data#.regression_coefficient.data
        
        
        
        # Do some additional sampling
        if ex == 0:
            
            outdict_31km = ut.mcsampler(sst_in,minlen,mciter,target_timeseries=[flx_in,],preserve_month=True,scramble_year=False)
            ts1 = outdict_31km['samples']
            ts2 = outdict_31km['other_sampled_timeseries'][0]
            for mc in tqdm.tqdm(range(mciter)):
                bb              = ut.calc_lag_regression_1d(ts2[mc,:],ts1[mc,:],lags)
                mcbetas[v,mc,:] = bb

            print(v)


#%% For each variable, plot the relationship with Nino 3.4

#for v in range(len(othervars)):

v = 1

fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(4.5,4.5))

for ex in range(nexps):
    
    plotvar = betas_all_othervars[v,ex,:]
    ax.plot(lags,plotvar,label=expnames_long[ex],color=expcols[ex],lw=2.5)

ax.set_xticks(lags[::2])
ax.set_xlim([lags[0],lags[-1]])

ax.axhline([0],lw=0.55,c='k')
ax.axvline([0],lw=0.55,c='k')
ax.set_title(vnames_other[v])
    

ax.set_ylim([-0.5,0.5])
ax.legend()
plt.show()

#%% Part III: Do Composites around ENSO events to understand their progression








# %% Load Global Means
