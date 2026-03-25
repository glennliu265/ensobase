#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

On 1x1 Regridded data, try to compute the SW feedback

Notes
Need to move over the lead lag code
Need to better rewrite xrfunc version



Created on Mon Feb 23 11:31:28 2026

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

#%% Indicate Paths

# Get the components
expnames = ["TCo319_ctl1950d", "TCo319_ssp585", "TCo1279-DART-1950", "TCo1279-DART-2090", "TCo2559-DART-1950C"]
dpath    = "/home/niu4/gliu8/projects/scrap/regrid_1x1/"
vnames   = ["sst","tscre"]

figpath      = "/home/niu4/gliu8/figures/bydate/2026-03-03/"
proc.makedir(figpath)

bbox_sep = [240,290,-50,10]
bbox     = bbox_sep
selmon   = [12,1,2]


#%% Load the variables

nexp   = len(expnames)
nvar   = len(vnames)
ds_all = []
for ex in tqdm.tqdm(range(nexp)):
    
    ds_exp = []
    for vv in range(nvar):
        
        ncname = "%s%s_%s_regrid1x1.nc" % (dpath,expnames[ex],vnames[vv])
        ds     = xr.open_dataset(ncname)[vnames[vv]].load()
        
        ds     = ut.standardize_names(ds)
        ds     = ut.remove_duplicate_times(ds)
        ds     = ut.varcheck(ds,vnames[vv],expnames[ex])
        
        ds     = proc.sel_region_xr(ds,bbox_sep)
        
        ds_exp.append(ds)
    
    #if expname != "TCo2259-DART-1950":
    #ds_exp = xr.merge(ds_exp)
    ds_all.append(ds_exp)


#%% Detrend, Deseason


sstanom_all   = [proc.xrdeseason(ds[0]) for ds in ds_all]
tscreanom_all = [proc.xrdeseason(ds[1]) for ds in ds_all]

#dsanom_all    = [proc.xrdeseason(ds) for ds in ds_all]

sst_dt         = [proc.xrdetrend(ds) for ds in sstanom_all]
tscre_dt      = [proc.xrdetrend(ds) for ds in tscreanom_all]

#%% Compute Pointwise Lag correlation



# def calc_lag_regression_1d(var1_lag,var2_base,lags,correlation=False): # CAn make 2d by mirroring calc_lag_covar_annn
#     # Calculate the regression where
#     # (+) lags indicate var1 lags  var2 (var 2 leads)
#     # (-) lags indicate var1 leads var2 (var 1 leads)
    
#     if np.any(np.isnan(var1_lag)) or np.any(np.isnan(var2_base)):
#         #print("NaN detected. Returning NaN...")
#         return np.nan * np.ones(len(lags*2))
    
#     ntime = len(var1_lag)
#     betalag = []
#     poslags = lags[lags >= 0]
#     for l,lag in enumerate(poslags):
#         varlag   = var1_lag[lag:]
#         varbase  = var2_base[:(ntime-lag)]
        
#         # Calculate correlation
#         if correlation:
#             beta = np.corrcoef(varbase,varlag)[0,1]
#         else:
#             beta = np.polyfit(varbase,varlag,1)[0]   
#         betalag.append(beta.item())
    
#     neglags = lags[lags < 0]
#     neglags_sort = np.sort(np.abs(neglags)) # Sort from least to greatest #.sort
#     betalead = []
    
#     for l,lag in enumerate(neglags_sort):
#         varlag   = var2_base[lag:] # Now Varlag is the base...
#         varbase  = var1_lag[:(ntime-lag)]
#         # Calculate correlation
#         if correlation:
#             beta = np.corrcoef(varlag,varbase)[0,1]
#         else:
#             beta = np.polyfit(varlag,varbase,1)[0]   
#         betalead.append(beta.item())
    
#     # Append Together
#     return np.concatenate([np.flip(np.array(betalead)),np.array(betalag)])


# # Copied from xrfunc in amv (based on leadlagcorr.
# def leadlagreg(dslag,dsbase,leadlags):
#     # Computes lead lag correlation between ds1 and ds2 over dimension ['time']
#     # over specified lead/lags [lags]. Loops over all other dimensions (lat,lon,depth,etc)
#     #calc_leadlag    = lambda x,y: proc.leadlag_corr(x,y,lags,corr_only=True)
#     # Note: I moved this to utilities module
    
#     def calc_leadlag_regr(varlag,varbase):
#         try:
#             llout = calc_lag_regression_1d(varlag,varbase,leadlags)
#             return llout
#         except:
#             return np.ones(len(leadlags)) * np.nan
    
#     llreg = xr.apply_ufunc(
#         calc_leadlag_regr,
#         dslag,
#         dsbase,
#         input_core_dims=[['time'],['time']],
#         output_core_dims=[['lags']],
#         vectorize=True,
#         )
#     #leadlags     = np.concatenate([np.flip(-1*lags)[:-1],lags],) 
#     llreg['lags'] = leadlags
#     return llreg

zero2nan       = lambda ds: xr.where(ds == 0.,np.nan,ds)
leadlags       = np.arange(-12,12)
leadlags_byexp = []
leadlags_byexp2 = []

for ex in tqdm.tqdm(range(nexp)):
    sst       = sst_dt[ex]
    tscre     = tscre_dt[ex]
    
    sst       = zero2nan(sst)
    tscre     = zero2nan(tscre)
    
    sst,tscre = proc.match_time_month(sst,tscre)
    
    sst['time'] = tscre['time']
    
    
    #llreg     = loop_leadlag(sst,tscre)
    #llreg    = leadlagreg(tscre,sst,leadlags)
    
    llreg = ut.calc_leadlagreg_pointwise(tscre,sst,leadlags)
    #leadlags_byexp2.append(llreg)
    
    leadlags_byexp.append(llreg)


#%% Only the Last One Didn't Work...
# Temp Save

outpathtemp = "/home/niu4/gliu8/projects/scrap/scratch/"
for ex in tqdm.tqdm(range(nexp)): 
    outname = "%s%s_sst_cre_leadlag_corr_12.nc" % (outpathtemp,expnames[ex])
    leadlags_byexp[ex].to_netcdf(outname)

#%% Look at the values, globally

#for ex in range(nexp):

ex      = 3

fig,ax  = ut.init_globalmap(1,1,)
plotvar = leadlags_byexp[ex].sel(lags=0)
pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,cmap='cmo.balance')
cb      = viz.hcbar(pcm,ax=ax)

plt.show()

#%% Download and Check Reference Values for CERES

obspath   = "/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_EBAF_ERA5_2001_2024/raw/"
cerespath = "/home/niu4/gliu8/share/CERES/processed/"
flxobs    = xr.open_dataset(obspath + "tscre.nc").tscre.load()
sflxobs   = xr.open_dataset(cerespath + "CERES_EBAF_sscre_2000-03_to_2025-09.nc").sscre.load()
sstobs    = xr.open_dataset(obspath + "sst.nc").sst.load()

invar_obs = [sstobs,flxobs,sflxobs]
invar_obs = [ut.standardize_names(ds) for ds in invar_obs]
invar_obs = [proc.sel_region_xr(ds,bbox_sep) for ds in invar_obs]

anom_obs        = [proc.xrdeseason(ds) for ds in invar_obs]
anom_dt_obs     = [proc.xrdetrend(ds) for ds in anom_obs]

sstobs,tscreobs = proc.match_time_month(anom_dt_obs[0],anom_dt_obs[1])
_,sscreobs = proc.match_time_month(anom_dt_obs[0],anom_dt_obs[2])

#%% Newer Calculation

leadlags       = np.arange(-12,12)
ds_llout       = ut.calc_leadlagreg_pointwise(tscreobs,sstobs,leadlags)

#%% Check Seasonality


def preprocess(ds):
    ds  = ut.standardize_names(ds)
    dsa = proc.xrdeseason(ds)
    dsa_dt = proc.xrdetrend(dsa)
    return dsa_dt
           
sstobs_anom   = preprocess(sstobs)
tscreobs_anom = preprocess(sflxobs)



sstobs_anom,tscreobs_anom = proc.match_time_month(sstobs_anom,tscreobs_anom)


sstobs_anom['time'] = tscreobs_anom['time']

#%%

sstobs_djf   = proc.selmon_ds(sstobs_anom,[12,1,2])
tscreobs_djf = proc.selmon_ds(tscreobs_anom,[12,1,2])




leadlags     = np.arange(-6,7,1)#[-6,6]

ds_llout_djf = ut.calc_leadlagreg_pointwise(tscreobs_djf,sstobs_djf,leadlags)


#%%

ds_llout_djf.sel(lags=0).plot(vmin=-15,vmax=15,cmap='cmo.balance'),plt.show()

bbox_mini = [-90+360,-75+360,-40,-15]



#%% Now Check the values

dsout_sep_djf = proc.sel_region_xr(ds_llout_djf,bbox_mini)

aavg_sep_djf = proc.area_avg_cosweight(dsout_sep_djf)


#%% Older Calculation / Loop Version
#llobs           = leadlagreg(sstobs,tscreobs,leadlags)

# leadlags       = np.arange(-12,12)
# ntime,nlat,nlon = sstobs.shape
# llout_all          = np.zeros((len(leadlags),nlat,nlon)) * np.nan
# llout_surf          = np.zeros((len(leadlags),nlat,nlon)) * np.nan
# for o in tqdm.tqdm(range(nlon)):
#     for a in range(nlat):
        
#         # Calculate for TOA
#         sstpt = sstobs.isel(lon=o,lat=a)
#         flxpt = tscreobs.isel(lon=o,lat=a)
#         sflxpt     = sscreobs.isel(lon=o,lat=a)
#         if np.any(np.isnan(sstpt)) or np.any(np.isnan(flxpt)) or np.any(np.isnan(sflxpt)):
#             continue
        
#         #llout = ut.calc_lag_regression_1d(sstpt,flxpt,leadlags)
#         # First variable is lagged!
#         llout = ut.calc_lag_regression_1d(flxpt,sstpt,leadlags)
#         llout_all[:,a,o] = llout.copy()
        
#         # Recalculate for surface too
        
#         #sllout = ut.calc_lag_regression_1d(sstpt,sflxpt,leadlags)
#         #sllout = ut.calc_lag_regression_1d(sflxpt,sstpt,leadlags)
#         #llout_surf[:,a,o] = sllout.copy()
        
        
# coords   = dict(lags=leadlags,lat=sstobs.lat,lon=sstobs.lon)
# ds_llout = xr.DataArray(llout_all,dims=coords,coords=coords,name='beta')
        
        
        
# ds_llout_surf = xr.DataArray(llout_surf,dims=coords,coords=coords,name='beta')
        

# dtmon = 1
# (ds_llout_surf.sel(lags=0)*dtmon).plot(),plt.show()

# (ds_llout_surf.sel(lags=1)*dtmon).plot(),plt.show()
    
#%%

def loop_leadlag(sst,flx):
    
    ntime,nlat,nlon     = sst.shape
    llout_all           = np.zeros((len(leadlags),nlat,nlon)) * np.nan
    llout_surf          = np.zeros((len(leadlags),nlat,nlon)) * np.nan
    for o in tqdm.tqdm(range(nlon)):
        for a in range(nlat):
            
            # Calculate for TOA
            sstpt = sst.isel(lon=o,lat=a)
            flxpt = flx.isel(lon=o,lat=a)
            if np.any(np.isnan(sstpt)) or np.any(np.isnan(flxpt)):
                continue
            
            #llout = ut.calc_lag_regression_1d(sstpt,flxpt,leadlags)
            # First variable is lagged!
            llout = calc_lag_regression_1d(flxpt,sstpt,leadlags)
            llout_all[:,a,o] = llout.copy()
    
    coords   = dict(lags=leadlags,lat=sstobs.lat,lon=sstobs.lon)
    ds_llout = xr.DataArray(llout_all,dims=coords,coords=coords,name='beta')
    return ds_llout

#%% Compute Pointwise lag regression

latf = -20
lonf = -80+360

sstpt   = proc.selpt_ds(sstobs,lonf,latf)
tscrept = proc.selpt_ds(tscreobs,lonf,latf)

llout   = ut.calc_lag_regression_1d(flxpt,sstpt,leadlags)
llout   = ut.calc_lag_regression_1d(sstpt,flxpt,leadlags)

np.polyfit(sstpt,tscrept,deg=1)


ds_llout_surf.sel(lags=0).plot(vmin=-15,vmax=15,cmap='cmo.balance'),plt.show()


#%% Check mean state variable


obsmean    = proc.sel_region_xr(flxobs,bbox_sep).mean('time')
modelmeans = [ds.mean('time').tscre for ds in ds_all]

fig,axs = plt.subplots(1,nexp,figsize=(20,4.5))

for ex in range(nexp):
    ax =axs[ex]
    
    plotvar = modelmeans[ex]
    
    pcm=ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,cmap='cmo.balance',vmin=-100,vmax=100)
    ax.set_title(expnames[ex])
    
cb = viz.hcbar(pcm,ax=axs.flatten())
plt.show()

#%%

fig,ax  = plt.subplots(1,1,figsize=(5,4.5))

plotvar = obsmean

pcm=ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,cmap='cmo.balance',vmin=-100,vmax=100)
ax.set_title("CERES")
cb = viz.hcbar(pcm,ax=ax)
plt.show()



#%% Plot Zero Lead Lag

bbox_mini    = [-90+360,-75+360,-40,-15] # This is reset by plot box, need to fix this
bbox_mini360 = [-90+360,-75+360,-40,-15]

proj    = ccrs.PlateCarree()
cints   = np.arange(-16,18,2)

fig,axs = plt.subplots(1,nexp,figsize=(20,4.5),constrained_layout=True,subplot_kw={'projection':proj})

for ex in range(nexp):
    ax      = axs[ex]
    ax      = viz.add_coast_grid(ax,bbox_sep,fill_color="k")
    
    plotvar = leadlags_byexp[ex].sel(lags=0)
    
    pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,cmap='cmo.balance',vmin=-15,vmax=15)
    cl      = ax.contour(plotvar.lon,plotvar.lat,plotvar,levels=cints,linewidths=0.55,colors="k")
    ax.clabel(cl)
    
    # Calculate Mean Feedback
    meanreg = proc.sel_region_xr(plotvar,bbox_mini360)
    meanreg = proc.area_avg_cosweight(meanreg)
    
    ax.set_title("%s, %.2f W/m2/K" % (expnames[ex],meanreg.data.item()))
    
    viz.plot_box(bbox_mini,ax=ax,color="k",linestyle='dashed',linewidth=0.75)
    
cb = viz.hcbar(pcm,ax=axs.flatten())
cb.set_label("TOA Shortwave Cloud Feedback (Lag 0, SST,TSCRE) $W m^{-2} K^{-1}$")


figname = "%sSEP_Shortwave_Cloud_Feedback_AWI-CM3.png" % (figpath)
plt.savefig(figname,dpi=150,bbox_inches='tight')

plt.show()

#%% Repeat Plot for Observations

fig,ax = plt.subplots(1,1,figsize=(5,4.5),constrained_layout=True,subplot_kw={'projection':proj})

# ----
ax      = viz.add_coast_grid(ax,bbox_sep,fill_color="k")

plotvar = ds_llout.sel(lags=0)

pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,cmap='cmo.balance',vmin=-15,vmax=15)
cl      = ax.contour(plotvar.lon,plotvar.lat,plotvar,levels=cints,linewidths=0.55,colors="k")
ax.clabel(cl)

# Calculate Mean Feedback
meanreg = proc.sel_region_xr(plotvar,bbox_mini360)
meanreg = proc.area_avg_cosweight(meanreg)

ax.set_title("%s, %.2f W/m2/K" % ("CERES EBAF and ERA5",meanreg.data.item()))

viz.plot_box(bbox_mini,ax=ax,color="k",linestyle='dashed',linewidth=0.75)
# ----

cb = viz.hcbar(pcm,ax=ax)
cb.set_label("TOA Shortwave Cloud Feedback (Lag 0, SST,TSCRE) $W m^{-2} K^{-1}$")


figname = "%sSEP_Shortwave_Cloud_Feedback_Obs.png" % (figpath)
plt.savefig(figname,dpi=150,bbox_inches='tight')


#%% Look at some area-average values



obs_llout    = proc.sel_region_xr(ds_llout,bbox_mini360)
model_llouts = [proc.sel_region_xr(ds,bbox_mini360) for ds in leadlags_byexp]


#%% Compare the mean

ytks = np.arange(-4,16,2)
fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(8,4.5))

expcols   = ["deepskyblue","pink","royalblue","red","darkslateblue"]
expmarker = ["v","v","o","o","^"]

for ex in range(nexp):
    expname = expnames[ex]
    
    if "TCo319" in expname:
        ls = "dotted"
        
    elif "TCo1279" in expname:
        ls = "dashed"
    else:
        ls = "solid"
    
    # if "1950" or "ctl" in expname:
    #     c  = "blue"
    # elif "2090" or "ssp585" in expname:
    #     c  = "red"
    
    plotvar = proc.area_avg_cosweight(model_llouts[ex])
    ax.plot(leadlags,plotvar,label=expname,
            lw=2,c=expcols[ex],ls=ls,marker=expmarker[ex],fillstyle='none')

plotvar = proc.area_avg_cosweight(obs_llout)
ax.plot(leadlags,plotvar,label="CERES-ERA5",color="k",lw=2,marker="d",fillstyle='none')

ax.legend()
ax = viz.add_axlines(ax)
ax.set_xticks(leadlags)
ax.set_xlabel(" <-- Flux Leads <--  | Lags (Month) | --> SST Leads -->")
ax.set_yticks(ytks)
ax.set_ylabel("TOA SW Cloud Feedback (W/m2/K)")
ax.grid(True,ls='dotted')
ax.set_xlim([-11,11])

figname = "%sSEP_Shortwave_Cloud_Feedback_Area_Average_LeadLag.png" % (figpath)
plt.savefig(figname,dpi=150,bbox_inches='tight')
plt.show()

#%% For each case, plot all the points

for ex in range(nexp):
    #ex     = 0
    fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(8,4.5))
    
    
    expname = expnames[ex]
    
    if "TCo319" in expname:
        ls = "dotted"
        
    elif "TCo1279" in expname:
        ls = "dashed"
    else:
        ls = "solid"
    
    incorr        = model_llouts[ex]
    nl,nlatr,nlonr = incorr.shape
    
    for a in range(nlatr):
        for o in range(nlonr):
            plotvar = incorr.isel(lat=a,lon=o)
            ax.plot(leadlags,plotvar,label="",alpha=0.05,
                    lw=2,c='gray',ls=ls)
    
    # Plot Area Average
    plotvar = proc.area_avg_cosweight(model_llouts[ex])
    ax.plot(leadlags,plotvar,label=expname,
            lw=2,c='blue',ls=ls,marker=expmarker[ex],fillstyle='none')
    
    
    # Plot Observationss
    plotvar = proc.area_avg_cosweight(obs_llout)
    ax.plot(leadlags,plotvar,label="CERES-ERA5",color="k",lw=2,marker="d",fillstyle='none')
    
    
    #ax.legend()
    ax.set_title(expnames[ex])
    ax = viz.add_axlines(ax)
    ax.set_xticks(leadlags)
    ax.set_xlabel(" <-- Flux Leads <--  | Lags (Month) | --> SST Leads -->")
    ax.set_yticks(ytks)
    ax.set_ylabel("TOA SW Cloud Feedback (W/m2/K)")
    ax.grid(True,ls='dotted')
    ax.set_xlim([-11,11])
    ax.set_ylim([ytks[0],ytks[-1]])
    
    figname = "%sSEP_Shortwave_Cloud_Feedback_Area_Average_LeadLag_%s.png" % (figpath,expname)
    plt.savefig(figname,dpi=150,bbox_inches='tight')

#%% Repeat same, but for observations

fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(8,4.5))


expname = expnames[ex]

if "TCo319" in expname:
    ls = "dotted"
    
elif "TCo1279" in expname:
    ls = "dashed"
else:
    ls = "solid"

incorr        = obs_llout
nl,nlatr,nlonr = incorr.shape

for a in range(nlatr):
    for o in range(nlonr):
        plotvar = incorr.isel(lat=a,lon=o)
        ax.plot(leadlags,plotvar,label="",alpha=0.05,
                lw=2,c='gray',ls=ls)

# # Plot Area Average
# plotvar = proc.area_avg_cosweight(model_llouts[ex])
# ax.plot(leadlags,plotvar,label=expname,
#         lw=2,c='blue',ls=ls,marker=expmarker[ex],fillstyle='none')


# Plot Observationss
plotvar = proc.area_avg_cosweight(obs_llout)
ax.plot(leadlags,plotvar,label="CERES-ERA5",color="k",lw=2,marker="d",fillstyle='none')


#ax.legend()
ax.set_title("CERES-ERA5")
ax = viz.add_axlines(ax)
ax.set_xticks(leadlags)
ax.set_xlabel(" <-- Flux Leads <--  | Lags (Month) | --> SST Leads -->")
ax.set_yticks(ytks)
ax.set_ylabel("TOA SW Cloud Feedback (W/m2/K)")
ax.grid(True,ls='dotted')
ax.set_xlim([-11,11])
ax.set_ylim([ytks[0],ytks[-1]])

figname = "%sSEP_Shortwave_Cloud_Feedback_Area_Average_LeadLag_%s.png" % (figpath,"Obs")
plt.savefig(figname,dpi=150,bbox_inches='tight')


#%%
# bboxplot  = bbox_sep
# proj      = ccrs.PlateCarree()
# fig,ax,_  = viz.init_orthomap(1,1,bboxplot,centlon=-90,centlat=-10,figsize=(24,12))
# ax        = viz.add_coast_grid(ax,bboxplot)


#%% Compare the calculations

latf = -20
lonf = -80+360


ex = 0
for ex in range(nexp):
    
    loopr  = proc.selpt_ds(leadlags_byexp[ex],lonf,latf)
    ufuncr = proc.selpt_ds(leadlags_byexp2[ex],lonf,latf)
    
    
    
    fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(8,4.5))
    
    
    
    # Plot Observationss
    plotvar = loopr
    ax.plot(leadlags,plotvar,label="Loooped",color="k",lw=2,marker="d",fillstyle='none')
    
    # Plot Observationss
    plotvar = ufuncr
    ax.plot(leadlags,plotvar,label="ufunc",color="r",lw=2,marker="d",fillstyle='none',ls='dashed')
    
    
    #ax.legend()
    ax.set_title("Loop vs Ufunc for exp %s" % expnames[ex])
    ax = viz.add_axlines(ax)
    ax.set_xticks(leadlags)
    ax.set_xlabel(" <-- Flux Leads <--  | Lags (Month) | --> SST Leads -->")
    ax.set_yticks(ytks)
    ax.set_ylabel("TOA SW Cloud Feedback (W/m2/K)")
    ax.grid(True,ls='dotted')
    ax.set_xlim([-11,11])
    ax.set_ylim([ytks[0],ytks[-1]])
    
    plt.show()
    
#%% Try Compute Lead Lag Regression of SEP-Averaged Values

sst_aavg = proc.area_avg_cosweight(proc.sel_region_xr(sstobs,bbox_mini360))
tscre_aavg = proc.area_avg_cosweight(proc.sel_region_xr(tscreobs,bbox_mini360))

aavg_ll = ut.calc_lag_regression_1d(tscre_aavg,sst_aavg,leadlags)
#%% Plot the Comparison




fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(8,4.5))

plotvar = proc.area_avg_cosweight(obs_llout)
ax.plot(leadlags,plotvar,label="Area-Average Last",color="k",lw=2,marker="d",fillstyle='none')

plotvar = aavg_ll
ax.plot(leadlags,plotvar,label="Area-Average First ",color="blue",lw=2,marker="d",fillstyle='none')

ax.legend()
ax = viz.add_axlines(ax)
ax.set_xticks(leadlags)
ax.set_xlabel(" <-- Flux Leads <--  | Lags (Month) | --> SST Leads -->")
ax.set_yticks(ytks)
ax.set_ylabel("TOA SW Cloud Feedback (W/m2/K)")
ax.grid(True,ls='dotted')
ax.set_xlim([-11,11])

figname = "%sSEP_Shortwave_Cloud_Feedback_Area_Average_LeadLag_obs_Comparison.png" % (figpath)
plt.savefig(figname,dpi=150,bbox_inches='tight')
plt.show()
