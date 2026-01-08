#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Analyze the relationship between ENSO and NAO in AWI-CM3

Created on Tue Dec  9 17:19:13 2025

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

#%%

figpath     = "/home/niu4/gliu8/figures/bydate/2025-12-AWI-Hackathon/"
proc.makedir(figpath)

bbox     = [-90,40,20,80] # Assumes with degrees west #[-90+360, 40, 20, 80]
expnames = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090",]#"TCo2559-DART-1950C"]

#%% Load Nino3.4 Indices

standardize = False
ensoid_name = "nino34"
ensoids     = [ut.load_ensoid(expname,ensoid_name,standardize=standardize) for expname in expnames]

#%% All Month


#%% Load NAO EOFs (All Month)

naopath     = "/home/niu4/gliu8/projects/scrap/nao_crop/"
ncname_eof  = naopath + "%s_nao_index.nc"
ds_naos     = [xr.open_dataset(ncname_eof % (expnames[ex])).load() for ex in range(4)]

#%% Look at Scatterplot

naoids = [ds.pcs for ds in ds_naos]

#%% 

for ex in range(4):
    ensoin = ensoids[ex]
    naoin  = naoids[ex].isel(mode=0)
    plt.scatter( ensoin, naoin)
    plt.title("%s, r= %.2f" % (expnames[ex],np.corrcoef(ensoin,naoin)[0,1]))
    plt.show()

#%% Calculate Lead Lag relationship

lags        = np.arange(-60,61,1)
nlags       = len(lags)
llcorr      = np.zeros([4,nlags]) * np.nan

for ex in range(4):
    ensoin = ensoids[ex]
    naoin  = naoids[ex].isel(mode=0)

    llcorr[ex,:] = ut.calc_lag_regression_1d(ensoin.data,naoin.data,lags,correlation=True)
    
#%%

fig,ax = plt.subplots(1,1)

for ex in range(4):
    ax.plot(lags,llcorr[ex,:],label=expnames[ex])
ax.legend()
plt.show()

#%% Subsample from future and past runs

minlen_control = 240
minlen_ssp585  = 120
mciter         = 1000

outdict_control = ut.mcsampler(ensoids[0],minlen_control,mciter,target_timeseries=[naoids[0].isel(mode=0),],preserve_month=True,scramble_year=False)
outdict_ssp585  = ut.mcsampler(ensoids[1],minlen_ssp585,mciter,target_timeseries=[naoids[1].isel(mode=0),],preserve_month=True,scramble_year=False)

llcorr_control = np.zeros((mciter,nlags)) * np.nan
llcorr_ssp585  = np.zeros((mciter,nlags)) * np.nan
for mc in tqdm.tqdm(range(mciter)):
    
    ts1 = outdict_control['samples']
    ts2 = outdict_control['other_sampled_timeseries'][0]
    bb  = ut.calc_lag_regression_1d(ts2[mc,:],ts1[mc,:],lags)
    llcorr_control[mc,:] = bb.copy()
    
    
    ts1 = outdict_ssp585['samples']
    ts2 = outdict_ssp585['other_sampled_timeseries'][0]
    bb  = ut.calc_lag_regression_1d(ts2[mc,:],ts1[mc,:],lags)
    llcorr_ssp585[mc,:] = bb.copy()
    

#%% Plot Control with significance

xtks    = lags[::3]
expcols = ["cornflowerblue","hotpink","mediumblue","firebrick"]
expls   = ['dashed','dashed','solid','solid']


fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(10,4.5))

# Plot the Full Lead-Lag
explot = [0,2]
for ex in explot:
    ax.plot(lags,llcorr[ex,:],label=expnames[ex],c=expcols[ex],ls=expls[ex])
ax.legend()

# Plot the range of laed lag
mu    = llcorr_control.mean(0)
sigma = llcorr_control.std(0)
ax.fill_between(lags,mu-sigma,mu+sigma,color=expcols[0],alpha=0.10,label="")
ax.plot(lags,mu,color='magenta',ls='dotted',lw=0.55,label="")

ax.axhline([0],lw=0.55,c="k")
ax.axvline([0],lw=0.55,c="k")

ax.set_ylim([-.5,.5])

ax.set_xticks(xtks)
ax.set_xlim([lags[0],lags[-1]])
ax.set_ylabel("Correlation")
ax.set_xlabel(" <-- NAO Leads | ENSO Leads -->")

figname = "%sAWI-CM3_NAO_ENSO_LeadLag_Control.png" % (figpath)
plt.savefig(figname,dpi=150)

#plt.show()
    
#%% Plot SSP585 with significance

xtks    = lags[::3]
expcols = ["cornflowerblue","hotpink","mediumblue","firebrick"]
expls   = ['dashed','dashed','solid','solid']

fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(10,4.5))

# Plot the Full Lead-Lag
explot = [1,3]
for ex in explot:
    ax.plot(lags,llcorr[ex,:],label=expnames[ex],c=expcols[ex],ls=expls[ex])
ax.legend()

# Plot the range of laed lag
mu    = llcorr_ssp585.mean(0)
sigma = llcorr_ssp585.std(0)
ax.fill_between(lags,mu-sigma,mu+sigma,color=expcols[0],alpha=0.10,label="")
ax.plot(lags,mu,color='gray',ls='dotted',lw=0.55,label="")

ax.axhline([0],lw=0.55,c="k")
ax.axvline([0],lw=0.55,c="k")

ax.set_ylim([-.5,.5])

ax.set_xticks(xtks)
ax.set_xlim([lags[0],lags[-1]])
ax.set_ylabel("Correlation")
ax.set_xlabel(" <-- NAO Leads | ENSO Leads -->")
ax.legend()

figname = "%sAWI-CM3_NAO_ENSO_LeadLag_SSP585.png" % (figpath)
plt.savefig(figname,dpi=150)

plt.show()

#%% Just Plotting the Indices to examine the timeseries

ex      = 1

for ex in range(4):
    
    fig,axs = plt.subplots(2,1,constrained_layout=True,figsize=(12,6))\
        
    ax      = axs[0]
    plotvar = ensoids[ex]
    ax.plot(plotvar.time,plotvar,color="k")
    
    timeplot = plotvar.time.data
    plotvar  = plotvar.data
    maskneg  = plotvar.squeeze() <0
    maskpos  = plotvar.squeeze() >=0
    #ax.fill_between(timeplot[maskneg],plotvar.squeeze()[maskneg],np.zeros(len(timeplot[maskneg])),label="ENSO-",color='cornflowerblue',alpha=0.25)
    #ax.fill_between(timeplot[maskpos],plotvar.squeeze()[maskpos],np.zeros(len(timeplot[maskpos])),label="ENSO+",color='tomato',alpha=0.25)
    #ax.bar(timeplot[maskneg],plotvar.squeeze()[maskneg],label="ENSO-",color='cornflowerblue',width=5,alpha=1)
    #ax.bar(timeplot[maskpos],plotvar.squeeze()[maskpos],label="ENSO+",color='tomato',width=5,alpha=1)
    
    ax.legend()
    
    #plt.show()
    
    #
    
    ax      = axs[1]
    plotvar = naoids[ex].isel(mode=0)
    ax.plot(plotvar.time,plotvar,color="k")
    
    timeplot = plotvar.time.data
    plotvar  = plotvar.data
    maskneg  = plotvar.squeeze()<0
    maskpos  = plotvar.squeeze()>=0
    #ax.bar(timeplot[maskneg],plotvar.squeeze()[maskneg],label="NAO-",color='cornflowerblue',width=1,alpha=1)
    #ax.bar(timeplot[maskpos],plotvar.squeeze()[maskpos],label="NAO+",color='tomato',width=1,alpha=1)
    
    ax.legend()
    
    for ax in axs:
        ax.axhline([0],lw=0.55,c="k")
        ax.set_xlim([timeplot[0],timeplot[-1]])
        ax.set_ylim([-4,4])
    plt.suptitle("%s, r=%.02i" % (expnames[ex],np.corrcoef(ensoids[ex],naoids[ex].isel(mode=0))[0,1]))
    figname = "%sAWI-CM3_NAO_ENSO_Timeseries_%s.png" % (figpath,expnames[ex])
    plt.savefig(figname,dpi=150)

plt.show()

#%% Visualize NAO Pattern

cints    = np.arange(-600,650,50)
proj     = ccrs.PlateCarree()

eof1     = [ds.eofs.isel(mode=0) for ds in ds_naos]
eof2     = [ds.eofs.isel(mode=1) for ds in ds_naos]
eof3     = [ds.eofs.isel(mode=2) for ds in ds_naos]

#%% Plot EOF1

bboxplot = [-90,40,20,80]

for ex in range(4):
    fig,ax,_ = viz.init_orthomap(1,1,bboxplot,centlon=-30)
    
    ax       = viz.add_coast_grid(ax,bboxplot)
    
    
    plotvar = eof1[ex]
    if ex == 3:
        plotvar = plotvar * -1
    
    pcm = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,cmap='cmo.balance',vmin=cints[0],vmax=cints[-1])
    ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,colors="k",levels=cints,linewidths=0.55)
    
    ax.set_title("%s EOF1, Var. Explained = %.2f%%" % (expnames[ex],ds_naos[ex].varexp.isel(mode=0).data.item()*100),)
    
    cb = viz.hcbar(pcm,ax=ax)
    
    figname = "%sAWI-CM3_NAO_Patterns_%s.png" % (figpath,expnames[ex])
    plt.savefig(figname,dpi=150)

#%% Plot EOF2


for ex in range(4):
    fig,ax,_ = viz.init_orthomap(1,1,bboxplot,centlon=-30)
    
    ax       = viz.add_coast_grid(ax,bboxplot)
    
    
    plotvar = eof2[ex]
    # if ex == 3:
    #     plotvar = plotvar * -1
    
    pcm = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,cmap='cmo.balance',vmin=cints[0],vmax=cints[-1])
    ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,colors="k",levels=cints,linewidths=0.55)
    
    ax.set_title("%s EOF2, Var. Explained = %.2f%%" % (expnames[ex],ds_naos[ex].varexp.isel(mode=1).data.item()*100),)
    
    cb = viz.hcbar(pcm,ax=ax)
    
    figname = "%sAWI-CM3_EAP_Patterns_%s.png" % (figpath,expnames[ex])
    plt.savefig(figname,dpi=150)

#%% EOF 3



for ex in range(4):
    fig,ax,_ = viz.init_orthomap(1,1,bboxplot,centlon=-30)
    
    ax       = viz.add_coast_grid(ax,bboxplot)
    
    
    plotvar = eof3[ex]
    # if ex == 3:
    #     plotvar = plotvar * -1
    
    pcm = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,cmap='cmo.balance',vmin=cints[0],vmax=cints[-1])
    ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,colors="k",levels=cints,linewidths=0.55)
    
    ax.set_title("%s EOF3, Var. Explained = %.2f%%" % (expnames[ex],ds_naos[ex].varexp.isel(mode=2).data.item()*100),)
    
    cb = viz.hcbar(pcm,ax=ax)
    
    figname = "%sAWI-CM3_EOF3_Patterns_%s.png" % (figpath,expnames[ex])
    plt.savefig(figname,dpi=150)



plt.show()
#ax.coastlines()

# ============================
#%% Part (2), Do DJF Data Load
# ============================

naopath     = "/home/niu4/gliu8/projects/scrap/nao_crop/"
ncname_eof  = naopath + "%s_nao_index_DJF.nc"
ds_nao_djf     = [xr.open_dataset(ncname_eof % (expnames[ex])).load() for ex in range(4)]

#%% Check the Sign


# Re-Flip Signs
spgbox     = [-60,20,60,80] # Shift NAO box North (for AWI output)
eapbox     = [-60,20,40,60] # Shift Box west for EAP
bbox_check = [spgbox,eapbox,]    
print("Flipping boxes based on [bbox_check]")
nmode_check = len(bbox_check)

for ex in range(4):
    
    eof_exp = ds_nao_djf[ex]
    
    for N in tqdm.tqdm(range(nmode_check)):
        chkbox = bbox_check[N]
        eofreg = proc.sel_region_xr(eof_exp.eofs.isel(mode=N),chkbox)
        regsum = eofreg.sum(('lat','lon'))
        if regsum.data > 0:
            print("Flipping sign for %s, mode %i" % (expnames[ex],N+1))
            
            eof_arr    = eof_exp.eofs.data
            pcs_arr     = eof_exp.pcs.data
            varexp_exp = eof_exp.varexp#.data
            
            eof_arr[N,:,:]  = eof_arr[N,:,:] * -1
            pcs_arr[N,:]     = pcs_arr[N,:] * -1
            
            coords_eof  = dict(mode=eof_exp.mode,lat=eof_exp.lat,lon=eof_exp.lon)
            coords_pc   = dict(mode=eof_exp.mode,time=eof_exp.time)
            #coords_vexp = dict(mode=eof_exp.mode)
            
            eof_new = xr.DataArray(eof_arr,coords=coords_eof,dims=coords_eof,name="eofs")
            pcs_new = xr.DataArray(pcs_arr,coords=coords_pc,dims=coords_pc,name="pcs")
            eof_exp_new = xr.merge([eof_new,pcs_new,varexp_exp])
            
            ds_nao_djf[ex] = eof_exp_new
            


#%% Visualize the Patterns

cints    = np.arange(-750,800,50)
proj     = ccrs.PlateCarree()
bboxplot = [-90,40,20,80]

for ex in range(4): # Experiment loop
    expname = expnames[ex]
    
    for nn in range(3): # Mmode loop
        
        plotvar = ds_nao_djf[ex].isel(mode=nn).eofs
        varexp  = ds_nao_djf[ex].isel(mode=nn).varexp
        
        
        fig,ax,_ = viz.init_orthomap(1,1,bboxplot,centlon=-30)
        ax       = viz.add_coast_grid(ax,bboxplot)
        
        # if ex == 3 and nn == 0:
        #     plotvar = plotvar * -1
        
        pcm = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,cmap='cmo.balance',vmin=cints[0],vmax=cints[-1])
        ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,colors="k",levels=cints,linewidths=0.55)
        
        ax.set_title("%s EOF%i, Var. Explained = %.2f%%" % (expnames[ex],nn+1,varexp*100),)
        
        cb = viz.hcbar(pcm,ax=ax)
        
        figname = "%sAWI-CM3_EOF%i_DJF_Patterns_%s.png" % (figpath,nn+1,expnames[ex])
        plt.savefig(figname,dpi=150)
        

#%% Compute NAO-ENSO Lag, but selecting by month

"""

Lag 0

"""

lags        = 18
selmon_base = np.array([12,1,2])
naobase     = ds_nao_djf[ex].pcs.isel(mode=0)

# For leads
leads = np.arange(lags) * np.nan
llags = leads.copy()
for lag in range(lags):
    
    # Calculate when ENSO leads NAO
    selmon_shift = selmon_base - lag
    print(selmon_shift)
    selmon_shift = [mm%12 for mm in selmon_shift]#selmon_shift % 12
    for mm in range(3): 
        if selmon_shift[mm] == 0:
            selmon_shift[mm] = 12
            
        # if selmon_shift[mm] < 0:
        #     selmon_shift[mm] = 12 + selmon_shift[mm]
    print("Lag %i" % lag)
    print("\t%s" % (selmon_shift))
    
    ensoin = proc.selmon_ds(ensoids[ex],selmon_shift)
    
    rr = np.corrcoef(naobase.data,ensoin.data,)[0,1]
    leads[lag] = rr
    
    # Calculate when ENSO lags NAO
    selmon_shift = selmon_base + lag
    print(selmon_shift)
    selmon_shift = [mm%12 for mm in selmon_shift]#selmon_shift % 12
    for mm in range(3): 
        if selmon_shift[mm] == 0:
            selmon_shift[mm] = 12
            
        # if selmon_shift[mm] < 0:
        #     selmon_shift[mm] = 12 + selmon_shift[mm]
    print("Lead %i" % lag)
    print("\t%s" % (selmon_shift))
    
    ensoin = proc.selmon_ds(ensoids[ex],selmon_shift)
    
    rr = np.corrcoef(naobase.data,ensoin.data,)[0,1]
    llags[lag] = rr
    
    
#%%

fig,ax = plt.subplots(1,1,figsize=(12.5,4))
ax.plot(np.flip(np.arange(lags)),leads)
ax.plot(np.arange(lags),llags)
plt.show()
leads = np.flip(leads)


#%% Look at seasonal correlation

ex           = 0

selmons      = [[12,1,2],[3,4,5],[6,7,8],[9,10,11],[12,1,2],[3,4,5],[6,7,8],[9,10,11],[12,1,2]] #["DJF","MAM","JJA","SON","DJF","MAM","JJA","SON","DJF"]
years        = [    -1  ,      -1,     -1,       -1,      0,      0,      0,        0,      0]

nloop        = len(selmons)


ymin         = naobase.time.dt.year[0].item()
ymax         = naobase.time.dt.year[-1].item()

naobase      = ds_nao_djf[ex].pcs.isel(mode=0)
rr_shiftenso = np.zeros(nloop) * np.nan

for ll in range(nloop):
    
    ensoin      = proc.selmon_ds(ensoids[ex],selmons[ll])
    nyr         = int(ensoin.shape[0]/3)
    ensoin      = ensoin.data.reshape(nyr,3)
    
    naoin       = naobase.data.reshape(nyr,3)
    
    if years[ll] < 0: # ENSO is Leading
        rr          = np.corrcoef(naoin[1:,:].flatten(),ensoin[:(-1),:].flatten())[0,1]
    else:
        rr          = np.corrcoef(naoin[:(-1),:].flatten(),ensoin[1:,:].flatten())[0,1]
        
    rr_shiftenso[ll] = rr.copy()
    
#%%


# selmons_str = [proc.]
# fig,ax   = plt.subplots(1,1,)

# plotlags = np.arange(-4,5,)
# ax.plot(plotlags,rr_shiftenso,)

# ax.set_xlabel("<-- ENSO Leads | --> NAO Leads")


# plt.show()
    
    
    
    
#     naoin  = naobase.sel(time=slice(ymin))
    
##naobase_annual = ds_nao_djf[ex].pcs.isel(mode=0).groupby('time.year').mean('time')




    
    #ensoin = proc.selmon_ds(ensoids[ex],selmon_shift)
    #naoin  = 
    
    

#%% Plot NAO Indices

for ex in range(4):
    
    fig,ax = plt.subplots(1,1,figsize=(12,4.5),constrained_layout=True)
    plotvar = ds_nao_djf[ex].pcs.isel(mode=0).groupby('time.year').mean('time')
    ax.plot(plotvar.year,plotvar,color='k')
    plt.show()
    
    
    

#%% Simple Plots:
    
# Plot Instantaneous
ex          = 0
vmax = 4
for ex in range(5):
    if ex <4:
        ensoin      = proc.selmon_ds(ensoids[ex],[12,1,2])
        naoin       = ds_nao_djf[ex].pcs.isel(mode=0)
        
        fig,ax      = plt.subplots(1,1,constrained_layout=True,figsize=(4.5,4.5))
        
        tcolor      = np.arange(len(ensoin))
        
        sc = ax.scatter(naoin,ensoin,c=tcolor,alpha=0.8,cmap='cmo.deep')
        #ax.scatter(naoin,ensoin,c=expcols[ex],alpha=0.5)
        ax.set_xlabel("NAO (DJF)")
        ax.set_ylabel("Nino3.4 (DJF)")
        ax.set_title("%s \n Instantaneous (r=%.2f)" % (expnames[ex],np.corrcoef(naoin,ensoin)[0,1]))
        
        ax.set_ylim([-vmax,vmax])
    
        ax.set_xlim([-vmax,vmax])
        ax.set_aspect(1)
        
        figname = "%sAWI-CM3_%s_DJFNAO-ENSO-Lag0.png" % (figpath,expnames[ex])
        plt.savefig(figname,dpi=150)
        
        
    else: # Plot Colorbar
        fig,ax      = plt.subplots(1,1,constrained_layout=True,figsize=(4.5,4.5))
        cb = fig.colorbar(sc,ax=ax,fraction=0.025,pad=0.01)
        figname = "%sAWI-CM3_DJFNAO-ENSO-Lag0_Colorbar.png" % (figpath)
        plt.savefig(figname,dpi=150)
        
        


plt.show()
    
#%%


    
# Plot Lag 1 
ex          = 0
for ex in range(4):
    ensoin      = proc.selmon_ds(ensoids[ex],[3,4,5])
    naoin       = ds_nao_djf[ex].pcs.isel(mode=0)
    fig,ax      = plt.subplots(1,1,constrained_layout=True,figsize=(4.5,4.5))
    tcolor      = np.arange(len(ensoin))
    ax.scatter(naoin,ensoin,c=tcolor,alpha=0.8,cmap='cmo.deep')
    #ax.scatter(naoin,ensoin,c=expcols[ex],alpha=0.5)
    ax.set_xlabel("NAO (DJF)")
    ax.set_ylabel("Nino3.4 (MAM)")
    ax.set_title("%s \n Lag 1 (r=%.2f)" % (expnames[ex],np.corrcoef(naoin,ensoin)[0,1]))
    
    ax.set_ylim([-3,3])
    
    ax.set_xlim([-3,3])
    ax.set_aspect(1)
    
    figname = "%sAWI-CM3_%s_DJFNAO-ENSO-Lag1.png" % (figpath,expnames[ex])
    plt.savefig(figname,dpi=150)

#%% Plot Lag 2

for ex in range(4):
    ensoin      = proc.selmon_ds(ensoids[ex],[6,7,8])
    naoin       = ds_nao_djf[ex].pcs.isel(mode=0)
    fig,ax      = plt.subplots(1,1,constrained_layout=True,figsize=(4.5,4.5))
    
    ax.scatter(naoin,ensoin,c=expcols[ex],alpha=0.5)
    ax.set_xlabel("NAO (DJF)")
    ax.set_ylabel("Nino3.4 (JJA)")
    ax.set_title("%s \n Lag 2 (r=%.2f)" % (expnames[ex],np.corrcoef(naoin,ensoin)[0,1]))
    
    ax.set_ylim([-3,3])
    
    ax.set_xlim([-3,3])
    ax.set_aspect(1)
    
    figname = "%sAWI-CM3_%s_DJFNAO-ENSO-Lag2.png" % (figpath,expnames[ex])
    plt.savefig(figname,dpi=150)



#%%  Lag -1
for ex in range(4):
    ensoin      = proc.selmon_ds(ensoids[ex],[9,10,11])
    naoin       = ds_nao_djf[ex].pcs.isel(mode=0)
    fig,ax      = plt.subplots(1,1,constrained_layout=True,figsize=(4.5,4.5))
    
    ax.scatter(naoin,ensoin,c=expcols[ex],alpha=0.5)
    ax.set_xlabel("NAO (DJF)")
    ax.set_ylabel("Nino3.4 (SON)")
    ax.set_title("%s \n Lag -1 (r=%.2f)" % (expnames[ex],np.corrcoef(naoin,ensoin)[0,1]))
    
    ax.set_ylim([-3,3])
    
    ax.set_xlim([-3,3])
    ax.set_aspect(1)
    
    figname = "%sAWI-CM3_%s_DJFNAO-ENSO-Lag-1.png" % (figpath,expnames[ex])
    plt.savefig(figname,dpi=150)
    

#%% Attempting Again, but Applying 1-year Lag

vmax = 4
for ex in range(4):
    
    # Apply 3 month lead/lag
    ensoin      = proc.selmon_ds(ensoids[ex],[12,1,2])[:(-3)]
    naoin       = ds_nao_djf[ex].pcs.isel(mode=0)[3:]
    
    tcolor      = np.arange(len(ensoin))
    
    # Apply 1-Year Lag
    #ystart     = ensoin.time[0].dt.year
    #yend       = ensoin.time[1].dt.year
    
    
    fig,ax      = plt.subplots(1,1,constrained_layout=True,figsize=(4.5,4.5))
    
    ax.scatter(naoin,ensoin,c=tcolor,alpha=0.8,cmap='cmo.deep')
    ax.set_xlabel("NAO (DJF,Year 1)")
    ax.set_ylabel("Nino3.4 (DJF, Year 0)")
    ax.set_title("%s \n 1-Year Lag (r=%.2f)" % (expnames[ex],np.corrcoef(naoin,ensoin)[0,1]))
    
    ax.set_ylim([-vmax,vmax])
    
    ax.set_xlim([-vmax,vmax])
    ax.set_aspect(1)
    
    figname = "%sAWI-CM3_%s_DJFNAO-ENSO-Lag1Year.png" % (figpath,expnames[ex])
    plt.savefig(figname,dpi=150)

plt.show()

#%% Next Part (3) Descriptive Analysis


eof1s       = [ds.eofs.isel(mode=0) for ds in ds_naos]
pc1s        = [ds.pcs.isel(mode=0) for ds in ds_naos]

nao_monstd  = [ds.groupby('time.month').std('time') for ds in pc1s]
enso_monstd = [ds.groupby('time.month').std('time') for ds in ensoids]


#%% Remake into Pandas Dataframe for plotting

import pandas as pd
import seaborn as sns


dfs = []
for ex in range(4):
    dfin = enso_monstd[ex].to_pandas()
    dfin = pd.DataFrame(
        {'month':dfin.index,
         'sd':dfin.values,
         'exp':expnames[ex],
         }
        )
    dfs.append(dfin)


#[group,col,val]
#[month,simulation,value]

#%% Plot Grouped Bar Plot (31km)

dfmerge = pd.concat([dfs[0], dfs[1]], ignore_index=True, sort=False) #pd.merge(dfs[0],dfs[1])

fig,ax  = plt.subplots(1,1,constrained_layout=True)
sns.barplot(data=dfmerge, x='month', y='sd', hue='exp',palette=[[.15,.15,.15],'tomato'],ax=ax)
ax.set_ylim([0,2.5])
#dfmerge.pivot(("month", "exp", "sd")).plot(kind='bar',ax=ax)
plt.show()

#%% Plot Grouped Bar Plot (9km)

dfmerge = pd.concat([dfs[2], dfs[3]], ignore_index=True, sort=False) #pd.merge(dfs[0],dfs[1])

fig,ax = plt.subplots(1,1,constrained_layout=True)
sns.barplot(data=dfmerge, x='month', y='sd', hue='exp',palette=[[.15,.15,.15],'tomato'],ax=ax)
ax.set_ylim([0,2.5])
#dfmerge.pivot(("month", "exp", "sd")).plot(kind='bar',ax=ax)
plt.show()

#%% Repeat for NAO


dfs = []
for ex in range(4):
    dfin = nao_monstd[ex].to_pandas()
    dfin = pd.DataFrame(
        {'month':dfin.index,
         'sd':dfin.values,
         'exp':expnames[ex],
         }
        )
    dfs.append(dfin)


dfmerge = pd.concat([dfs[0], dfs[1]], ignore_index=True, sort=False) #pd.merge(dfs[0],dfs[1])

fig,ax = plt.subplots(1,1,constrained_layout=True)
sns.barplot(data=dfmerge, x='month', y='sd', hue='exp',palette=[[.15,.15,.15],'tomato'],ax=ax)
ax.set_ylim([0,2.5])
#dfmerge.pivot(("month", "exp", "sd")).plot(kind='bar',ax=ax)
plt.show()



dfmerge = pd.concat([dfs[2], dfs[3]], ignore_index=True, sort=False) #pd.merge(dfs[0],dfs[1])

fig,ax = plt.subplots(1,1,constrained_layout=True)
sns.barplot(data=dfmerge, x='month', y='sd', hue='exp',palette=[[.15,.15,.15],'tomato'],ax=ax)
ax.set_ylim([0,2.5])
#dfmerge.pivot(("month", "exp", "sd")).plot(kind='bar',ax=ax)
plt.show()

#%% Calculate Sliding variance and Correlation

winlen  = 240
ex      = 0

calc_r1 = lambda ts: np.corrcoef(ts[:(-1)],ts[1:])[0,1]

ds_sliding = []
for ex in range(2):
    
    ensoin = ensoids[ex]
    naoin  = pc1s[ex]
    
    ntime     = len(naoin)
    
    starttime = []
    ensovar   = []
    naovar    = []
    nccorr    = []
    ensor1    = []
    naor1     = []
    
    for i in tqdm.tqdm(range(ntime-winlen)):
        
        trange  = np.arange(i,i+winlen)
        
        naosub  = naoin[trange]
        ensosub = ensoin[trange]
        
        starttime.append(naosub.time[0].data)
        ensovar.append(ensosub.std('time'))
        naovar.append(naosub.std('time'))
        nccorr.append(np.corrcoef(naosub,ensosub)[0,1])
        ensor1.append(calc_r1(ensosub))
        naor1.append(calc_r1(naosub))
    
    coords      = dict(time=starttime)
    
    
    da_ensovar  = xr.DataArray(ensovar,coords=coords,dims=coords,name='enso_sd')
    da_naovar   = xr.DataArray(naovar,coords=coords,dims=coords,name='nao_sd')
    da_corr     = xr.DataArray(nccorr,coords=coords,dims=coords,name='corr')
    da_ensor1   = xr.DataArray(ensor1,coords=coords,dims=coords,name='enso_r1')
    da_naor1    = xr.DataArray(naor1,coords=coords,dims=coords,name='nao_r1')
    
    da_out      = xr.merge([da_ensovar,da_naovar,da_corr,da_ensor1,da_naor1])
    
    ds_sliding.append(da_out)





#%% Make Some Plots of the Changes
ex      = 1
lw      = 2.5


if ex == 1:
    ex2 = 3
elif ex == 0:
    ex2 = 2
    
cols2 = ["","",[.15,.15,.15],'tomato']

fig,axs = plt.subplots(3,1,constrained_layout=True,figsize=(12.5,6))

ax      = axs[0]
plotvar = ds_sliding[ex].enso_sd
ax.plot(plotvar.time,plotvar,label="ENSO SD (31km)",lw=lw,color="hotpink")
for ex2 in [2,3]:
    ex2std = ensoids[ex2].std('time')
    ax.axhline(ex2std,ls='dotted',color=cols2[ex2],lw=2,label=expnames[ex2])

ax = axs[1]
plotvar = ds_sliding[ex].nao_sd
ax.plot(plotvar.time,plotvar,label="NAO SD (31km)",lw=lw,color=[.15,.15,.15])
for ex2 in [2,3]:
    ex2std = pc1s[ex2].std('time')
    print(ex2std)
    ax.axhline([ex2std],ls='dotted',color=cols2[ex2],lw=2,label=expnames[ex2])

ax = axs[2]
plotvar = ds_sliding[ex].corr
ax.plot(plotvar.time,plotvar,label="ENSO-NAO Correlation (31km)",lw=lw,color='violet')
for ex2 in [2,3]:
    ex2std = np.corrcoef(pc1s[ex2],ensoids[ex2])[0,1]#pc1s[ex2].std('time')
    if ex2==3:
        ex2std*=-1 #no sign corr to nao
    ax.axhline([ex2std],ls='dotted',color=cols2[ex2],lw=2,label=expnames[ex2])

for ax in axs:
    ax.legend()
    ax.set_xlim([plotvar.time[0],plotvar.time[-1]])
plt.show()

#%%
ex      = 1
lw      = 2.5


if ex == 1:
    ex2 = 3
elif ex == 0:
    ex2 = 2
    
cols2 = ["","",[.15,.15,.15],'tomato']

fig,axs = plt.subplots(5,1,constrained_layout=True,figsize=(12.5,8))

ax      = axs[0]
plotvar = ds_sliding[ex].enso_sd
ax.plot(plotvar.time,plotvar,label="ENSO SD (31km)",lw=lw,color="hotpink")
for ex2 in [2,3]:
    ex2std = ensoids[ex2].std('time')
    ax.axhline(ex2std,ls='dotted',color=cols2[ex2],lw=2,label=expnames[ex2])

ax = axs[1]
plotvar = ds_sliding[ex].nao_sd
ax.plot(plotvar.time,plotvar,label="NAO SD (31km)",lw=lw,color=[.15,.15,.15])
for ex2 in [2,3]:
    ex2std = pc1s[ex2].std('time')
    print(ex2std)
    ax.axhline([ex2std],ls='dotted',color=cols2[ex2],lw=2,label=expnames[ex2])

ax = axs[2]
plotvar = ds_sliding[ex].corr
ax.plot(plotvar.time,plotvar,label="ENSO-NAO Correlation (31km)",lw=lw,color='violet')
for ex2 in [2,3]:
    ex2std = np.corrcoef(pc1s[ex2],ensoids[ex2])[0,1]#pc1s[ex2].std('time')
    if ex2==3:
        ex2std*=-1 #no sign corr to nao
    ax.axhline([ex2std],ls='dotted',color=cols2[ex2],lw=2,label=expnames[ex2])

ax = axs[3]
plotvar = ds_sliding[ex].enso_r1
ax.plot(plotvar.time,plotvar,label="",lw=lw,color='violet')
for ex2 in [2,3]:
    ex2std = calc_r1(ensoids[ex2])#np.corrcoef(pc1s[ex2],ensoids[ex2])[0,1]#pc1s[ex2].std('time')
    if ex2==3:
        ex2std*=-1 #no sign corr to nao
    ax.axhline([ex2std],ls='dotted',color=cols2[ex2],lw=2,label=expnames[ex2])
    

ax = axs[3]
plotvar = ds_sliding[ex].nao_r1
ax.plot(plotvar.time,plotvar,label="",lw=lw,color='violet')
for ex2 in [2,3]:
    ex2std = calc_r1(s[ex2])#pc1s[ex2].std('time')
    if ex2==3:
        ex2std*=-1 #no sign corr to nao
    ax.axhline([ex2std],ls='dotted',color=cols2[ex2],lw=2,label=expnames[ex2])
    



for ax in axs:
    ax.legend()
    ax.set_xlim([plotvar.time[0],plotvar.time[-1]])
plt.show()


    
#%%


if ex == 1:
    ex2 = 3
elif ex == 0:
    ex2 = 2

cols2 = ["","",[.15,.15,.15],'tomato']

fig,axs = plt.subplots(5,1,constrained_layout=True,figsize=(12.5,8))

ax      = axs[0]
plotvar = ds_sliding[ex].enso_sd
ax.plot(plotvar.time,plotvar,label="ENSO SD (31km)",lw=lw,color="hotpink")
for ex2 in [2,3]:
    ex2std = ensoids[ex2].std('time')
    ax.axhline(ex2std,ls='dotted',color=cols2[ex2],lw=2,label=expnames[ex2])
# # also plot range of natural variabilie
# mu       = ds_sliding[0].enso_sd.mean('time').data
# sigma    = 2*ds_sliding[0].enso_sd.std('time').data
# tplot    = plotvar.time
# ntime    = len(tplot)
# ax.fill_between(plotvar.time,np.ones(ntime) * mu - sigma, np.ones(ntime)*mu + sigma,color='gray',alpha=0.25)




ax = axs[1]
plotvar = ds_sliding[ex].nao_sd
ax.plot(plotvar.time,plotvar,label="NAO SD (31km)",lw=lw,color=[.15,.15,.15])
for ex2 in [2,3]:
    ex2std = pc1s[ex2].std('time')
    print(ex2std)
    ax.axhline([ex2std],ls='dotted',color=cols2[ex2],lw=2,label=expnames[ex2])
# # also plot range of natural variabilie
# mu       = ds_sliding[0].nao_sd.mean('time').data
# sigma    = 2*ds_sliding[0].nao_sd.std('time').data
# tplot    = plotvar.time
# ntime    = len(tplot)
# ax.fill_between(plotvar.time,np.ones(ntime) * mu - sigma, np.ones(ntime)*mu + sigma,color='gray',alpha=0.25)


ax = axs[2]
plotvar = ds_sliding[ex].corr
ax.plot(plotvar.time,plotvar,label="ENSO-NAO Correlation (31km)",lw=lw,color='violet')
for ex2 in [2,3]:
    ex2std = np.corrcoef(pc1s[ex2],ensoids[ex2])[0,1]#pc1s[ex2].std('time')
    if ex2==3:
        ex2std*=-1 #no sign corr to nao
    ax.axhline([ex2std],ls='dotted',color=cols2[ex2],lw=2,label=expnames[ex2])
# # also plot range of natural variabilie
# mu       = ds_sliding[0].corr.mean('time').data
# sigma    = 2*ds_sliding[0].corr.std('time').data
# tplot    = plotvar.time
# ntime    = len(tplot)
# ax.fill_between(plotvar.time,np.ones(ntime) * mu - sigma, np.ones(ntime)*mu + sigma,color='gray',alpha=0.25)




ax = axs[3]
plotvar = ds_sliding[ex].enso_r1
ax.plot(plotvar.time,plotvar,label="ENSO R1",lw=lw,color='hotpink')
for ex2 in [2,3]:
    ex2std = calc_r1(ensoids[ex2])#np.corrcoef(pc1s[ex2],ensoids[ex2])[0,1]#pc1s[ex2].std('time')

    ax.axhline([ex2std],ls='dotted',color=cols2[ex2],lw=2,label=expnames[ex2])
# # also plot range of natural variabilie
# mu       = ds_sliding[0].enso_r1.mean('time').data
# sigma    = 2*ds_sliding[0].enso_r1.std('time').data
# tplot    = plotvar.time
# ntime    = len(tplot)
# ax.fill_between(plotvar.time,np.ones(ntime) * mu - sigma, np.ones(ntime)*mu + sigma,color='gray',alpha=0.25)



ax = axs[4]
plotvar = ds_sliding[ex].nao_r1
ax.plot(plotvar.time,plotvar,label="NAO R1",lw=lw,color=[.15,.15,.15])
for ex2 in [2,3]:
    ex2std = calc_r1(pc1s[ex2])#pc1s[ex2].std('time')

    ax.axhline([ex2std],ls='dotted',color=cols2[ex2],lw=2,label=expnames[ex2])
# # also plot range of natural variabilie
# mu       = ds_sliding[0].nao_r1.mean('time').data
# sigma    = 2*ds_sliding[0].nao_r1.std('time').data
# tplot    = plotvar.time
# ntime    = len(tplot)
# ax.fill_between(plotvar.time,np.ones(ntime) * mu - sigma, np.ones(ntime)*mu + sigma,color='gray',alpha=0.25)






for ax in axs:
    ax.legend(ncol=3)
    ax.set_xlim([plotvar.time[0],plotvar.time[-1]])
plt.show()






#%%



#%% Make Some Exploratory Plots