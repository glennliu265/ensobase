#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Composite ENSO events by threshold

- Copied from enso_event_identifier

Created on Wed Oct 15 16:37:08 2025

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

#%% Indicate paths

figpath         = "/home/niu4/gliu8/figures/bydate/2025-10-21/"
proc.makedir(figpath)

datpath         = "/home/niu4/gliu8/projects/scrap/TP_crop/"

expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C","glorys"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950","GLORYS"]



vname           = "sst"#"Dmaxgrad" #"sst"#"str"

vnames          = ['sst','ssr','str','Dmaxgrad',"D20"]
nvars           = len(vnames)

ninoid_name     = "nino34"


ninopath = "/home/niu4/gliu8/projects/scrap/nino34/"

nexps    = len(expnames)

#%% Some Helper Functions  (Copied from calculate enso response)

def preprocess_enso(ds):
    # Remove Mean Seasonal Cycle and the Quadratic Trend
    dsds   = proc.xrdeseason(ds)
    
    # Compute Spectra
    def detrend_quadratic(ds):
        x = np.arange(len(ds))
        y = ds.data
        if np.any(np.isnan(y)):
            return np.ones(y.shape)*np.nan
        ydetrended,model=proc.detrend_poly(x,y,2)
        return ydetrended
        
    st = time.time()
    dsanom = xr.apply_ufunc(
        detrend_quadratic,  # Pass the function
        dsds,  # The inputs in order that is expected
        # Which dimensions to operate over for each argument...
        input_core_dims=[['time'],],
        output_core_dims=[['time'],],  # Output Dimension
        vectorize=True,  # True to loop over non-core dims
        )
    print("Detrended in %.2fs" % (time.time()-st))
    #dsanom = proc.xrdetrend_1d(ds,2) 
    return dsanom

def swap_rename(ds,chkvar,newvar):
    if chkvar in list(ds.coords):
        print("Renaming [%s] to [%s]" % (chkvar,newvar))
        ds = ds.rename({chkvar:newvar})
    return ds

def standardize_names(ds):
    
    ds = swap_rename(ds,'time_counter','time')
    ds = swap_rename(ds,"TIME_COUNTER",'time')
    ds = swap_rename(ds,"LON","lon")
    ds = swap_rename(ds,"LAT","lat")
    ds = swap_rename (ds,"LAT232_409","lat")
    return ds

def init_tp_map():
    bbplot = [120, 290, -20, 20]
    fix_lon = np.hstack([np.arange(120,190,10),np.arange(-180,-60,10)])
    proj   = ccrs.PlateCarree(central_longitude=180)
    projd  = ccrs.PlateCarree()
    fig,ax = plt.subplots(1,1,figsize=(12.5,4.5),subplot_kw={'projection':proj})
    ax.set_extent(bbplot)
    ax     = viz.add_coast_grid(ax,bbox=bbplot,fill_color='k',
                                proj=ccrs.PlateCarree(),fix_lon=fix_lon,ignore_error=True)
    return fig,ax



#%% Load ENSO Indices


ds_enso = []

ensoids = []
for ex in range(nexps):
    ds = ut.load_ensoid(expnames[ex],ninoid_name,standardize=False)
    ds_enso.append(ds)
    ds = ds * ds['std'].data.item()
    ensoids.append(ds)


#%% Identify Events

tol       = 12
ninodicts = []
ninadicts = []
for ex in tqdm.tqdm(range(nexps)):
    
    ensoin  = ensoids[ex]
    sigma   = ensoin.std('time')
    id_nino = ensoin >= sigma
    id_nina = ensoin <= -sigma
    
    oout = ut.combine_events(ensoin,id_nino,tol=tol,verbose=False)
    aout = ut.combine_events(ensoin,id_nina,tol=tol,verbose=False)
    ninodicts.append(oout)
    ninadicts.append(aout)
    print("For %s:" % expnames_long[ex])
    print("\tIdentified %i Nino Events" % (len(oout['event_combine'])))
    print("\tIdentified %i Nina Events" % (len(aout['event_combine'])))

#%% Get Stacked Events

stacked_events = []
for ninotype in ["nino","nina"]:
        
    if ninotype == "nino":
        indict = ninodicts
    elif ninotype == "nina":
        indict = ninadicts
    
    exstack = []
    for ex in range(nexps):
        expdict     = indict[ex]
        print(ex)
        ibefore     = 12
        iafter      = 36
        plotlags    = np.hstack([np.flip((np.arange(0,ibefore+1) * -1)),np.arange(1,iafter+1,1)])
        target_var  = ensoids[ex].data
        eventids    = expdict['center_ids']
        
        exstack.append(ut.stack_events(target_var,eventids,ibefore,iafter))
    stacked_events.append(exstack)
    
#%% Plot Stacked Events and select threshold, subset events


thres       = 1.0
ex          = 0
ninotype    = "nino"
xtkslag     = np.arange(-12,39,3)

# Select the Experiment
expdict     = indict[ex]
if ninotype == "nino":
    indict = ninodicts
    ninon  = 0
elif ninotype == "nina":
    indict = ninadicts
    ninon  = 1
eventids           = expdict['center_ids']
nevents            = len(eventids)
stacked_events_in  = stacked_events[ninon][ex]
eamplitudes        = expdict['event_max']

fig,ax      = plt.subplots(1,1,constrained_layout=True,figsize=(8,4.5))
#ax.set_prop_cycle(color=list(mpl.colormaps['Accent'].colors))

eventsel_comp = []
eventsel_id   = []
sigma         = ensoids[ex]['std'].item()
for ie in range(nevents):
    
    if ie == 0:
        lab = "Individual Event"
    else:
        lab = ""
    
    plotvar = stacked_events_in[ie,:]
    if eamplitudes[ie] >= thres:
        eventsel_comp.append(plotvar)
        eventsel_id.append(eventids[ie])
        plotc = 'red'
    else:
        plotc = 'gray'
    ax.plot(plotlags,plotvar,lw=.80,alpha=0.55,label=lab,c=plotc)
    #ax.plot([0,],eamplitudes[ie],marker="o",color=plotc)
    
# Plot Events Above Threshold    
ax.axhline([thres],ls='solid',color="firebrick",lw=1.5,label="Threshold=%.2f $\degree C$" % thres)    
ax.plot(plotlags,np.array(eventsel_comp).mean(0),c='firebrick',label="Selected Events (n=%i)" % len(eventsel_id))

ax.plot(plotlags,np.nanmean(stacked_events_in,0),color="k",lw=1,alpha=1,label='Mean')

ax.legend()

ax.axvline([0],ls='solid',color="k",lw=0.50)
ax.axhline([0],ls='solid',color="k",lw=0.50)
ax.axhline([sigma],ls='dotted',color="red",lw=0.75)
ax.axhline([-sigma],ls='dotted',color="blue",lw=0.75)

ax.set_xlabel("Lags (Months)")
ax.set_ylabel("SST Anomaly $\degree C$")
if expnames[ex] == "TCo319_ssp585":
    ax.set_ylim([-7.5,7.5])
else:
    ax.set_ylim([-4.5,4.5])
ax.set_xticks(xtkslag)
ax.set_xlim([plotlags[0],plotlags[-1]])

ax.set_title("%s Composites (%s)\n$\sigma=$%.2f" %(ninotype,expnames_long[ex],sigma))

figname = "%s%s_%s_Spaghetti_thres%.2f.png" % (figpath,expnames[ex],ninoid_name,thres)
plt.savefig(figname,dpi=150,bbox_inches='tight')
plt.show()

#%% Make the composites

#%% Load variable of choice

vname           = "sst"
datpath_anom    = datpath + "anom_detrend2/"
# Load Dataset
ncname = "%s%s_%s_anom.nc" % (datpath_anom,expnames[ex],vname)

try:
    ds = xr.open_dataset(ncname).load()[vname]
except:

    print("Could not find %s %s" % (vname,expnames[ex]))

#%% Get event ids

dssubset,timessubset = ut.stack_events_2d(ds,eventsel_id,12,18,times_da=ds.time)

#%% Plot selected lags (quick plot)

sellags      = [0,3,6,9,12]
compositesel = dssubset.mean('eventid')

for ii in range(len(sellags)):
    
    compositesel.sel(lag=sellags[ii]).plot(vmin=-1.5,vmax=1.5,cmap='cmo.balance')
    plt.show()
    
#%% Make Nicer plots (in the style of the regression plots from before)

sellags      = [0,3,6,9,12]
compositesel = dssubset.mean('eventid')

sstdict = dict(
    vname      = "sst",
    longname   = "sea surface temperature",
    units      = "$\degree C$",
    cints_mean = np.arange(18,35.5,.5),
    cmap_mean  = "cmo.balance",
    cints_enso = np.arange(-2.5,2.6,0.1),
    )

indict = sstdict
vname  = 'sst'

cints   = np.arange(-1.5,1.6,0.1)#indict['cints_enso']#np.arange(-1,1.1,0.1)#np.arange(-5,5.5,0.5)#np.arange(-24,26,2)
cints_over = np.arange(1.6,4.4,0.4)
vmax    = cints[0]  
cmap    = indict['cmap_mean']
projd   = ccrs.PlateCarree()
#indict  = vdicts[vname]


for ii in range(len(sellags)):
    
    fig,ax  = ut.init_tp_map()
    plotvar = compositesel.sel(lag=sellags[ii])
    itlag    = plotvar.lag.data.item()
    
    # Plot the variable
    pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,
                            transform=projd,cmap=cmap,vmin=-vmax,vmax=vmax)
    
    cl      = ax.contour(plotvar.lon,plotvar.lat,plotvar,
                         linewidths=0.55,levels=cints,
                            transform=projd,colors='k')
    
    clo      = ax.contour(plotvar.lon,plotvar.lat,plotvar,
                         linewidths=0.55,levels=cints_over,
                            transform=projd,colors='w')
    
    
    ax.clabel(cl)
    ax.clabel(clo)
    
    cb      = fig.colorbar(pcm,ax=ax,fraction=0.01,pad=0.01)
    cb.set_label("%s [%s per $1\sigma$ %s]" % (vname,indict['units'],ninoid_name))
    ax.set_title("AWI-CM3 (%s) %s Regression, Lag %02i" % (expnames_long[ex],indict['longname'],itlag))
    
    figname = "%sENSO_Composite_%s_%s_thres%.2f_lag%02i.png" % (figpath,vname,expnames[ex],thres,itlag)
    plt.savefig(figname,dpi=150,bbox_inches='tight')

    plt.close()

