#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Identify ENSO events within a timeseries

- Copied Upper Section of calculate_enso_response

Created on Mon Sep 29 15:30:15 2025

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

figpath         = "/home/niu4/gliu8/figures/bydate/2025-10-14/"
proc.makedir(figpath)

datpath         = "/home/niu4/gliu8/projects/scrap/TP_crop/"

expnames        = ["TCo2559-DART-1950C","TCo319_ctl1950d","TCo319_ssp585",]#"TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C","glorys"]
expnames_long   = ["5km 1950","31km Control","31km SSP585",]#,"9km 1950","9km 2090","5km 1950","GLORYS"]

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
    
    # Old version (load without function) ---------
    # ninonc      = "%s%s_nino34.nc" % (ninopath,expnames[ex])
    # ds          = xr.open_dataset(ninonc).load()
    
    # ds_enso.append(ds)
    # ensoids.append(ds.sst * ds['std'].data.item())
    # ------------------
    
# =========
#%% Identify Events (Scrap)

#ensoin = ensoids[-1]


"""

1): Identify Events (+/- 1stdev)
    
2): Combine like events (take maximum as "center")

3): Make into DataArray

"""

#%%
# ensoin  = ensoids[-1]
# sigma   = ensoin.std('time')
# id_nino = ensoin >= sigma
# id_nina = ensoin <= -sigma



# print("Identified %i events..." % (nevents))




# # num2next     = np.roll(eventid,-1) - eventid # Time to next event
# # num2next[-1] = (eventid[-1] - eventid[-2]).item()
# # consecutive     = (num2next == 1) # Identify Consecutive Events
# # separate_events = eventid[~consecutive]


# var_in          = ensoin#[-1]
# id_in           = id_nino

# # Separate into discrete events
# nevents        = id_in.data.sum().item()
# eventids       = np.where(id_in)[0]
# event_combine = []
# for ii in range(nevents+1): #event_id:
#     if ii == (nevents):
#         print("Merging last event: %s" % event_merge)
#         event_combine.append(event_merge)
#         continue
#     ievent = eventids[ii].item()#.data
    
#     if ii == 0:
#         prev_id = 0
#         event_merge = [ievent]
#         continue
    
#     if (ievent - prev_id) == 1: # Consecutive Event
#         event_merge.append(ievent)
#         print("%i is consecutive to previous events (%s)" % (ievent,event_merge))
#     else: # Otherwise, just add event and merge
#         event_combine.append(event_merge)
#         event_merge = [ievent,] # Make a new one
#         print("Making new event sequence at %i" % (ievent))
#     prev_id = ievent
# print("Identified %i events!" % len(event_combine))

# # Identify Event Centers

# ncevent         = len(event_combine)
# event_time      = []
# center_ids      = []
# for ie in range(ncevent):
#     idsin = event_combine[ie]
#     if len(idsin) == 1: # Only 1 step
#         event_time.append(var_in.time.isel(time=idsin[0]))
#         center_ids.append(idsin[0])
#     else:
#         amplitudes = var_in.isel(time=idsin) #.argmax()
#         idmax      = np.argmax(np.abs(amplitudes.data)).item()
#         event_time.append(var_in.time.isel(time=idsin[idmax]))
#         center_ids.append(idsin[idmax])

#%% Make above into function

def combine_events(var_in,id_in,tol=1,return_dict=False,verbose=True):
    
    # Separate into discrete events
    nevents        = id_in.data.sum().item()
    eventids       = np.where(id_in)[0]
    event_combine = []
    for ii in range(nevents+1): #event_id:
        if ii == (nevents):
            if verbose:
                print("Merging last event: %s" % event_merge)
            event_combine.append(event_merge)
            continue
        
        ievent = eventids[ii].item()#.data
        
        if ii == 0:
            prev_id     = ievent
            event_merge = [ievent,]
            continue
        
        if (ievent - prev_id) <= tol: # Consecutive Event
            event_merge.append(ievent)
            if verbose:
                print("%i is consecutive to previous events (%s)" % (ievent,event_merge))
        else: # Otherwise, just add event and merge
            event_combine.append(event_merge)
            event_merge = [ievent,] # Make a new one
            if verbose:
                print("Making new event sequence at %i" % (ievent))
        prev_id = ievent
    if verbose:
        print("Identified %i events!" % len(event_combine))

    # Identify Event Centers
    ncevent         = len(event_combine)
    event_time      = []
    center_ids      = []
    event_max       = []
    for ie in range(ncevent):
        idsin = event_combine[ie]
        if len(idsin) == 1: # Only 1 step
            event_time.append(var_in.time.isel(time=idsin[0]))
            center_ids.append(idsin[0])
        else:
            amplitudes = var_in.isel(time=idsin) #.argmax()
            idmax      = np.argmax(np.abs(amplitudes.data)).item()
            event_max.append(np.nanmax(np.abs(amplitudes.data)).item())
            event_time.append(var_in.time.isel(time=idsin[idmax]))
            center_ids.append(idsin[idmax])
    
    if return_dict:
        
        def get_mon(ds):
            return [da.time.dt.month.data.item() for da in ds]
        
        def get_duration(ds):
            return [len(ts) for ts in ds]
        
        durations = get_duration(event_combine)
        months    = get_mon(event_time)
        
        outdict = dict(event_time=event_time,
                       center_ids=center_ids,
                       event_combine=event_combine,
                       event_max=event_max,
                       durations=durations,
                       eventmonths=months)
        return outdict
    return event_time,center_ids,event_combine

def stack_events(target_var,eventids,ibefore,iafter):
    # Stack events between -ibefore and +iafter months
    
    nevents        = len(eventids)
    plotlags       = np.hstack([np.flip((np.arange(0,ibefore+1) * -1)),np.arange(1,iafter+1,1)])
    stacked_events = np.zeros((nevents,len(plotlags))) * np.nan
    ntime          = len(target_var)
    for ie in range(nevents):
        
        ievent    = eventids[ie]
        istart    = ievent-ibefore
        iend      = ievent+iafter
        
        if (istart >=0) and (iend < ntime):
            stacked_events[ie,:] = target_var[istart:(iend+1)]
            
        elif iend >= ntime:
            filler = np.zeros( (iend-ntime+1)) * np.nan
            subset = np.hstack([target_var[istart:],filler])
            stacked_events[ie,:] = subset
        elif istart < 0: # Note havent tested this
            filler  = np.zeros(np.abs(istart)) * np.nan
            subset  = np.hstack([filler,target_var[:(iend+1)],])
            stacked_events[ie,:] = subset
    return stacked_events
            


# ninotimes,center_ids_nino,nino_combine = combine_events(ensoin,id_nino)
# ninatimes,center_ids_nina,nina_combine = combine_events(ensoin,id_nina)

# #%% Plot Identified Events


    

# fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(12.5,4.5))
# ax.plot(ensoin.time,ensoin)

# sigma  = ensoin.std('time')

# ax.axhline([0],ls='solid',lw=0.55,c='k')
# ax.axhline([sigma],ls='dashed',lw=0.55,c='r')
# ax.axhline([-sigma],ls='dashed',lw=0.55,c='b')


# ninos = ensoin.isel(time=center_ids_nino)
# ax.plot(ninos.time,ninos,c='r',marker="o",markersize=20,linestyle='none')

# ninas = ensoin.isel(time=center_ids_nina)
# ax.plot(ninas.time,ninas,c='b',marker="d",markersize=20,linestyle='none')


# plt.show()

# =========================================
#%% Now repeat identification for each case
# =========================================

tol       = 12
ninodicts = []
ninadicts = []
for ex in tqdm.tqdm(range(nexps)):
    
    ensoin  = ensoids[ex]
    sigma   = ensoin.std('time')
    id_nino = ensoin >= sigma
    id_nina = ensoin <= -sigma
    
    oout = combine_events(ensoin,id_nino,tol=tol,return_dict=True,verbose=False)
    aout = combine_events(ensoin,id_nina,tol=tol,return_dict=True,verbose=False)
    ninodicts.append(oout)
    ninadicts.append(aout)
    print("For %s:" % expnames_long[ex])
    print("\tIdentified %i Nino Events" % (len(oout['event_combine'])))
    print("\tIdentified %i Nina Events" % (len(aout['event_combine'])))
    
#%% Identify Frequency by Month + Duration

def get_mon(ds):
    return [da.time.dt.month.data.item() for da in ds]

def get_duration(ds):
    return [len(ts) for ts in ds]

ninomons = [get_mon(ds['event_time']) for ds in ninodicts]
ninamons = [get_mon(ds['event_time']) for ds in ninadicts]

# Also get # of years for reference 
expyrs = [int(len(ts)/12) for ts in ensoids]

# Count length of the event
ninodurations = [get_duration(ds['event_combine']) for ds in ninodicts]
ninadurations = [get_duration(ds['event_combine']) for ds in ninadicts]

# Print it
print("Average Duration for El Nino Events:")
for ex in range(nexps):
    print("\t%s    \t\t: \t%i months" % (expnames[ex],np.nanmean(ninodurations[ex])))
print("Average Duration for La Nina Events:")
for ex in range(nexps):
    print("\t%s    \t\t: \t%i months" % (expnames[ex],np.nanmean(ninadurations[ex])))



# ==================================================
#%% Temprarily commenting below to sppeed up visualize ation of events
    
# #%% Plot Identified ENSO Events

# for ex in range(nexps):

#     ensoin = ensoids[ex]
    
#     fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(12.5,4.5))
    
#     ax.plot(ensoin.time,ensoin,c='k',lw=2.5)
#     sigma  = ensoin.std('time')
    
#     ax.axhline([0],ls='solid',lw=0.55,c='k')
#     ax.axhline([sigma],ls='dashed',lw=0.55,c='r')
#     ax.axhline([-sigma],ls='dashed',lw=0.55,c='b')
    
#     ninos = ensoin.isel(time=ninodicts[ex]['center_ids'])
#     ax.plot(ninos.time,ninos,c='r',marker="o",markersize=10,linestyle='none')
    
#     ninas = ensoin.isel(time=ninadicts[ex]['center_ids'])
#     ax.plot(ninas.time,ninas,c='b',marker="d",markersize=10,linestyle='none')
    
    
#     ax.set_xlabel("Time (Years)")
#     ax.set_ylabel(r"Ni$\tilde{n}$o3.4 Index [$\degree$C]")
#     title = "AWI-CM3 (%s)\n 1$\sigma$=%.4f, " % (expnames_long[ex],sigma) + r"#$Ni\tilde{n}o$: %i, #$Ni\tilde{n}a$: %i" % (len(ninos),len(ninas))
#     ax.set_title(title)
    
#     if expnames[ex] == "TCo319_ssp585":
#         ax.set_ylim([-7.5,7.5])
#     else:
#         ax.set_ylim([-4.5,4.5])
#     ax.grid(True,which='both',ls='dotted',lw=0.50)
    
#     ax.set_xlim([ensoin.time.isel(time=0),ensoin.time.isel(time=-1)])
    
#     # Save the figure
#     figname = "%sNino_Events_Identified_%s_tol%02imon.png" % (figpath,expnames[ex],tol)
#     plt.savefig(figname,dpi=150,bbox_inches='tight')
#     #plt.show()


# #%% Plot Frequency by Month

# mons3    = proc.get_monstr()
# fsz_axis = 12

# binedges = np.arange(0.5,13.5,1)
# imons =  np.arange(1,13,1)
# ex = 0


# for ex in range(4):
#     fig,axs = plt.subplots(2,1,figsize=(6,8))
    
#     # Plot Nino
#     ax      = axs[0]
#     plotvar = ninomons[ex]
#     _,_,bpnino  = ax.hist(plotvar,bins =binedges,color='firebrick',alpha=0.5,edgecolor='w')
#     ax.set_xticks(imons,labels=mons3)
#     ax.grid(True,which='both',ls='dotted',lw=0.50,c='gray')
#     ax.bar_label(bpnino,fmt="%02i",c='lightgray',fontsize=fsz_axis,label_type='center')
#     ax.set_ylabel(r"Count (El $Ni\tilde{n}o$)")
#     if ex > 1:
#         ax.set_ylim([0,4])
#     else:
#         ax.set_ylim([0,13])
    
#     # Plot Nina
#     ax      = axs[1]
#     plotvar = ninamons[ex]
#     _,_,bpnina = ax.hist(plotvar,bins =binedges,color='cornflowerblue',alpha=0.5,edgecolor='w')
#     ax.set_xticks(imons,labels=mons3)
#     ax.grid(True,which='both',ls='dotted',lw=0.50,c='gray')
#     ax.bar_label(bpnina,fmt="%02i",c='w',fontsize=fsz_axis,label_type='center')
#     #ax.set_title(#$Ni\tilde{n}a$: %i" % (len(ninos),len(ninas)))
#     ax.set_ylabel(r"Count (La $Ni\tilde{n}a$)")
#     if ex > 1:
#         ax.set_ylim([0,4])
#     else:
#         ax.set_ylim([0,13])
    
#     plt.suptitle("AWI-CM3 (%s), nyears=%i" % (expnames_long[ex],expyrs[ex]),y=.90)
    
    
#     # Save the figure
#     figname = "%sNino_Events_MonthDistribution_%s_tol%02imon.png" % (figpath,expnames[ex],tol)
#     plt.savefig(figname,dpi=150,bbox_inches='tight')

# #%% Make Spaghetti Plot Centered around ENSO events

# ninotype = "nina"

# for ninotype in ["nino","nina"]:
    
#     if ninotype == "nino":
#         indict = ninodicts
#     elif ninotype == "nina":
#         indict = ninadicts
        
#     for ex in range(nexps):
#         expdict     = indict[ex]
        
#         ibefore = 12
#         iafter  = 36
        
#         plotlags    = np.hstack([np.flip((np.arange(0,ibefore+1) * -1)),np.arange(1,iafter+1,1)])
        
        
#         target_var  = ensoids[ex].data
        
        
        
#         eventids    = expdict['center_ids']
#         nevents     = len(eventids)
        
        
#         # plotlags       = np.hstack([np.flip((np.arange(0,ibefore+1) * -1)),np.arange(1,iafter+1,1)])
#         # stacked_events = np.zeros((nevents,len(plotlags))) * np.nan
#         # ntime          = len(target_var)
#         # for ie in range(nevents):
            
#         #     ievent    = eventids[ie]
#         #     istart    = ievent-ibefore
#         #     iend      = ievent+iafter
        
            
#         #     if (istart >=0) and (iend < ntime):
#         #         stacked_events[ie,:] = target_var[istart:(iend+1)]
                
#         #     elif iend >= ntime:
#         #         filler = np.zeros( (iend-ntime+1)) * np.nan
#         #         subset = np.hstack([target_var[istart:],filler])
#         #         stacked_events[ie,:] = subset
#         #     elif istart < 0: # Note havent tested this
#         #         filler  = np.zeros(np.abs(istart)) * np.nan
#         #         subset  = np.hstack([filler,target_var[:(iend+1)],])
#         #         stacked_events[ie,:] = subset
#         stacked_events = stack_events(target_var,eventids,ibefore,iafter)
        
#         #%% Now Plot Each one
        
#         xtkslag     = np.arange(-12,39,3)
        
#         fig,ax      = plt.subplots(1,1,constrained_layout=True,figsize=(8,4.5))
        
#         #ax.set_prop_cycle(color=list(mpl.colormaps['Accent'].colors))
        
#         sigma       = ensoids[ex]['std'].item()
        
#         for ie in range(nevents):
            
#             if ie == 0:
#                 lab = "Individual Event"
#             else:
#                 lab = ""
                
#             plotvar = stacked_events[ie,:]
#             ax.plot(plotlags,plotvar,lw=.80,alpha=0.55,label=lab)
            
#         ax.plot(plotlags,np.nanmean(stacked_events,0),color="k",lw=1,alpha=1,label='Mean')
#         ax.legend()
        
#         ax.axvline([0],ls='solid',color="k",lw=0.50)
#         ax.axhline([0],ls='solid',color="k",lw=0.50)
#         ax.axhline([sigma],ls='dotted',color="red",lw=0.75)
#         ax.axhline([-sigma],ls='dotted',color="blue",lw=0.75)
        
#         ax.set_xlabel("Lags (Months)")
#         ax.set_ylabel("SST Anomaly $\degree C$")
#         if expnames[ex] == "TCo319_ssp585":
#             ax.set_ylim([-7.5,7.5])
#         else:
#             ax.set_ylim([-4.5,4.5])
#         ax.set_xticks(xtkslag)
#         ax.set_xlim([plotlags[0],plotlags[-1]])
        
#         ax.set_title("%s Composites (%s)\n$\sigma=$%.2f" %(ninotype,expnames_long[ex],sigma))
        
        
#         # Save the figure
#         figname = "%s%s_Events_Spaghetti_%s_tol%02imon.png" % (figpath,ninotype,expnames[ex],tol)
#         plt.savefig(figname,dpi=150,bbox_inches='tight')
        
#         #plt.show()

# #%% Get Stacked Events

# stacked_events = []
# for ninotype in ["nino","nina"]:
        
#     if ninotype == "nino":
#         indict = ninodicts
#     elif ninotype == "nina":
#         indict = ninadicts
    
#     exstack = []
#     for ex in range(nexps):
#         expdict     = indict[ex]
#         print(ex)
#         ibefore = 12
#         iafter  = 36
#         plotlags    = np.hstack([np.flip((np.arange(0,ibefore+1) * -1)),np.arange(1,iafter+1,1)])
#         target_var  = ensoids[ex].data
#         eventids    = expdict['center_ids']
        
#         exstack.append(stack_events(target_var,eventids,ibefore,iafter))
#     stacked_events.append(exstack)
        
# #%% Plot the average Response


# fig,ax = plt.subplots(1,1,constrained_layout=True)
# for nn in range(2):
    
#     for ex in range(nexps):
        
#         plotvar = np.nanmean(stacked_events[nn][ex],0)
#         ax.plot(plotlags,plotvar,label=expnames[ex],lw=2.5)
# ax.legend()

# ax.axvline([0],ls='solid',color="k",lw=0.50)
# ax.axhline([0],ls='solid',color="k",lw=0.50)
# #ax.axhline([sigma],ls='dotted',color="red",lw=0.75)
# #ax.axhline([-sigma],ls='dotted',color="blue",lw=0.75)

# ax.set_xlabel("Lags (Months)")
# ax.set_ylabel("SST Anomaly $\degree C$")
# ax.set_ylim([-4.25,4.25])
# ax.set_xticks(xtkslag)
# ax.set_xlim([plotlags[0],plotlags[-1]])
    
# plt.show()

# #%% Plot the standard Deviations


# expcols = ["cornflowerblue",'lightcoral',
#            "slateblue","firebrick",
#            "midnightblue","k"]

# ensostdevs = np.array([es['std'].data.item() for es in ensoids])

# fig,ax = plt.subplots(1,1,sharey=True,constrained_layout=True,
#                       figsize=(8,4.5))


# plotexps = [0,1,2,3,4,5]

# braw = ax.bar(np.array(expnames_long)[plotexps],
#         ensostdevs[plotexps],color=np.array(expcols)[plotexps]
#         )

# ax.bar_label(braw,fmt="%.02f",c='white',fontsize=fsz_axis,label_type='center')

# ax.set_ylabel(r"1$\sigma$ $Ni\tilde{n}o3.4$ [$\degree C$]")
# ax.set_ylim([0,2])

# # ax       = axs[1]
# # plotexps = [1,3]
# # ax.bar(np.array(expnames_long)[plotexps],
# #        ensostdevs[plotexps],color=np.array(expcols)[plotexps]
# #        )
# ax.spines[['right', 'top']].set_visible(False)

# figname = "%sENSO_Thresholds_%s.png" % (figpath,ninoid_name)
# plt.savefig(figname,dpi=150,bbox_inches='tight')

# plt.show()
# ==================================================
#%% Temprarily commenting above to visualize events

#%% Examine the duration of ENSO events and how they have changed


    
    
    
    


# def get_ids_stacker():
    

# def stack_events(eventids,targetarr,ibefore,iafter,axis=0,):
    
#     nevents        = len(eventids)
#     plotlags       = np.hstack([np.flip((np.arange(0,ibefore+1) * -1)),np.arange(1,iafter+1,1)])
#     stacked_events = np.zeros((nevents,len(plotlags))) * np.nan
#     ntime          = targetarr.shape[axis]
    
#     for ie in range(nevents):
        
#         ievent    = eventids[ie]
#         istart    = ievent-ibefore
#         iend      = ievent+iafter
        
#         # Corrections for end-case indices
#         if istart < 0:
#             print("istart is at %s" % istart)
#             insert_start = np.abs(istart)#(lags - np.abs(istart)).item()
#             istart       = 0
#         else:
#             insert_start = 0
        
#         if iend > ntime:
#             print("iend is at %s" % (iend))
#             insert_end = lags + (ntime-ievent) #(lags*2+1) - (iend - ntime)
#             iend = ntime
#         else:
#             insert_end = lags*2+1
            
        
#         indices_in = np.arange(insert_start,insert_end)
#         temp_var[ie,indices_in,:,:] = invar[istart:(iend+1),:,:]
#         event_times[ie,indices_in]  = times_da[istart:(iend+1)]
        
    
    

    
# Part 2 ====================================================================== 
#%% Lets Load a variable to analyze

vname           = 'sst'#"Dmaxgrad"
ds_var          = []

for ex in tqdm.tqdm(range(nexps)):
    
    
    if vname == "tx_sur":
        vname_file = "tx_surf_1m"
        ncname = "%s%s_%s.nc" % (datpath,expnames[ex],vname_file)
    else:
        ncname = "%s%s_%s.nc" % (datpath,expnames[ex],vname)
    ds = xr.open_dataset(ncname).load()
    
    if vname.upper() in list(ds.keys()):
        print("Renaming %s to %s" % (vname.upper(),vname))
        ds = ds.rename({vname.upper():vname})
    
    ds_var.append(ds)

ds_var      = [ut.standardize_names(ds) for ds in ds_var]
if vname in ["D20","Dmaxgrad"]:
    ds_anoms = [preprocess_enso(ds['nz1']) for ds in ds_var]
else:
    ds_anoms = [preprocess_enso(ds[vname]) for ds in ds_var]

#%% Examine Lag Compsites of choice variable before and after an event

# User Settings
plot_single_events  = True # Set to True to make plots for individual events
lags                = 18    # Number of lead + lags to composite over
vmax_event          = 2.5    # Colorbar Limits for Single Events
vmax_composite      = 2.5     # Colorbar limits for composites
sel_mons            = None  # [12,1,2] # Indicate which central months to include in the composite
plot_composite      = True  # Set to True to make plots for composite

# Adjust vmax based on variable
if vname == "ssr" or vname == "str":
    vmax_event      = 50
    vmax_composite  = 40
elif vname == "tx_sur":
    vmax_event      = 0.02
    vmax_composite  = 0.01
elif vname == "D20" or vname == "Dmaxgrad":
    vmax_event = 30
    vmax_composite=30
elif vname == "eis":
    vmax_event     = 7.5
    vmax_composite = 5
elif vname == "sst":
    vmax_event     = 2.5
    vmax_composite = 2.5


#Make some necessary variables
leadlags            = np.arange(-lags,lags+1)
proj                = ccrs.PlateCarree()

composites_byexp    = []
for ex in range(nexps):#[0,1,2,3]:
    
    # Open the Nino/Nina Indices and Variables
    ninoids   = ninodicts[ex]['center_ids']
    ninotimes = ninodicts[ex]['event_time']
    #ninomon   = [ds.time.dt.month.item() for ds in ninotimes]
    
    ninaids   = ninadicts[ex]['center_ids']
    ninatimes = ninadicts[ex]['event_time']
    #ninamon   = [ds.time.dt.month.item() for ds in ninatimes]
    
    invar    = ds_anoms[ex].transpose('time','lat','lon')
    times_da = invar.time.data.astype('datetime64[ns]')
    
    # Select/align qualifying Events and prepare to composite
    
    composites_bytype = []
    for eventid_type in ['nino','nina']:
        
        if eventid_type == "nina":
            eventids_in = ninaids
        elif eventid_type == "nino":
            eventids_in = ninoids
            
        
        #%% PREALLOCATE AND Subset Events
        ntime,nlat,nlon = invar.shape
        nevents     = len(eventids_in)
        temp_var    = np.zeros((nevents,lags*2+1,nlat,nlon)) * np.nan
        event_times = np.zeros((nevents,lags*2+1)) * np.nan
                               
        for ie in range(nevents):
            
            ievent    = eventids_in[ie]
            istart    = ievent-lags
            iend      = ievent+lags
            
            # Corrections for end-case indices
            if istart < 0:
                print("istart is at %s" % istart)
                insert_start = np.abs(istart)#(lags - np.abs(istart)).item()
                istart       = 0
            else:
                insert_start = 0
            
            if iend > ntime:
                print("iend is at %s" % (iend))
                insert_end = lags + (ntime-ievent) #(lags*2+1) - (iend - ntime)
                iend = ntime
            else:
                insert_end = lags*2+1
                
            
            indices_in = np.arange(insert_start,insert_end)
            temp_var[ie,indices_in,:,:] = invar[istart:(iend+1),:,:]
            event_times[ie,indices_in]  = times_da[istart:(iend+1)]
        
        
        
        coords             = dict(eventid=eventids_in,lag=leadlags,lat=invar.lat,lon=invar.lon)
        invar_subset       = xr.DataArray(temp_var,coords=coords,dims=coords)
        
        coords_time        = dict(eventid=eventids_in,lag=leadlags)
        event_times_subset = xr.DataArray(event_times.astype('datetime64[ns]'),coords=coords_time,dims=coords_time)
        
        #% Plot a sequence to test (generate all frames)----------------------
        if plot_single_events:
            plotlags        = leadlags #[-18,-12,-6,0,6,12,18]
            vmax            = vmax_event
            
            # Looping for each Event
            for ie in range(nevents):
                
                # Select Event and make folder
                eventnum    = invar_subset.isel(eventid=ie).eventid.item()
                eventtime   = str(event_times_subset.sel(lag=0).isel(eventid=ie).data)[:7]
                eventdir    = "%s/pre_composites/%s/%s/%s/event_%04i_%s/" % (figpath,expnames[ex],eventid_type,vname,eventnum,eventtime)
                proc.makedir(eventdir)
                
                for il in tqdm.tqdm(range(len(plotlags))):
                    lag     = plotlags[il]
                    
                    fig,ax  = init_tp_map()
                    
                    seltime = str(event_times_subset.sel(lag=lag).isel(eventid=ie).data)[:7]
                    plotvar = invar_subset.isel(eventid=ie).sel(lag=lag)
                    ax.set_title("AWI-CM3 %s Event %i, Lag %i (%s)" % (expnames_long[ex],plotvar.eventid.data.item(),lag,seltime))
                    
                    pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,vmin=-vmax,vmax=vmax,cmap='cmo.balance',transform=proj)
                    cb      = fig.colorbar(pcm,ax=ax,fraction=0.010,pad=0.01)
                    
                    # Save the figure
                    figname = "%sframe%0i.png" % (eventdir,il)
                    plt.savefig(figname,dpi=150,bbox_inches='tight')
                    
                    plt.close()
                    #plt.show()
        
        # Next Part: Actually make the composites =============================
        if eventid_type == "nino":
            eventmons        = ninomons[ex]
        elif eventid_type == "nina":
            eventmons        = ninamons[ex]
        if sel_mons is None:
            
            selected_events  = np.arange(nevents)#eventmons # Just Take All Months

        else:
            # Note this might be wrong
            selected_events  = [m in sel_mons for m in eventmons]
        
        # Make the Composite
        composite_events = temp_var[selected_events,:,:,:]
        composite_out    = np.nanmean(composite_events,0)
        
        composites_bytype.append(composite_out)
        # Optionally plot them as well =============================
        #plotlags = leadlags #[-18,-12,-6,0,6,12,18]
        if plot_composite:
            
            # Prepare and make composite directory
            lon,lat  = ds_anoms[ex].lon,ds_anoms[ex].lat
            eventdir  = "%s/pre_composites/%s/%s/%s/composites/" % (figpath,expnames[ex],eventid_type,vname)
            proc.makedir(eventdir)
            
            vmax = vmax_composite

            for il in tqdm.tqdm(range(len(plotlags))):
                lag     = plotlags[il]
                ilag    = np.where(leadlags == lag)[0][0].item()
                
                
                fig,ax  = init_tp_map()
                
                plotvar = composite_out[il,:,:]
                ax.set_title("AWI-CM3 %s %s Composite\nEvent Month = %s, # Events = %i, Lag %i" % (expnames_long[ex],eventid_type,sel_mons,composite_events.shape[0],lag))
                
                pcm     = ax.pcolormesh(lon,lat,plotvar,vmin=-vmax,vmax=vmax,cmap='cmo.balance',transform=proj)
                cb      = fig.colorbar(pcm,ax=ax,fraction=0.010,pad=0.01)
                
                #plt.show()
                
                
                
                # Save the figure
                figname = "%sENSO_Composite_%s_%s_frame%0i.png" % (eventdir,expnames[ex],eventid_type,il)
                plt.savefig(figname,dpi=150,bbox_inches='tight')
                
                plt.close()
                #plt.show()
    composites_byexp.append(composites_bytype)
            
#%% Save the Output

outpath = "/home/niu4/gliu8/projects/scrap/TP_crop/composites/"
for ex in range(nexps):
    
    outname = "%sENSO_Composites_%s_%s.nc" % (outpath,vname,expnames[ex])
    
    coords  = dict(event=["nino","nina"],lags=leadlags,lat=ds_anoms[ex].lat,lon=ds_anoms[ex].lon)
    arrout  = np.array(composites_byexp[ex])
    da_out  = xr.DataArray(arrout,coords=coords,dims=coords,name=vname)
    
    encoding=proc.make_encoding_dict(da_out)
    da_out.to_netcdf(outname,encoding=encoding)

# # ============================================================================================================================================================
# #%%Scrap Below


# #%% First lets make a composite agnostic to month (i know... not ideal)

# # Run the above to extract the events




# print(len(eventmons) == len(ninoids))






# #%% Try PLotting the sequence of selected lags



# #for il in range(len(plotlags)):













# #%% Check Work





#     # if iend > ntime:
#     #     iend   = ntime
#     # if istart < 0:
#     #     istart = 0
        
    
#     varbefore = invar.data[istart:ievent,:,:]
#     varafter  = invar.data[ievent:(iend+1),:,:]
    
#     temp_var[ie,:,nlat,nlon] = np.concatenate([varbefore,varafter],axis=0)
    


# #%% Try to figure things out
    
# edummy = 2 
# ntime  = 16
# lags   = 18

# dummyvar = np.zeros((lags*2+1,1,1)) * np.nan
# loopy    = np.arange(ntime)#np.random.normal(0,1,ntime)

# for ie in range(nevents):
    
#     ievent    = edummy
#     istart    = ievent-lags
#     iend      = ievent+lags
    
#     if istart < 0:
#         print("istart is at %s" % istart)
#         insert_start = np.abs(istart)#(lags - np.abs(istart)).item()
#         istart       = 0
#     else:
#         insert_start = 0
    
#     if iend > ntime:
#         print("iend is at %s" % (iend))
#         insert_end = lags + (ntime-ievent) #(lags*2+1) - (iend - ntime)
#         iend = ntime
#     else:
#         insert_end = lags*2+1
        
#     indices_in = np.arange(insert_start,insert_end)
#     dummyvar[indices_in,[0],[0]] = loopy[istart:(iend+1)]
        
        
        
# #%%

# lags = np.arange(-18,19)
# lagid = np.arange(len(lags))

# fig,ax = plt.subplots(2,1)
# ax[0].set_xticks(lags)
# ax[1].set_xticks(lagid)

# plt.show()


#%%
    
# for ie in range(ncevent):
#     idtime = center_ids[ie]
#     ax.axvline(ensoin.time.isel(time=idtime),c='r',ls='dotted')




    
    
#     if (ii == 0) or (new_event is True):
#         event_merge = [] # Make Array for New Event
    
    
    
#     if (ievent - prev_id) > 1 
    
        
        
#     event_end = False
    
    
#     if ii > 0:
#         if (ievent - prev_id) > 1:
#             event_combine_
        
    


# num2next = np.roll(eventid,-1) - eventid # Time to next event
# num2next[-1] = 0


    



#%% 


