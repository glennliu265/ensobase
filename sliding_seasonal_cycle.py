#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Analyze Phase Change of seasonal cycle in variables

Check First for Nino3.4 in SSP5.85 Simulation

Created on Fri May 15 11:08:55 2026

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
from tqdm import tqdm

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut


#%% Load ENSO Index

ninoid_name = "nino34"
expname     = "TCo319_ssp585"
dsnino      = ut.load_ensoid(expname,ninoid_name=ninoid_name,standardize=False)

#%% Calculating Sliding Metric


def sliding_calc_xr(dsraw,nyr_window,func,detrend=False,convert_arr=False,preproc=True):
    # Adapted from sliding_pointwise_spectra
    # Take in dsraw, subset into periods, detrend, then compute spectra
    
    subsets,tranges = ut.generate_periods(dsraw,nyr_window)
    nperiods        = len(subsets) # Get Number of Periods
    tcenters        = [ut.get_center_time(t) for t in tranges]
    
    # Preprocess by Period
    if preproc:
        subsets_anom   = ut.preprocess_byperiod(subsets,verbose=False,detrend=detrend) # Don't Detrend...
    else:
        subsets_anom   = subsets
    # Calculate Power Spectra for each period
    calc_byperiod  = []
    for nw in range(nperiods):
        sample_in            = subsets_anom[nw]#.data
        if convert_arr:
            sample_in        = sample_in.data
        calcout              = func(sample_in) 
        calc_byperiod.append(calcout)
    
    # Concatenate By period (works)
    calc_byperiod             = xr.concat(calc_byperiod,dim='period')
    
    # Make Into DataArray
    coords     = dict(period=np.arange(nperiods))
    #da_out     = xr.DataArray(np.array(r2_byperiod),dims=coords,coords=coords,name=vname)
    da_tcenter = xr.DataArray(tcenters,coords=coords,dims=coords,name="tcenter")
    da_out     = xr.merge([calc_byperiod,da_tcenter])
    
    return da_out




calc_scycle = lambda ds : ds.groupby('time.month').mean('time')
calc_monstd = lambda ds : ds.groupby('time.month').std('time')


#%%

nyr_window = 30
func       = calc_monstd
dsraw      = dsnino

slideout = sliding_calc_xr(dsraw,nyr_window,func)

#%% Make an early/late period plot

nperiods = len(slideout.period)

cm.YlGnBu



#%% Plot Early and Late

mons3 = proc.get_monstr()
pp    = 10

fig,ax = viz.init_monplot(1,1)

plotvar = slideout.sst.isel(period=np.arange(-pp,0)).mean('period')
ax.plot(mons3,plotvar,c="hotpink",label="Latest %i Periods" % pp)

plotvar = slideout.sst.isel(period=np.arange(0,pp)).mean('period')
ax.plot(mons3,plotvar,c="midnightblue",label="Earliest %i Periods" % pp)

plotvar = slideout.sst.mean('period') #isel(period=-1)
ax.plot(mons3,plotvar,c="k",label="All Period Average")

ax.legend()
ax.set_ylabel("Monthly Stdev, Nino3.4 Index ($\degree$C)")


plt.show()


#%% Load AWI-CM3 Variables


expname = "TCo319_ssp585"
vnames  = ["sst","cre","tscre","ttcre"]
vunits  = ["\degree C","W m^{-2}","W m^{-2}","W m^{-2}"]
datpath = "/home/niu4/gliu8/projects/scrap/regrid_1x1/"
nvars   = len(vnames)

dsvars  = []
for vv in range(nvars):
    vname   = vnames[vv]
    ncname  = "%s%s_%s_regrid1x1.nc" % (datpath,expname,vname)
    ds      = xr.open_dataset(ncname)[vname].load()
    ds      = ut.standardize_names(ds)
    ds      = ut.varcheck(ds,vname,expname)
    dsvars.append(ds)


#%%
bbox_sep    = [-90+360,-75+360,-40,-15]
sepavg = lambda ds: proc.area_avg_cosweight(proc.sel_region_xr(ds,bbox_sep))


dsreg = [sepavg(ds) for ds in dsvars]

vv    = 2
dsin = dsreg[vv]

#%% Do pointwise sinusoidal fit

t = np.arange(len(dsin.time))

ssout = proc.fit_sinfunc(t,dsin,fix_freq=2*np.pi/12)



def sliding_sinfit(dsraw,nyr_window,detrend=False,fix_freq=2*np.pi/12):
    # Adapted from sliding_pointwise_spectra
    # Take in dsraw, subset into periods, detrend, then compute spectra
    
    subsets,tranges = ut.generate_periods(dsraw,nyr_window)
    nperiods        = len(subsets) # Get Number of Periods
    tcenters        = [ut.get_center_time(t) for t in tranges]
    
    # Calculate By Period
    calc_byperiod  = []
    for nw in range(nperiods):
        
        sample_in            = subsets[nw]#.data
        t                    = np.arange(len(sample_in))
        
        calcout = proc.fit_sinfunc(t,sample_in,fix_freq=fix_freq)
        calc_byperiod.append(calcout)
    
    # Concatenate By period (works)
    #calc_byperiod             = xr.concat(calc_byperiod,dim='period')
    
    # Make Into DataArray
    coords     = dict(period=np.arange(nperiods))
    #da_out     = xr.DataArray(np.array(r2_byperiod),dims=coords,coords=coords,name=vname)
    da_tcenter = xr.DataArray(tcenters,coords=coords,dims=coords,name="tcenter")

    return calc_byperiod,subsets,da_tcenter


slideout_sinfit = sliding_sinfit(dsin,nyr_window,)

calc_byperiod,subsets,tcenters = slideout_sinfit

#%% Plot Change in Period

calc_byperiod = slideout_sinfit[0]
phases        = np.array([dd['phase'] for dd in calc_byperiod])
amplitudes    = np.array([dd['amplitude'] for dd in calc_byperiod])
frequencies   = np.array([dd['frequency'] for dd in calc_byperiod])
offset        = np.array([dd['offset'] for dd in calc_byperiod])


#%%

# Plot Change in Sine Fit
fig,axs = plt.subplots(2,1)

# Plot Offset
ax  = axs[0]
ax.plot(tcenters,offset,c='cornflowerblue')
ax.set_ylabel("Offset (D)")
ax  = viz.change_axcol('left','cornflowerblue',ax=ax)
ax1 = ax.twinx()
ax1.plot(tcenters,amplitudes,c='red')
ax1.set_ylabel("Amplitude ($A$)")
ax1 = viz.change_axcol('right','red',ax=ax1)

# 
ax = axs[1]
ax.plot(tcenters,np.rad2deg(phases),c='forestgreen')
ax.set_ylabel("Phase Shift (C)")
ax  = viz.change_axcol('left','forestgreen',ax=ax)
ax1 = ax.twinx()
ax1.plot(tcenters,np.rad2deg(frequencies),c='orange')
ax1.set_ylabel("Frequency (B)")
ax1 = viz.change_axcol('right','orange',ax=ax1)

# 
for ax in axs:
    ax.set_xticks(tcenters[::10])
plt.suptitle("A * Sin (Bt + C) + D")

plt.show()

ax.plot(tcenters,,label="Phase",marker="o")
plt.show()


