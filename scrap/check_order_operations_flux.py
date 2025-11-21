#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Check Calculation Sensitivity 

Examine the following at a point

- ttr/tsr CRE from anomaly vs. from full term
- CRE from difference all and clear sky vs from adding ttcre and tscre


Created on Tue Oct 28 10:03:43 2025

@author: 
    
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

#%%
lonf    = 330
latf    = 50

#%%

dpath  = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_control/"
vnames = ['tsr','tsrc','ttr','ttrc']

dspts = []
for vname in vnames:
    nclist = glob.glob(dpath + "*_%s*.nc" % vname)[0]
    ds     = xr.open_dataset(nclist)[vname]
    dspt   = proc.selpt_ds(ds,lonf,latf)
    dspts.append(dspt)
    
    
dspts = [ds.load() for ds in dspts]
    

    
dspts            = [ut.standardize_names(ds) for ds in dspts]
dspts            = [ut.remove_duplicate_times(ds) for ds in dspts]

scycles          = xr.merge([ds.groupby('time.month').mean('time') for ds in dspts])


dsanoms          = xr.merge([proc.xrdeseason(ds) for ds in dspts],join='override')
dspts            = xr.merge(dspts,join='override')

allsky           = dspts.ttr + dspts.tsr
clearsky         = dspts.ttrc + dspts.tsrc




#%% Part 1

"""
Anomalize First, then compute allsky vs. compute allsky first

It seems that there is a slight difference

- Winter values are more positive
- Summer values are more negative...


My instinct would thus be to compute allsky first, then perform any preprocessing...

"""


allsky_anomfirst = dsanoms.ttr + dsanoms.tsr
allsky_first     = proc.xrdeseason(allsky)


scycle_asf  = proc.calc_clim(allsky.data,0)
asprime     = proc.deseason(allsky.data,dim=0)

allskyrs    = allsky.data.reshape(285,12)
scycle      = allskyrs.mean(0)
allskyanom  = allskyrs - scycle[None,:]

allskyanom2 = allskyanom - allskyanom.mean(0)[None,:]



#%% Plot the Difference

fig,ax = plt.subplots(1,1,figsize=(12.5,4.5),constrained_layout=True)

plotvar = allsky_anomfirst
ax.plot(plotvar.time,plotvar,color="blue",label="Anomalize First, then sum")

plotvar = allsky_first
ax.plot(plotvar.time,plotvar,color="red",label="Sum first, then anomalize",ls='dashed')

ax.legend()

#%% 

diff = allsky_anomfirst - allsky_first
plt.plot(diff),plt.show()

diff.groupby('time.month').mean('time').plot(),plt.show()

allsky_anomfirst.groupby('time.month').mean('time').plot(),plt.show()

allsky_first.groupby('time.month').mean('time').plot(),plt.show()

#%% Part 2

"""
Look at Differences in CRE
(compute for TTRE or TSCRE first, then add)
(or just compute separately)

Focus first on the total terms

"""

netcre           = allsky - clearsky
ttcre            = dspts.ttr - dspts.ttrc
tscre            = dspts.tsr - dspts.tsrc
sumcre           = ttcre + tscre

#%% Plot the Diffeerences


fig,ax = plt.subplots(1,1,figsize=(12.5,4.5),constrained_layout=True)

plotvar = netcre
ax.plot(plotvar.time,plotvar,color="blue",label="CRE (allsky + clearsky)")

plotvar = sumcre
ax.plot(plotvar.time,plotvar,color="red",label="CRE (ttcre + tscre)",ls='dashed')

ax.legend()
#ax.set_xlim([0,100])

plt.show()

#%%

diffc = netcre - sumcre
diffc.groupby('time.month').mean('time').plot(),plt.show()

#%% 

net_scycle = netcre.groupby('time.month').mean('time')
sum_scycle = sumcre.groupby('time.month').mean('time')


fig,ax = viz.init_monplot(1,1)
ax.plot(net_scycle,color='c',label="allsky - clearsky")
ax.plot(sum_scycle,color='r',label="ttcre + tscre",ls='dashed')
ax.legend()
plt.show()

#%% Just plot the seasonal cycle of each variable to be sure

cols = ['blue','blue','red','red']
lss  = ['solid','dashed','solid','dashed']

fig,ax = viz.init_monplot(1,1)
for v,vname in enumerate(vnames):
    ax.plot(scycles[vname],label=vname,c=cols[v],ls=lss[v])
    
ax.plot(ttcre.groupby('time.month').mean('time'),label='ttcre',c='red',ls='dotted')
ax.plot(tscre.groupby('time.month').mean('time'),label='tscre',c='blue',ls='dotted')

ax.plot(allsky.groupby('time.month').mean('time'),label='allsky',c='k')
ax.plot(clearsky.groupby('time.month').mean('time'),label='clearsky',c='cornflowerblue')
ax.plot(netcre.groupby('time.month').mean('time'),label='CRE Net',c='gray')

ax.legend()
plt.show()

