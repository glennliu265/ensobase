#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Estimated Inversion strength in AWI-CM3 simulations
using the function from climlab

Created on Thu Oct  9 10:47:15 2025

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

#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% Datpaths

expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"]


# Just look at the 4th experiment
ex              = 4
# timecrops       = [[1950,2100],None,None,None,None]

# Only available for the 5km simulation
datpath         = "/home/niu4/gliu8/projects/scrap/TP_crop/"
vnames          = ["skt","pl_t_700"]
expname         = "TCo2559-DART-1950C"



#%% First, just try for the tropical pacific
# need to load in T700 and T0 (skt is ok?)

dsall = []
for v in range(2):
    ncname = "%s%s_%s.nc" % (datpath,expname,vnames[v])
    ds = xr.open_dataset(ncname).load()
    dsall.append(ds)

skt,t700    = dsall
skt         = skt.skt
t700        = t700.t.squeeze()


print(np.any(t700 < 0))
print(np.any(skt < 0))

#%% Check they are in Kelvin

eis = climlab.thermo.EIS(skt.data,t700.data)

coords = dict(time=skt.time_counter.data,lat=skt.lat,lon=skt.lon)#skt.coords
ds_eis = xr.DataArray(eis,coords=coords,dims=coords,name="eis")

edict  = proc.make_encoding_dict(ds_eis)

outname = "%s%s_eis_useskt.nc" % (datpath,expname)
ds_eis.to_netcdf(outname,encoding=edict)


#%% Reload if already calculated

outname = "%s%s_eis_useskt.nc" % (datpath,expname)
ds_eis = xr.open_dataset(outname).eis.load()




#%%
eismean = ds_eis.mean('time')

#%% Visualize the time-mean eis

proj    = ccrs.PlateCarree()
cints   = np.arange(-8,8.5,0.5)
fig,ax  = ut.init_tp_map()
plotvar = eismean
pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,vmin=cints[0],vmax=cints[-1],cmap='cmo.balance')
cl      = ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=cints,
                     colors="k",linewidths=0.75)

cb      = viz.hcbar(pcm,ax=ax,fraction=0.045,pad=0.1)
cb.set_label("Mean Estimated Inversion Strength (%s)" % expnames[ex])
ax.clabel(cl,levels=cints[::2],fontsize=8)

plt.show()

#%% Load in ENSO Indices

ensoid   = ut.load_ensoid(expname,"nino34",standardize=False)


