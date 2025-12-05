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

TP_crop         = False # Set True to calculate using Tropical Pacific Values
expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"]


# Just look at the 4th experiment
ex              = 0
# timecrops       = [[1950,2100],None,None,None,None]

# Only available for the 5km simulation
datpath         = "/home/niu4/gliu8/projects/scrap/TP_crop/"
datpath_glob    = "/home/niu4/gliu8/projects/scrap/processed_global/"
vnames          = ["skt","pl_t_700"]
expname         = "TCo319_ctl1950d" #"TCo1279-DART-2090" #"TCo1279-DART-1950" #"TCo2559-DART-1950C"

processed_path = False

if expname == "TCo2559-DART-1950C":
    vnames      = ["skt","pl_t_700"]
    
elif expname == "TCo1279-DART-1950":
    vnames          = ["skt","t700"]
    
elif expname == "TCo319_ctl1950d" or expname == "TCo1279-DART-2090":
    processed_path = True
    
    vnames      = ["sst","t700"] # Use SST as it is the same thing as skin temperature and you want ocean-only points
    datpath     = '/home/niu4/gliu8/projects/scrap/processed_global/'
    datpath_glob = '/home/niu4/gliu8/projects/scrap/processed_global/'
    
    

    
    
    

#%% First, just try for the tropical pacific
# need to load in T700 and T0 (skt is ok?)

chunkdict = dict(
    lat  = 'auto',#5120/40,
    lon  = 'auto',#10256/16,#16
    time_counter = 'auto',#120/30 
    )

dsall = []
if TP_crop: # Load just the tropical version
    for v in range(2):
        ncname = "%s%s_%s.nc" % (datpath,expname,vnames[v])
        ds = xr.open_dataset(ncname).load()
        dsall.append(ds)
else:
    
    
    
    for v in range(2):
        vname  = vnames[v]
        if processed_path:
            ncname = "%s%s_%s.nc" % (datpath_glob,expname,vname)
        else:
            ncname = ut.get_rawpath_awi(expname,vname,ensnum=None)[0]
        ds = xr.open_dataset(ncname)#chunks=chunkdict)#.load()
        dsall.append(ds)


sktname     = vnames[0]
t700name    = vnames[1]

dsall       = [ds.chunk(chunkdict) for ds in dsall]
skt,t700    = dsall
skt         = skt[sktname]
t700        = t700[t700name].squeeze()

# Make sure time matches up
skt,t700 = proc.match_time_month(skt,t700,timename='time_counter')


print(np.any(t700 < 0))
print(np.any(skt < 0))

#%% Check they are in Kelvin

# Calculate EIS
st     = time.time()
eis    = climlab.thermo.EIS(skt.data,t700.data)

# Read back into DataArray
coords = dict(time=skt.time_counter.data,lat=skt.lat,lon=skt.lon)#skt.coords
ds_eis = xr.DataArray(eis,coords=coords,dims=coords,name="eis")
edict  = proc.make_encoding_dict(ds_eis)

# Load and do computation (552 Seconds for Global)
if not TP_crop:
    ds_eis = ds_eis.compute()
print("Calculated in %2.fs" % (time.time()-st))

# Set Output and Write 1235 sec
st = time.time()
if TP_crop:
    outname = #"%s%s_eis_useskt.nc" % (datpath,expname)
else:
    outname = "%s%s_eis.nc" % (datpath_glob,expname) #"%s%s_eis_useskt.nc" % (datpath_glob,expname)
ds_eis.to_netcdf(outname,encoding=edict)
print("Saved in %2.fs" % (time.time()-st))

# #%% Reload if already calculated

# outname = "%s%s_eis_useskt.nc" % (datpath,expname)
# ds_eis = xr.open_dataset(outname).eis.load()

# #%%
# eismean = ds_eis.mean('time')

# #%% Visualize the time-mean eis

# proj    = ccrs.PlateCarree()
# cints   = np.arange(-8,8.5,0.5)
# fig,ax  = ut.init_tp_map()
# plotvar = eismean
# pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,vmin=cints[0],vmax=cints[-1],cmap='cmo.balance')
# cl      = ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=cints,
#                      colors="k",linewidths=0.75)

# cb      = viz.hcbar(pcm,ax=ax,fraction=0.045,pad=0.1)
# cb.set_label("Mean Estimated Inversion Strength (%s)" % expnames[ex])
# ax.clabel(cl,levels=cints[::2],fontsize=8)

# plt.show()

# #%% Load in ENSO Indices

# ensoid   = ut.load_ensoid(expname,"nino34",standardize=False)


