#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calc EIS ERA5

calculate EIS in era5

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

datpath = '/home/niu4/gliu8/share/ERA5/processed/'
outname = datpath+"eis_1979_2024.nc"
vnames  = ['skt','t700']
ncnames = [
    "skt_1979_2024.nc",
    "t700_1979_2024.nc",
    ]

timename = 'valid_time'
latname  = 'latitude'
lonname  = 'longitude'


    

#%% First, just try for the tropical pacific
# need to load in T700 and T0 (skt is ok?)

dsall = []
for v in range(2):
    vname  = vnames[v]
    ncname = datpath+ncnames[v]#ut.get_rawpath_awi(expname,vname,ensnum=None)[0]
    ds = xr.open_dataset(ncname)[vname]#chunks=chunkdict)#.load()
    dsall.append(ds)


#dsall       = [ds.chunk(chunkdict) for ds in dsall]
dsall       = [ds.rename({timename: 'time'}) for ds in dsall]

if len(dsall[0].time) != len(dsall[0].time):
    dsall       = proc.match_time_month(dsall[0],dsall[1])
    rename_time=True
    timename_old = timename
    timename = 'time'
else:
    rename_time=False
    
skt,t700    = dsall
skt         = skt.squeeze()#.skt
t700        = t700.squeeze()



print(np.any(t700 < 0))
print(np.any(skt < 0))


#%% Check they are in Kelvin



# Calculate EIS
st     = time.time()
eis    = climlab.thermo.EIS(skt.data,t700.data)

# Read back into DataArray
coords = {timename : skt[timename].data,
          latname  : skt[latname].data,
          lonname :  skt[lonname].data}

ds_eis = xr.DataArray(eis,coords=coords,dims=coords,name="eis")
ds_eis = ds_eis.rename({timename:timename_old})
edict  = proc.make_encoding_dict(ds_eis)

# # Load and do computation (552 Seconds for Global)
# if not TP_crop:
#     ds_eis = ds_eis.compute()
# print("Calculated in %2.fs" % (time.time()-st))

st = time.time()
ds_eis.to_netcdf(outname,encoding=edict)
print("Saved in %2.fs" % (time.time()-st))

