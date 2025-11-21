#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

AWI Grid Checker

Created on Thu Nov  6 15:31:13 2025

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

expnames            = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
msktest     = [ut.load_land_mask_awi(expname) for expname in expnames]

#%%

def check_grid(msk,dsin,return_grid=False):
    lon = dsin.lon.data
    lat = dsin.lat.data
    
    lonref = msk.lon
    latref = msk.lat
    
    chklon = lon == lonref
    chklat = lat == latref
    
    if np.all(chklon):
        print("\tLon ok (size=%i)" % len(chklon))
    else:
        print("\tMismatch on lon: %s" % chklon)
        
    if np.all(chklat):
        print("\tLat ok (size=%i)" % len(chklat))
    else:
        print("\tMismatch on lat: %s" %  chklat)
    if return_grid:
        llout = [lon,lat]
        return chklon,chklat,llout
    return chklon,chklat

datpath = "/home/niu4/gliu8/projects/scrap/processed_global/"

#for expname in expnames:
ex       = 2
expname = expnames[ex]
ncsearch = "%s%s*.nc" % (datpath,expname)
nclist   = glob.glob(ncsearch)

bad_nc  = []
lonlats = []
for nc in nclist:
    ds = xr.open_dataset(nc)
    out = check_grid(msktest[ex],ds,return_grid=True)
    clon,clat,lonlat = out
    if ~np.any(clon.data) or ~np.any(clat.data):
        print("%s has mismatched lat.lon" % nc)
        bad_nc.append(nc)
        lonlats.append(lonlat)
        

lonref    = msktest[ex].lon
latref    = msktest[ex].lat
    
targetlon = lonlats[0][0]
targetlat = lonlats[0][1]

difflon   = lonref.data - targetlon.data
difflat   =  targetlat.data - latref.data

    
fig,ax    = plt.subplots(1,1)
ax.plot(lonref,color="k",marker="o")
ax.plot(lonlats[0][0],color="red",marker="o")
plt.show()

#%%

fig,ax = plt.subplots


    #if np.any()
    
    
    #print("Target Array Size: %s" % dsin.shape)




#outpath = "/home/niu4/gliu8/projects/common_data/awi_cm3/"
#datpath = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_control/"



