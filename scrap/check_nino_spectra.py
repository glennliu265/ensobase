#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Visualize power spectra for Nino3.4 Indices (or Index of Choice)



Created on Tue Oct 14 14:02:56 2025

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


#%% Niu Paths

scmpath = "/home/niu4/gliu8/scripts/commons/stochmod/model"
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

sys.path.append(scmpath)
import scm

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut


#%% Commons Variables


#%% Shared variable names and experiment
# Copied from visualize composites

#datpath         = "/Users/gliu/Downloads/02_Research/01_Projects/07_ENSO/01_Data/TP_Crop/composites/"


# Simulation Names -----
expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"]
expcols         = ["cornflowerblue","red","cornflowerblue","red","cornflowerblue"]
expls           = ["solid","solid",'dashed','dashed','dotted','dotted']
emks            = ["v","v","o","o","^","^"]

# Initial Variable Analysis -----
vnames          = ["sst","ssr","str","tx_sur","D20","Dmaxgrad"]
vunits          = [r"$\degree C$",r"$\frac{W}{m^2}$",r"$\frac{W}{m^2}$",r"$\frac{m}{s^2}$","m","m"]
vnames_long     = ["SST","Surface Shortwave","Surface Longwave","Zonal Wind Stress","Thermocline (20$\degree$ Isotherm)","Thermocline (Max Vertical Gradient)"]
vmaxes          = [2,40,20,0.02,20,20]

# ENSO Names -----
ninoname        = [r"$El$ $Ni\tilde{n}o$",r"$La$ $Ni\tilde{n}a$"]
ninocol         = ["cornflowerblue","firebrick"]
ninoshort       = ['nino','nina']

# Conversion for STR and SSR considering 3h Accumulation -----
conversion      = 1/(3*3600) # 3 h accumulation time...? #1/(24*30*3600)
# https://forum.ecmwf.int/t/surface-radiation-parameters-joule-m-2-to-watt-m-2/1588

# Bounding Boxes from Jin et al. 2020 Eqn. 6.6  -----
bbox_cep        = [150      , -130+360 , -5, 5]   # Central Equatorial Pacific, for [tau_x], 
bbox_nino3      = [-150+360 , -90+360  , -5, 5]  # Nino 3 Box: For SST, <tau_x>
bbox_nino34     = [-170+360 , -120+360 , -5, 5]
bbox_epac       = [-155+360 , -80+360  , -5, 5]  # Eastern Pacific (for h_e calculation)
bbox_wpac       = [120      , -155+360 , -5, 5]  # Western Pacific (for h_w calculation)

bboxes          = [bbox_cep,bbox_nino3,bbox_nino34,bbox_epac,bbox_wpac]
bbnames_long    = ["Central Equatorial Pacific",r"$Ni\tilde{n}o3$",r"$Ni\tilde{n}o3.4$","Tropical Eastern Pacific","Tropical Western Pacific"]
bbnames         = ["CEP","nino3","nino34","EPac","WPac"]

expcols = ["cornflowerblue",'lightcoral',
           "slateblue","firebrick",
           "midnightblue","k"]

#%% Indicate paths

figpath         = "/home/niu4/gliu8/figures/bydate/2025-10-14/"
proc.makedir(figpath)

datpath         = "/home/niu4/gliu8/projects/scrap/TP_crop/"
datpath_anom    = datpath + "anom_detrend2/"



ninopath        = "/home/niu4/gliu8/projects/scrap/nino34/"

nexps           = len(expnames)


#%% Load ENSO Indices

ninoid_name   = "nino34"
unstandardize = True
ds_enso       = []
ensoids       = []
for ex in range(nexps):
    
    # ninonc      = "%s%s_%s.nc" % (ninopath,expnames[ex],ninoid_name)
    # ds          = xr.open_dataset(ninonc).load()
    
    ds = ut.load_ensoid(expnames[ex],ninoid_name,standardize=True)
    
    ds_enso.append(ds)
    ensoids.append(ds * ds['std'].data.item())
    
if unstandardize:
    ensoids = [ds * ds['std'].item() for ds in ensoids]
    
#%%

ntimes   = [ts.time.shape[0] for ts in ensoids]
nsmooths = [1,1,1,1,1]
in_ts    = [ds.data for ds in ensoids]

pct      = 0.10
specout  = scm.quick_spectrum(in_ts,nsmooths,pct,return_dict=True)

#%% Plot the spectra

secday = 3600*24
secyr  = 3600*24*365
secmon = 3600*24*30

fig,ax = plt.subplots(1,1,figsize=(12.5,4.5),constrained_layout=True)

for ex in range(nexps):
    
    plotspec = specout['specs'][ex] #/ secmon #/ dtplot
    plotfreq = specout['freqs'][ex] #* secmon#* dtplot
    
    ax.plot(plotfreq,plotspec,label=expnames_long[ex],c=expcols[ex],lw=2.5,marker="o")
    


ax.axvline([1/(8*secmon)],label="8 Mon")
ax.axvline([1/(12*secmon)],label="1 Year",c=[0.9,0.9,0.9])
ax.axvline([1/(15*secmon)],label="15 Mon",c=[0.7,0.7,0.7])
ax.axvline([1/(3*secyr)],label="3 Year",c=[0.5,0.5,0.5])
ax.axvline([1/(5*secyr)],label="5 Year",c=[0.3,0.3,0.3])
ax.axvline([1/(8*secyr)],label="8 Year",c=[0,0,0])


ax.legend()

ax.set_xlim([1/(20*secyr),(1/(6*secmon))])

plt.show()


    