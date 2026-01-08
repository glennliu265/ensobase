#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Check Nino34 indices from DimaSein Runs

Created on Mon Dec 29 18:11:26 2025

@author: gliu

"""


import sys
import time
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import xarray as xr
import sys
import glob 
import scipy as sp
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
from scipy.io import loadmat
import matplotlib as mpl

import importlib
from tqdm import tqdm


#%%

# local device (currently set to run on Astraeus, customize later)
amvpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/" # amv module
scmpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/"

sys.path.append(amvpath)
sys.path.append(scmpath)

from amv import proc,viz
import scm
import amv.loaders as dl
import cvd_utils as cvd

ensopath = "/Users/gliu/Downloads/02_Research/01_Projects/07_ENSO/03_Scripts/ensobase/"
sys.path.append(ensopath)
import utils as ut

#%%


dpath         = "/Users/gliu/Downloads/02_Research/03_Code/awi_hackathon/scrap/tropical_crop/"
expnames      = ["N39","N43","N44"]
expnames_long = ["N39 (1950-1979)","N43 (2015-2100)", "N44 (1980-2014)"]

ncname        = "%sDimaSein_%s_sst_anom.nc" # %(dpath,expname,)


ninoid_name        = 'nino34'#'nino34' # 

bbox_nino34   = [-170+360,-120+360,-5,5]
bbox_nino3    = [-150+360, -90+360 , -5, 5]  # Nino 3 Box: For SST, <tau_x>

if ninoid_name == "nino34":
    bbox = bbox_nino34
elif ninoid_name == 'nino3':
    bbox = bbox_nino3
    
    
#%%

dsall   = []
ninoids = []
for ex in expnames:
    ncname_in  = ncname % (dpath,ex,)
    ds = xr.open_dataset(ncname_in).sst
    ds = proc.sel_region_xr(ds,bbox).load()
    ninoids.append(proc.area_avg_cosweight(ds))
    dsall.append(ds)
    
    
#%% 

fig,ax = plt.subplots(layout='constrained',figsize=(12.5,4))

for ex in range(len(expnames)):
    plotvar = ninoids[ex]
    ax.plot(plotvar.time_counter,plotvar,label=expnames[ex],lw=1.5)
    
ax.legend()
ax.set_xlabel("Time")
ax.set_ylabel(r"$Ni\tilde{n}o$3.4 [$\degree C$]")
#ax = viz.add_axlines(ax)


#%% Look at Monthly Variance

monstd = [ds.groupby('time_counter.month').std('time_counter') for ds in ninoids]

#%%

x           = np.arange(12)  # the label locations
width       = 0.25  # the width of the bars
multiplier  = 0

mons3 = proc.get_monstr()

fig,ax = plt.subplots(1,1,layout='constrained')

for im in range(12):
    offset = width * multiplier
    
    
plotex = [0,2,1]

for ex in range(3):
    offset  = width * ex
    
    
    
    if ex == 0:
        color_in = 'tab:blue'
    elif ex == 1:
        color_in = "tab:green"
    elif ex == 2:
        color_in = "orange"
    
    measurement = monstd[plotex[ex]].data
    
    rects  = ax.bar(x + offset, measurement, width, label=expnames_long[plotex[ex]],
                    color=color_in,linewidth=0.5)
    
    ax.set_xticks(x,mons3)
    #ax.bar_label(rects, padding=3,fmt="%.02f")
    multiplier += 1

ax.set_ylabel("SST S.D. [$\degree C$]")
ax.set_xlabel("Month")

ax.set_ylim([0,2.5])
    
    

    
ax.legend()



