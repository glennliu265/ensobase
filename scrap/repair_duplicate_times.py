#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 16:46:18 2026

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

from tqdm.notebook import tqdm

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut


#%% TRy TSR

dir    = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_control/"
dstsrc = xr.open_dataset(dir+"ctl1950d_atm_remapped_1m_tsrc_1850-2134.nc").load()
dstsr  = xr.open_dataset(dir+"ctl1950d_atm_remapped_1m_tsr_1850-2134.nc").load() #THis was the iossue case, has 3780 timesteps instead of 3480

#dstsrc,dstsr


#%%

timedict  = dict(time_counter='time')
timedictr = dict(time='time_counter')
dstsr_fix = ut.remove_duplicate_times(dstsr.rename(timedict))

#%%
dstsr_fix = dstsr_fix.rename(timedictr)

outname   = "/home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_tsr.nc"
dstsr_fix.to_netcdf(outname)



