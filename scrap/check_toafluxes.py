#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualize Mean States of TOA Fluxes in ERA5 for compariso

Created on Wed Nov 26 15:08:31 2025

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
import importlib

from sklearn.linear_model import LinearRegression
import sklearn

#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%%
datpath = '/home/niu4/gliu8/projects/common_data/ERA5/regrid_1x1/'

dstsr   = xr.open_dataset(datpath + 'tsr_1979_2024.nc').load()
dtday   = 3600*24
tsr   = dstsr.tsr/dtday

tsr.mean('valid_time').plot(vmin=0,vmax=160,cmap='cmo.solar'),plt.show()

tsr.mean('valid_time').plot(vmin=50,vmax=350,cmap='cmo.solar'),plt.show()


#%% ttr
dsttr   = xr.open_dataset(datpath + 'ttr_1979_2024.nc').load()
ttr   = dsttr.ttr/dtday

ttr.mean('valid_time').plot(cmap='cmo.dense_r',vmin=-300,vmax=-140),plt.show()


#%% All Sky

vname   = "allsky"
ds      = xr.open_dataset(datpath + '%s_1979_2024.nc' % vname).load()
dsvar   = ds[vname]/dtday

dsvar.mean('valid_time').plot(cmap='cmo.balance',vmin=-100,vmax=100),plt.show()

#%% TTRC

vname   = "ttrc"
ds      = xr.open_dataset(datpath + '%s_1979_2024.nc' % vname).load()
dsvar   = ds[vname]/dtday

dsvar.mean('valid_time').plot(cmap='cmo.dense_r',vmin=-300,vmax=-140),plt.show()

#%% Tsrc

vname   = "tsrc"
ds      = xr.open_dataset(datpath + '%s_1979_2024.nc' % vname).load()
dsvar   = ds[vname]/dtday

dsvar.mean('valid_time').plot(vmin=50,vmax=350,cmap='cmo.solar'),plt.show()


#%% Clearsky

vname   = "clearsky"
ds      = xr.open_dataset(datpath + '%s_1979_2024.nc' % vname).load()
dsvar   = ds[vname]/dtday

dsvar.mean('valid_time').plot(cmap='cmo.balance',vmin=-100,vmax=100),plt.show()

#%% CRE

vname   = "cre"
ds      = xr.open_dataset(datpath + '%s_1979_2024.nc' % vname).load()
dsvar   = ds[vname]/dtday

dsvar.mean('valid_time').plot(cmap='cmo.balance',vmin=-100,vmax=100),plt.show()

#%% Tscre (apparent it has not been renamed)

vname   = "tscre"
ds      = xr.open_dataset(datpath + '%s_1979_2024.nc' % vname).load()
dsvar   = ds['tsr']/dtday

dsvar.mean('valid_time').plot(vmin=-100,vmax=100,cmap='cmo.balance'),plt.show()

#%% ttcre 

vname   = "ttcre"
ds      = xr.open_dataset(datpath + '%s_1979_2024.nc' % vname).load()
dsvar   = ds[vname]/dtday


dsvar.mean('valid_time').plot(vmin=-100,vmax=100,cmap='cmo.balance'),plt.show()