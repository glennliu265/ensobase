#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Make Regridding File
Based on regrid_re1x1.nc from Ray's cdo scripts

Created on Mon Oct 27 11:19:19 2025

@author: gliu

"""


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
from scipy.io import loadmat


#%% Load custom modules

# local device (currently set to run on Astraeus, customize later)
amvpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/" # amv module
scmpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/"

sys.path.append(amvpath)
sys.path.append(scmpath)

from amv import proc,viz
import scm
import amv.loaders as dl
import cvd_utils as cvd



#%% 

regrid_path     = "~/"
regrid_example  = "regrid_re1x1.nc"
#regrid_example = "regrid_CESM1CAM.nc"
ds_example      = xr.open_dataset(regrid_path + regrid_example)

#%% Check Regridded 

regridpath      = "/home/niu4/gliu8/projects/scrap/regrid_1x1/"
ncname          = regridpath + "TCo2559-DART-1950C_ttr_regrid1x1.nc"
ds              = xr.open_dataset(ncname)

#%% 