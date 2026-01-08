#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Do Significance Testing by Subsampling a dataset and repeating the calculation N Times


General Setup

<Part 1: Set-up>
- Load Dataset to Sample For
- Indicate Sample Size


<Part 2: Perform Calculations>



<Part 3: Save Output Statistics...>




Created on Fri Dec 19 11:20:30 2025

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

#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% User Edits

inpath = ""
innc   = ""

# MonteCarlo Options
mciter        = 1000
sample_length = 10*12 # In Months

# Other Toggles
"""
preprocess_method
    None            : No Preprocessing Applied 
    before_chunking : Preprocess whole timeseries before chunking
    by_chunk        : Preprocess each chunk
"""
preprocess_method = None

#%% Define Functions

def preprocess_ds(ds):
    return ds_preprocessed


#%% Set up Functions to do calculations


