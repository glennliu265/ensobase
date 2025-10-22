#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Global Mean of selected variable for AWI-CM3

Note: It appears that the cdo fldmean is identical to taking the cosine 
weighted mean in python (but much faster). Thus I will use fldmean to do so

example:
    
    cdo fldmean /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/TCo319_ssp585_sst_1m_2015-2114_1x1regrid.nc /home/niu4/gliu8/projects/scrap/global_mean/TCo319_ssp585_sst_global_mean.nc
    
    

Created on Wed Oct 15 10:55:03 2025

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


#%% Niu Paths for Custom Modules

scmpath = "/home/niu4/gliu8/scripts/commons/stochmod/model"
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

sys.path.append(scmpath)
import scm

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% Other User Edits

gmeanpath = "/home/niu4/gliu8/projects/scrap/global_mean/"

#%% Part (1) Check CDO fldmean function

ncname = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/TCo319_ssp585_sst_1m_2015-2114_1x1regrid.nc"
ds     = xr.open_dataset(ncname).load()


#st = time.time()
sstmean_unw = ds.sst.mean(('lat','lon'))

sstmean_cw  = proc.area_avg_cosweight(ds.sst)

# Load CDO fldmean
nc2 = gmeanpath + "TCo319_ssp585_sst_global_mean.nc"
dscdo = xr.open_dataset(nc2).load()

#%% Plot the three to check

fig,ax = plt.subplots(1,1)

plotvar = sstmean_unw
ax.plot(plotvar.time_counter,plotvar,label="Unweighted Global Mean")

plotvar = sstmean_cw
ax.plot(plotvar.time_counter,plotvar,label="Global Mean (Cosweight)",ls='solid')

plotvar = dscdo.sst.squeeze()
ax.plot(plotvar.time_counter,plotvar,label="CDO FldMean",ls='dashed')

ax.legend()
plt.show()



