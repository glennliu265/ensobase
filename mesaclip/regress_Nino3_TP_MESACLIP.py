#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Perform regression of Nino3 Indices onto Regridded MESACLIP Variables
based on segment from `Nino3_Analysis_MESACLIP.ipynb
Created on Thu Jul 23 11:49:13 2026

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

#%% User Edits

# Figure Path
figpath   = "/home/niu4/gliu8/figures/bydate/2026-07-29-CGC/"
proc.makedir(figpath)

# Output Path
outpath      = "/home/niu4/gliu8/projects/mesaclip/indices/"    # Path to Nino3 Indices
outpath_regr = "/home/niu4/gliu8/projects/mesaclip/regressions" # Regression Output
expnames  = ["lores","hires"]
ensnums   = [22,10] # Max # of Ensemble Members
vname     = "TS"


bbsel = [120,290,-30,30] # Tropical Pacific Subset

# Plotting Variables
bbox_cep        = [150      , -130+360 , -5, 5]   # Central Equatorial Pacific, for [tau_x]
bbox_nino3      = [-150+360, -90+360 , -5, 5]  # Nino 3 Box: For SST, <tau_x>
proj            = ccrs.PlateCarree()

#%% Load Indices

# (1) Load Nino3
nino3_byexp = []
for ex in range(2):
    expname = expnames[ex]
    nens    = ensnums[ex]
    ncname  = "%s/CESM1_%s_TS_nino3_DetrendEnsMean_nens%02i.nc" % (outpath,expname,nens)
    ds      = xr.open_dataset(ncname).load()
    nino3_byexp.append(ds)
    
#%%
    
for ex in range(2): # Loop by Experiment
    
    expname  = expnames[ex]
    nens     = ensnums[ex]
    
    for e in tqdm(range(nens)):
        ensnum   = e+1
        
        # Set Data Path and Search String
        datpath  = "/home/niu4/gliu8/share/CESM1/MESACLIP/RDA/%s/BHIST/month_1/regrid_1x1/ens%03i/" % (expname,ensnum)
        ncsearch = datpath + "%s_*_regrid1x1.nc" % vname
        nclist   = glob.glob(ncsearch)
        nclist.sort()
        print("Found %i files for %s (ens=%03i)" % (len(nclist),expname,ensnum))
        
        # Subset to Tropical Pacific and Load
        dsall = xr.open_mfdataset(nclist,concat_dim='time',combine='nested')
        dsreg = proc.sel_region_xr(dsall,bbsel).load()
        
        # Deseason and Detrend (Lienar)
        dsreg_anom    = proc.xrdeseason(dsreg[vname])
        dsreg_anom_dt = proc.xrdetrend(dsreg_anom)
        
        # Do Instantaneous Regression
        ninoin = nino3_byexp[0].isel(ens=0).TS
        ninoin,dsanomin = proc.match_time_month(ninoin,dsreg_anom_dt)
        fitout = proc.pointwise_polyfit(ninoin,dsanomin,1,ttest=True)
        
        # Save the Output
        outname      = "%s/%s_Nino3_Regression_TPac_%s_regrid1x1_ens%03i.nc" % (outpath_regr,vname,expnames[ex],ensnum)
        fitout       = fitout.drop_vars(['model','residual'])
        fitout.to_netcdf(outname)
        print(outname)
        
                
        # # Make a Plot =========================================================
        # cints           = np.arange(-1.8,2.0,0.2) #* 1e-2
        # fig,ax          = viz.init_tp_map(1,1)
        
        # plotvar         = fitout.coefficients_by_degree.isel(coeff=1) * -1 * 1e2 # Additional Multipliers for Wind Stress
        # pcm             = ax.contourf(plotvar.lon,plotvar.lat,plotvar,transform=proj,
        #                         levels=cints,cmap='cmo.balance')
        # cl              = ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,
        #                         levels=cints,colors="k",linewidths=0.75)
        
        # # Plot the Boxes
        # viz.plot_box(bbox_cep,color="yellow",ax=ax,proj=proj,linewidth=3,linestyle='dashed')
        # viz.plot_box(bbox_nino3,color="magenta",ax=ax,proj=proj,linewidth=3,linestyle='dotted')
        
        # # Plot Significance
        # plotmask = fitout.significance
        # viz.plot_mask(plotmask.lon,plotmask.lat,plotmask,reverse=False,ax=ax,proj=proj,geoaxes=True,markersize=0.8,color='gray')
        
        # # Set Labels
        # ax.clabel(cl,levels=cints)
        # ax.set_title("Zonal Wind Stress Regression to Niño3, Ens %03i (%s)" % (ensnum,expnames_fancy[ex]))
        # cb      = viz.vcbar(pcm,ax=ax)
        # cb.set_label("$N m^{-2} x 10^{-2}$")
        
        # outname = "%s/TAUX_Nino3_Regression_TPac_%s_ens%0i3.png" % (figpath,expnames[ex],ensnum)
        # print(outname)
        # plt.savefig(outname)
        
        
        
        
        