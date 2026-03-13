#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


Calculate Regression of selected response variable (sst, lcc, pr) to a selected ENSO mode
for a given experiment.


Copied from ``



Created on Mon Mar  9 10:55:05 2026

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

from tqdm import tqdm

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% Additional Functions


# Copied from calculate-Leadlag_ENSO_EOF_toa_fluxes.py
def load_ccf_radiation(expname,flxname,datpath=None,seasonal=False):
    
    if datpath is None:
        datpath = "/home/niu4/gliu8/projects/ccfs/radiative_components/regrid_1x1/"
    
    ccfs = ["sst","eis","Tadv","r700","w700","ws10"]
    
    dscomps = []
    for cc in range(6):
        
        ccf_var  = ccfs[cc]
        if seasonal:
            ncname = "%s%s/%s_%s_component_seasonal.nc" % (datpath,expname,flxname,ccf_var)
        else:
            ncname   = "%s%s/%s_%s_component.nc" % (datpath,expname,flxname,ccf_var)
        dsccf    = xr.open_dataset(ncname)[ccf_var].load()

        dscomps.append(dsccf)
    return dscomps

def load_enso_eof(expname,datpath=None,apply_smoothing=True,sep_mon=False):
    # Load EOF-based ENSO computed via `calc_EOF_ENSO.py`
    if datpath is None:
        datpath = "/home/niu4/gliu8/projects/scrap/nino34/"
    ninonc        = "%s%s_enso_eof_rotated.nc" % (datpath,expname)
    if sep_mon:
        ninonc = proc.addstrtoext(ninonc,"_sepmon",adjust=-1)
    dsnino        = xr.open_dataset(ninonc).load()
    
    # Apply 1-2-1 filter
    if apply_smoothing:
        filter_coeffs = [0.25,0.5,0.25]
        ep     = np.convolve(dsnino.ep,filter_coeffs,mode='same')
        cp     = np.convolve(dsnino.cp,filter_coeffs,mode='same')
        tcoord = dict(time=dsnino.time)
        ep     = xr.DataArray(ep,coords=tcoord,dims=tcoord,name='ep')
        cp     = xr.DataArray(cp,coords=tcoord,dims=tcoord,name='cp')
    else:
        ep     = dsnino.ep
        cp     = dsnino.cp
    return ep,cp

def generate_periods(ds,winlen):
    
    tstart = ds.time[0].dt.year.data.item()
    tend   = ds.time[-1].dt.year.data.item()

    nperiods = tend-tstart-winlen+1
    #nperiods
    tranges = []
    subsets = []
    for ii in range(nperiods):
        trange = ["%i-01-01" % (tstart+ii), "%i-01-01" % (tstart+winlen+ii)]
        subset = ds.sel(time=slice(trange[0],trange[1]))
        tranges.append(trange)
        subsets.append(subset)
    return subsets,tranges

def match_latlon(dstarget,dsref):
    dstarget['lon'] = dsref['lon']
    dstarget['lat'] = dsref['lat']
    return dstarget

def load_input(vname,expname,anom=False,regrid=True):
    inpath = "/home/niu4/gliu8/projects/ccfs/input_data/"
    if regrid:
        inpath += "regrid_1x1"
    if anom:
        procname = "anom_detrend1"
    else:
        procname = "raw"
    ncname = "%s/%s/%s/%s.nc" % (inpath,expname,procname,vname)
    return xr.open_dataset(ncname)[vname]

def sliding_regr_nino_flx(flxwindows,ninoid,leadlags):
    nwin = len(flxwindows)
    
    # ----------
    regr_bywin = []
    for nw in tqdm(range(nwin)):
        
        flx_segment        = flxwindows[nw]
        #ninoin             = ninos_in[nn][ex]
        flx_segment,ninoin = proc.match_time_month(flx_segment,ninoid,verbose=False)
        llout              = ut.calc_leadlag_regression_2d(ninoin,flx_segment,leadlags)
    
        regr_bywin.append(llout)
    # Also Merge into DataArray, concatenating by start date
    ds_winreg     = xr.concat(regr_bywin,dim="trange")
    
    # Add starting date for sliding window
    coords_trange        = dict(trange=np.arange(nwin),pos=['start','end'])
    da_trange            = xr.DataArray(tranges,dims=coords_trange,coords=coords_trange,name='tranges')
    ds_winreg['tranges'] = da_trange
    
    # Add center of starting date
    center_dates=[get_center_date(tt) for tt in tranges]
    coords_tcenter = dict(trange=np.arange(nwin))
    da_tcenter=xr.DataArray(center_dates,dims=coords_tcenter,coords=coords_tcenter,name='tcenters')
    ds_winreg['tcenter'] = da_tcenter
    
    
    return ds_winreg

def preprocess_byperiod(dswins,verbose=False):
    nwin    = len(dswins)
    dsanoms = []
    for nw in range(nwin):
        dsin   = dswins[nw].squeeze()
        dsanom = proc.xrdeseason(dsin,verbose=verbose)
        dsanom = proc.xrdetrend_nd(dsanom,1,verbose=verbose)
        dsanoms.append(dsanom)
    return dsanoms

def get_center_date(trange):
    # Gets Center date given YYYY-01-01 string for [ystart,yend]
    ystart  = int(trange[0][:4])
    yend    = int(trange[1][:4])
    dt      = yend-ystart
    ymiddle = int(ystart + np.ceil(dt/2))
    center_date = "%04i-01-01" % ymiddle
    return center_date



#%% User Edits

expnames   = [
    "TCo319_ctl1950d","TCo319-DART-ssp585d-gibbs-charn", # Target Experiment
]

# Selections for Window Calculations
viz_ccf           = False            # Set to True to visualize flux by CCFs
viz_totalflx      = True            # Visualize the total TOA fluxes (flxnames)
leadlags          = np.arange(-9,10) # Leadlags to compute
winlen            = 30               # Sliding Window Length in Years
outpath           = "/home/niu4/gliu8/projects/ccfs/metrics/regrid_1x1/scrap/sliding_leadlag/"
deseason_byperiod = True # Set to true to deseason each period individually


#metrics_path    = "/home/niu4/gliu8/projects/ccfs/metrics/regrid_1x1/scrap/" # Input Data?

# Get Lengths and set other variables
nino_names   = ["CP","EP"]
flxnames     = ["tscre","ttcre","cre"] # "cre",
ccf_vars     = ["sst","eis","Tadv","r700","w700","ws10"]

nexp         = len(expnames)
nvar         = len(flxnames)
nccfs        = len(ccf_vars)

if deseason_byperiod:
    print("Using raw variables and performing detrend/deseason by period")
    # Using raw variables and performing detrend/deseason by period
    outpath = outpath + "deseason_byperiod/"
    proc.makedir(outpath)


# Additional Loading for visualization
proj    = ccrs.PlateCarree()
figpath = "/home/niu4/gliu8/figures/bydate/2026-03-09/"
proc.makedir(figpath)

#%% Also Load Land Mask

dsmask=ut.load_land_mask_awi("TCo319",regrid=True,)
dsmask.plot()

for ex in range(len(expnames)):
    #%% Load ENSO
    
    ep_indices = []
    cp_indices = []
    for expname in expnames:
        if expname == "CERES_EBAF_ERA5_2001_2024":
            expname_in = "ERA5_1979_2024"
            eraflag    = True # Need to Crop to Dataset Timeperiod Later
        else:
            expname_in = expname
            eraflag    = False
        ep,cp = load_enso_eof(expname_in)
    
        ep_indices.append(ep)
        cp_indices.append(cp)
    
    ninos_in    = [cp_indices,ep_indices]
    
    #%% Load CCFs by Flux
    
    if viz_ccf: # Load CCFs
        ccf_byexp = [] # [exp][flx][ccf]
        for ex in range(nexp):
            expname = expnames[ex]
            
            byflx = []
            for vv in tqdm(range(nvar)):
                flxname = flxnames[vv]
                
                ccfrads = load_ccf_radiation(expname,flxname,seasonal=False)
                byflx.append(ccfrads)
        
            ccf_byexp.append(byflx)
        
    #%% Load total fluxes
    
    # Load total TOA fluxes for both experiments
    # Copied from `awi_cm3_toa_leadlag_analysis_area_avg
    if viz_totalflx:
        flx_byexp = [] # [exp][flx]
        for ex,expname in tqdm(enumerate(expnames)):
            dsflxs = []
            for flxname in flxnames:
                if deseason_byperiod:
                    load_anom=False
                else:
                    load_anom=True
                dsvar = load_input(flxname,expname,anom=load_anom,regrid=True).load()
                dsvar = ut.standardize_names(dsvar)
                if "TCo" in expname:
                    dsvar = ut.varcheck(dsvar,flxname,expname)
                dsvar = ut.remove_duplicate_times(dsvar)
                
                dsvar,_ = proc.match_time_month(dsvar,ep_indices[ex])
                dsflxs.append(dsvar)
        
                
            #dsflxs = xr.merge(dsflxs)
            flx_byexp.append(dsflxs)
    
    #%% Data Has Been Loaded past this point
    
    
    
            
    
    #%% Loop and Calculate the values
    
    ex              = 0 # Set to Zero for the non-reference, can make this into an outer experiment loop eventually...
    
    for vv in tqdm(range(nvar)): # Loop by Flux
    
        
        # Start by Visualizing Total Fluxes
        if viz_totalflx:
            st = time.time()
            print("Calc Regression by Total Fluxes for %s" % flxnames[vv])
            
            # Get Flux and set output file
            input_flux = flx_byexp[ex][vv].squeeze()
    
            # Generate Periods
            flxwindows,tranges = generate_periods(input_flux,winlen)
            
            # Detrend by Period
            if deseason_byperiod:
                st2 = time.time()
                flxwindows = preprocess_byperiod(flxwindows)
                print("\tCompleted period-wise detrend in %.2fs" % (time.time()-st2))
            
            for nn in range(2): # Loop by ENSO Mode
            
            # Get ENSO index
                ninoin             = ninos_in[nn][ex]
                outname            = "%sSliding_LeadLag_ByPeriod_%s_%sENSO_%s_winlen%02i_lagmax%02i.nc" % (outpath,expnames[ex],nino_names[nn],flxnames[vv],winlen,leadlags[-1])
                
                if proc.checkfile(outname):
                    continue
                
                # Perform Regression
                ds_out = sliding_regr_nino_flx(flxwindows,ninoin,leadlags)
                ds_out.to_netcdf(outname)
                
                print("\tCompleted %s in %.2fs" % (flxnames[vv],time.time()-st))
            
        
        if viz_ccf:
            print("Calc Regression by by CCFs for %s" % flxnames[vv])
            
            for cc in tqdm(range(nccfs)):
                st = time.time()
                print("\tCCF: %s" % ccf_vars[cc])
                
                # Get Flux and set output file
                input_flux = ccf_byexp[ex][vv][cc].squeeze()
                #outname    = "%sSliding_LeadLag_ByPeriod_%s_%sENSO_%s_CCF%s_winlen%02i_lagmax%02i.nc" % (outpath,expnames[ex],nino_names[nn],flxnames[vv],ccf_vars[cc],winlen,leadlags[-1])
                
                # Generate Periods
                flxwindows,tranges = generate_periods(input_flux,winlen)
                
                # Detrend by Period
                if deseason_byperiod:
                    st2 = time.time()
                    flxwindows = preprocess_byperiod(flxwindows)
                    print("\tCompleted period-wise detrend in %.2fs" % (time.time()-st2))
                
                # Loop for ENSO mode
                for nn in range(2):
                    ninoin  = ninos_in[nn][ex]
                    outname = "%sSliding_LeadLag_ByPeriod_%s_%sENSO_%s_CCF%s_winlen%02i_lagmax%02i.nc" % (outpath,expnames[ex],nino_names[nn],flxnames[vv],ccf_vars[cc],winlen,leadlags[-1])
                    if proc.checkfile(outname):
                        continue
                    
                    ds_out = sliding_regr_nino_flx(flxwindows,ninoin,leadlags)
                    ds_out.to_netcdf(outname)
                
                print("\tCompleted %s CCF %s in %.2fs" % (flxnames[vv],ccf_vars[cc],time.time()-st))
                    
    


