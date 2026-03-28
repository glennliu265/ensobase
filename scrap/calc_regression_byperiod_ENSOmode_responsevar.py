#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


Calculate Regression of selected response variable (sst, lcc, pr) to a selected ENSO mode
for a given experiment.


Copied from `calc_regression-byperiod_ensomode`



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

# def load_enso_eof(expname,datpath=None,apply_smoothing=True,sep_mon=False,by_period=False,winlen=None):
#     # Load EOF-based ENSO computed via `calc_EOF_ENSO.py`
#     # by_period loads for sliding period [winlen], currently does not support [sep_mon] `calc_EOF_ENSO_sliding.py`
    
#     if by_period:
#         print("Loading ENSO computed for sliding windows of %i-years" % winlen)
#         if datpath is None:
#             datpath       = "/home/niu4/gliu8/projects/ccfs/enso_eof/"
#         ninonc        = "%s%s_ENSO_EOF_slidingwinlen%02i.nc" % (datpath,expname,winlen)

#     else: # Load for calculation over whole timeseries
#         if datpath is None:
#             datpath = "/home/niu4/gliu8/projects/scrap/nino34/"
#         ninonc        = "%s%s_enso_eof_rotated.nc" % (datpath,expname)
#         if sep_mon:
#             ninonc    = proc.addstrtoext(ninonc,"_sepmon",adjust=-1)
#     dsnino        = xr.open_dataset(ninonc).load()
    
#     # Apply 1-2-1 filter
#     if apply_smoothing:
#         if by_period:
#             filter_coeffs = [0.25,0.5,0.25]
#             smooth121 = lambda timeseries: np.convolve(timeseries,filter_coeffs,mode='same')
#             ep = xr.apply_ufunc(smooth121, dsnino.ep, input_core_dims=[["timeindex"]],output_core_dims=[["timeindex"]],vectorize=True)
#             cp = xr.apply_ufunc(smooth121, dsnino.cp, input_core_dims=[["timeindex"]],output_core_dims=[["timeindex"]],vectorize=True)
#             ninotimes = [make_ninotime(trange,dsnino.timeindex) for trange in dsnino.trange] 
#         else:
#             ep     = np.convolve(dsnino.ep,filter_coeffs,mode='same')
#             cp     = np.convolve(dsnino.cp,filter_coeffs,mode='same')
#             tcoord = dict(time=dsnino.time)
#             ep     = xr.DataArray(ep,coords=tcoord,dims=tcoord,name='ep')
#             cp     = xr.DataArray(cp,coords=tcoord,dims=tcoord,name='cp')
        
#     else:
#         ep     = dsnino.ep
#         cp     = dsnino.cp
#     if by_period:
#         return ep,cp,ninotimes
#     return ep,cp

# def generate_periods(ds,winlen):
    
#     tstart = ds.time[0].dt.year.data.item()
#     tend   = ds.time[-1].dt.year.data.item()

#     nperiods = tend-tstart-winlen+1
#     #nperiods
#     tranges = []
#     subsets = []
#     for ii in range(nperiods):
#         trange = ["%i-01-01" % (tstart+ii), "%i-01-01" % (tstart+winlen+ii)]
#         subset = ds.sel(time=slice(trange[0],trange[1]))
#         tranges.append(trange)
#         subsets.append(subset)
#     return subsets,tranges

# def make_ninotime(trange,timeindex):
#     ntime = len(timeindex)
#     times = xr.date_range(start=trange[0].data.item(),end=trange[1].data.item(),periods=ntime)
#     return times # Note this does not return the exact middle of the month
    
def match_latlon(dstarget,dsref):
    dstarget['lon'] = dsref['lon']
    dstarget['lat'] = dsref['lat']
    return dstarget

def load_input_regrid(vname,expname,anom=False,regrid=True):
    inpath = "/home/niu4/gliu8/projects/scrap/regrid_1x1/"
    if anom: # Load Anomalies 
        inpath = "/home/niu4/gliu8/projects/scrap/regrid_1x1/anom_detrend1/"
    ncname = "%s%s_%s_regrid1x1.nc" % (inpath,expname,vname)
    return xr.open_dataset(ncname)[vname]

def sliding_regr_nino_flx(flxwindows,ninoid,leadlags,nino_byperiod=False,nino_times=None):
    nwin = len(flxwindows)
    
    # ----------
    regr_bywin = []
    tranges    = []
    for nw in tqdm(range(nwin)):
        
        flx_segment        = flxwindows[nw]
        if nino_byperiod:
            nino_segment = ninoid[nw]
            nino_segment = nino_segment.rename(dict(timeindex='time'))
            nino_segment['time'] = nino_times[nw]

            if nw == 0 and (nwin != len(ninoid)):
                print("Warning! Length of variable (%i) != Length of Nino Segments (%i)" % (nw,len(ninoid)))
        else:
            nino_segment = ninoid
        
        flx_segment,ninoin = proc.match_time_month(flx_segment,nino_segment,verbose=False)
        llout              = ut.calc_leadlag_regression_2d(ninoin,flx_segment,leadlags)
        
        tstart = str(flx_segment.time[0].data)[:10]
        tend   = str(flx_segment.time[-1].data)[:10]
        tranges.append([tstart,tend])
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

# def preprocess_byperiod(dswins,verbose=False):
#     nwin    = len(dswins)
#     dsanoms = []
#     for nw in range(nwin):
#         dsin   = dswins[nw].squeeze()
#         dsanom = proc.xrdeseason(dsin,verbose=verbose)
#         dsanom = proc.xrdetrend_nd(dsanom,1,verbose=verbose)
#         dsanoms.append(dsanom)
#     return dsanoms

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
    "EC-Earth3",
    "E3SM-1-0",
    "E3SM-1-1-ECA",
    #"TCo319_ssp585",#"TCo319_ctl1950d","TCo319-DART-ssp585d-gibbs-charn", # Target Experiment
]

# Selections for Window Calculations
leadlags          = np.arange(-9,10) # Leadlags to compute
winlen            = 30               # Sliding Window Length in Years
outpath           = "/home/niu4/gliu8/projects/ccfs/metrics/regrid_1x1/scrap/sliding_leadlag/"
deseason_byperiod = True # Set to true to deseason each period individually
enso_eof_byperiod = True # Set to True to load ENSO-EOF conducted over each period...
overwrite         = True # Set to True to overwrite existing files

#metrics_path    = "/home/niu4/gliu8/projects/ccfs/metrics/regrid_1x1/scrap/" # Input Data?

# Get Lengths and set other variables
nino_names   = ["CP","EP"]
varnames     = ["tscre","sst"]#"lcc","sst","pr","ttcre","cre"]#"tscre",]# ["lcc","sst"]#"pr","lcc"] # "cre",

nexp         = len(expnames)
nvar         = len(varnames)

if deseason_byperiod:
    print("Using raw variables and performing detrend/deseason by period")
    # Using raw variables and performing detrend/deseason by period
    outpath = outpath + "deseason_byperiod/"
    proc.makedir(outpath)
    
if enso_eof_byperiod:
    print("Using raw variables and performing detrend/deseason by period AND ENSO EOF by period")
    
    outpath = outpath + "deseason_enso_eof_byperiod/"
    proc.makedir(outpath)


# Additional Loading for visualization
proj    = ccrs.PlateCarree()
figpath = "/home/niu4/gliu8/figures/bydate/2026-03-09/"
proc.makedir(figpath)

#%% Also Load Land Mask

dsmask = ut.load_land_mask_awi("TCo319",regrid=True,)
dsmask.plot()

for ex in range(len(expnames)): # Loop by Experiment
    
    expname = expnames[ex]
    
    #%% Load ENSO
    if expname == "CERES_EBAF_ERA5_2001_2024":
        expname_in = "ERA5_1979_2024"
        eraflag    = True # Need to Crop to Dataset Timeperiod Later
    else:
        expname_in = expname
        eraflag    = False
    
    if enso_eof_byperiod:
        ep_index,cp_index,ninotimes = ut.load_enso_eof(expname_in,by_period=enso_eof_byperiod,winlen=winlen)
    else:
        ep_index,cp_index= ut.load_enso_eof(expname_in,by_period=enso_eof_byperiod,winlen=winlen)
    ninos_in          = [cp_index,ep_index]
    

        
    #%% Load response variable
    for vv,varname in tqdm(enumerate(varnames)): # Loop by variable
        
        st = time.time()
        print("Calc Regression with %s for %s" % (varnames[vv],expname))
        
        if deseason_byperiod:
            load_anom = False
        else:
            load_anom = True
        
        dsvar = load_input_regrid(varname,expname,anom=load_anom,regrid=True).load()
        dsvar = ut.standardize_names(dsvar)
        if "TCo" in expname:
            dsvar = ut.varcheck(dsvar,varname,expname)
        dsvar   = ut.remove_duplicate_times(dsvar)
        
        if enso_eof_byperiod is False:
            # Align to whoe period
            dsvar,_ = proc.match_time_month(dsvar,ep_index)
        
        
        # Generate Periods
        varwindows,tranges = ut.generate_periods(dsvar,winlen)
        
        # Detrend by Period
        if deseason_byperiod:
            st2 = time.time()
            varwindows = ut.preprocess_byperiod(varwindows)
            print("\tCompleted period-wise detrend in %.2fs" % (time.time()-st2))
        
        # Get ENSO index
        for nn in range(2):

            ninoin             = ninos_in[nn]
            
            
            outname            = "%sSliding_LeadLag_ByPeriod_%s_%sENSO_%s_winlen%02i_lagmax%02i.nc" % (outpath,expname,nino_names[nn],varname,winlen,leadlags[-1])
            
            if proc.checkfile(outname) and (overwrite is False):
                continue
            
            # Perform Regression
            ds_out = sliding_regr_nino_flx(varwindows,ninoin,leadlags,nino_byperiod=enso_eof_byperiod,nino_times=ninotimes)
            ds_out.to_netcdf(outname)
            
            print("\tCompleted %s in %.2fs" % (varname,time.time()-st))
            
        
        # Close unnecesary variables
        dsvar.close()
        
                    
    


