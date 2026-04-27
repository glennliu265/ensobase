#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copied from Calc_EOF_ENSO.py

Created on Mon Mar 16 20:35:09 2026

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

#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%%


def calc_enso(invar,lon,lat,pcrem,bbox=None,sep_mon=True):
    
    """
    Calculate ENSO for input SST (use longitude 0-360)!
    
    Parameters
    ----------
    invar : ARRAY [time x lat x lon360]
        monthly SST variable longitude is 0-360, lat -90 to 90
    lon : [lon]
        Longitude values
    lat : [lat]
        Latitude values
    pcrem : Int
        Number of PCs to alculate
    bbox : List [lonW,lonE,latS,latN], optional
        Bounding box for calculation.
        Default is [120, 290, -20, 20]. The default is None.

    Returns
    -------
    eofall : ARRAY [lat x lon x month x pc]
        EOF Patterns.
    pcall : ARRAY [time x month x pc]
        Principle components.
    varexpall : [month x pc]
        Variance Explained.
        
    Custom dependencies: check_enso_sign, proc.eof_simple
    """
    if bbox is None:
        bbox = [120, 290, -20, 20]
    
    #  Apply Area Weight
    _,Y = np.meshgrid(lon,lat)
    wgt = np.sqrt(np.cos(np.radians(Y))) # [lat x lon]
    ts = invar * wgt[None,:,:]
    
    # Reshape for ENSO calculations
    ntime,nlat,nlon = ts.shape 
    ts = ts.reshape(ntime,nlat*nlon) # [time x space]
    ts = ts.T #[space x time]
    
    # Remove NaN points
    okdata,knan,okpts = proc.find_nan(ts,1) # Find Non-Nan Points
    oksize = okdata.shape[0]
    
    # Calcuate monthly anomalies
    if sep_mon:
        okdata = okdata.reshape(oksize,int(ntime/12),12) # [space x yr x mon]
        manom  = okdata.mean(1)
        tsanom = okdata - manom[:,None,:]
        nyr    = tsanom.shape[1]
        
        # Preallocate
        eofall = np.zeros((nlat*nlon,12,pcrem)) *np.nan# [space x month x pc]
        pcall  = np.zeros((nyr,12,pcrem)) *np.nan# [year x month x pc]
        varexpall  = np.zeros((12,pcrem)) * np.nan #[month x pc]
    
        # Compute EOF!!
        for m in range(12):
            
            # Perform EOF
            st = time.time()
            eofsok,pcs,varexp= proc.eof_simple(tsanom[:,:,m],pcrem,1)
            #print("Performed EOF in %.2fs"%(time.time()-st))
        
            # Place back into full array
            eofs = np.zeros((nlat*nlon,pcrem)) * np.nan
            eofs[okpts,:] = eofsok   
        
            # Correct ENSO Signs
            eofs,pcs = check_ENSO_sign(eofs,pcs,lon,lat,verbose=True)
            
            
            # Save variables
            eofall[:,m,:] = eofs.copy()
            pcall[:,m,:] = pcs.copy()
            varexpall[m,:] = varexp.copy()
            
            print("Completed month %i in %.2fs"%(m+1,time.time()-st))
    
    
        # Reshape Data
        eofall = eofall.reshape(nlat,nlon,12,pcrem) # [lat x lon x month x pc]
    else: # Just do all months together
        manom  = okdata.mean(1)
        tsanom = okdata - manom[:,None] # Remove Mean
        
        # Perform EOF
        st = time.time()
        eofsok,pcs,varexpall= proc.eof_simple(tsanom,pcrem,1)
        
        # Place back into full array
        eofall = np.zeros((nlat*nlon,pcrem)) * np.nan
        eofall[okpts,:] = eofsok   
        
        # Correct ENSO Signs
        eofall,pcall = check_ENSO_sign(eofall,pcs,lon,lat,verbose=True)
        
        # Reshape Data
        eofall = eofall.reshape(nlat,nlon,pcrem) # [lat x lon x month x pc]
        print("Completed ENSO calculation in %.2fs"%(time.time()-st))
    
    return eofall,pcall,varexpall


def check_ENSO_sign(eofs,pcs,lon,lat,verbose=True):
    """
    checks sign of EOF for ENSO and flips sign by check the sum over
    conventional ENSO boxes (see below for more information)
    
    checks to see if the sum within the region and flips if sum < 0
    
    inputs:
        1) eofs [space (latxlon), PC] EOF spatial pattern from eof_simple
        2) pcs  [time, PC] PC timeseries calculated from eof_simple
        3) lat  [lat] Latitude values
        4) lon  [lon] Longitude values
    
    outputs:
        1) eofs [space, PC]
        2) pcs  [time, PC]
    
    """   
    # Set EOF boxes to check over (west,east,south,north)
    eofbox1 = [190,240,-5,5] #EOF1 sign check Lat/Lon Box (Currently using nino3.4)
    eofbox2 = [190,240,-5,5] #EOF2 (Frankignoul et al. 2011; extend of positive contours offshore S. Am > 0.1)
    eofbox3 = [200,260,5,5] #EOF3 (Frankignoul et al. 2011; narrow band of positive value around equator)
    chkboxes = [eofbox1,eofbox2,eofbox3]
    
    # Find dimensions and separate out space
    nlon = lon.shape[0]
    nlat = lat.shape[0]
    npcs = eofs.shape[1] 
    eofs = eofs.reshape(nlat,nlon,npcs)
    eofs = eofs.transpose(1,0,2) # [lon x lat x npcs]
        
    for n in range(npcs):
        if n >= len(chkboxes):
            print("Not checking EOF %i" %(n+1))
            continue
        chk = proc.sel_region(eofs,lon,lat,chkboxes[n],reg_sum=1)[n]
        
        if chk < 0:
            if verbose:
                print("Flipping EOF %i because sum is %.2f"%(n,chk))
    
            eofs[:,:,n] *= -1
            pcs[:,n] *= -1
    
    eofs = eofs.transpose(1,0,2).reshape(nlat*nlon,npcs) # Switch back to lat x lon x pcs
    
    return eofs,pcs

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

#%%

# Calculate for AWI-CM3
expname   = 'TCo319_ssp585'
mergefile = False
infile    = '/home/niu4/gliu8/projects/scrap/regrid_1x1/TCo319_ssp585_sst_regrid1x1.nc'
outpath   = '/home/niu4/gliu8/projects/ccfs/enso_eof/'
winlen    = 24
savename  = "%s%s_ENSO_EOF_slidingwinlen%02i.nc" % (outpath,expname,winlen)

# Dataset Information
vname    = 'sst'
timename = 'time_counter'
latname  = 'lat'
lonname  = 'lon'

# Calculation Options
tstart = '2015-01-01'
tend   = '2024-12-31'


#%%

# Calculate for E3SM-1-0
#expname   = 'E3SM-1-0'

# Loop for CMIP6 Experiments (where merging is needed)
for expname in ["E3SM-1-1-ECA","EC-Earth3"]:
    
    mergefile = True
    infile    = [
        "/home/niu4/gliu8/projects/scrap/regrid_1x1/%s_historical_sst.nc" % expname, 
        "/home/niu4/gliu8/projects/scrap/regrid_1x1/%s_ssp585_sst.nc" % expname,
        ]
    outpath   = '/home/niu4/gliu8/projects/ccfs/enso_eof/'
    winlen    = 30
    savename  = "%s%s_ENSO_EOF_slidingwinlen%02i.nc" % (outpath,expname,winlen)
    
    # Dataset Information
    vname    = 'sst'
    timename = 'time'
    latname  = 'lat'
    lonname  = 'lon'
    
    # Calculation Options
    tstart = '2014-01-01'
    tend   = '2099-12-31'
    # -----------------
    
    #%% addiitonal settings
    
    
    bbox_takahashi = [120,360-70,-10,10]
    pcrem          = 3
    sep_mon        = False
    
    if sep_mon:
        savename_in = proc.addstrtoext(savename,"_sepmon",adjust=-1)
    else:
        savename_in = savename
    
    #%% Part 1: Load the file
    
    if mergefile is False:
        dsraw  = xr.open_dataset(infile)
        dsreg  = proc.sel_region_xr(dsraw,bbox_takahashi).load()
        sstraw = ut.standardize_names(dsreg.sst)
        sstraw = sstraw.transpose('time','lat','lon')
    else:
        
        dsraws = []
        for ii in range(len(infile)):
            
            dsraw  = xr.open_dataset(infile[ii])
            
            
            dsraws.append(dsraw)
        dsraw = xr.concat(dsraws,dim='time')
        
        dsreg  = proc.sel_region_xr(dsraw,bbox_takahashi).load()
        sstraw = ut.standardize_names(dsreg.sst)
        sstraw = sstraw.transpose('time','lat','lon')
        
            
    
    #%% Split into time periods
    
    
    
    #%%  Generate Periods (takes awhile... 10 sec?)
    
    # Generate Periods
    varwindows,tranges = ut.generate_periods(sstraw,winlen)
    
    # Detrend by Period
    st2        = time.time()
    varwindows = ut.preprocess_byperiod(varwindows)
    print("\tCompleted period-wise detrend in %.2fs" % (time.time()-st2))
    
    #%% Looping for each period, calculate the CP and EP ENSO indices
    
    # Calculations ----------------------------------------------------------------
    
    nperiods      = len(varwindows)
    lat           = sstraw.lat
    lon           = sstraw.lon
    
    eof_byper     = []
    pc_byper      = []
    varexp_byper  = []
    trange_byper  = []
    tcenter_byper = []
    
    for pp in tqdm.tqdm(range(nperiods)):
        
        sstper  = varwindows[pp]
        sstper  = sstper.data
    
        invar   = sstper.data
        times   = varwindows[pp].time#.data
    
        ensoout = calc_enso(invar,lon,lat,pcrem,bbox=bbox_takahashi,sep_mon=sep_mon)
        eofall,pcall,varexpall = ensoout
        
        # Append Variables and Store
        eof_byper.append(eofall)
        pc_byper.append(pcall)
        varexp_byper.append(varexpall)
        
        # Get Additional Time Variables
        trange  = [str(times[0].data)[:10],str(times[-1].data)[:10],]
        tcenter = get_center_date(trange)
        trange_byper.append(trange)
        tcenter_byper.append(tcenter)
    
    
    toarr         = lambda x : np.array(x)
    eof_byper     = toarr(eof_byper)    # [Period x Lat x Lon x Mode]
    pc_byper      = toarr(pc_byper)     # [Period x Time x Mode]
    varexp_byper  = toarr(varexp_byper) # [Period x Mode]
    trange_byper  = toarr(trange_byper) # [Period x Bound]
    tcenter_byper = toarr(tcenter_byper) # [Period]
    
    # Calculate Rotated EOFs
    pc1           = pc_byper[:,:,0]
    pc2           = pc_byper[:,:,1]
    ep_byper      = (pc1 - pc2) / np.sqrt(2)
    cp_byper      = (pc1 + pc2) / np.sqrt(2)
    
    # Make into DataArrays/DataSet and Output -------------------------------------
    
    # Get Necessary Indices
    period    = np.arange(1,nperiods+1)
    pcnums    = np.arange(1,pcrem+1)
    timeindex = np.arange(len(times))
    
    # Make Coords
    coords_eofs   = dict(period=period,lat=lat,lon=lon,pc=pcnums) # 
    coords_pcs    = dict(period=period,timeindex=timeindex,pc=pcnums)
    coords_epcp   = dict(period=period,timeindex=timeindex)
    coords_varexp = dict(period=period,pc=pcnums)
    coords_trange  = dict(period=period,bounds=['start','end'])
    coords_tcenter = dict(period=period,)
        
    # Original Variables
    da_eofs       = xr.DataArray(eof_byper,coords=coords_eofs,dims=coords_eofs,name='eofs')
    da_pcs        = xr.DataArray(pc_byper,coords=coords_pcs,dims=coords_pcs,name='pcs')
    da_varexp     = xr.DataArray(varexp_byper,coords=coords_varexp,dims=coords_varexp,name='varexp')
    
    # Rotated EOFs
    da_ep         = xr.DataArray(ep_byper,coords=coords_epcp,dims=coords_epcp,name='ep')
    da_cp         = xr.DataArray(cp_byper,coords=coords_epcp,dims=coords_epcp,name='cp')
    
    # Sliding Analysis DA
    da_trange     = xr.DataArray(trange_byper,coords=coords_trange,dims=coords_trange,name='trange')
    da_tcenter    = xr.DataArray(tcenter_byper,coords=coords_tcenter,dims=coords_tcenter,name='tcenter')
    
    ds_out = xr.merge([da_eofs,da_pcs,da_varexp,da_ep,da_cp,da_trange,da_tcenter])
    
    edict = proc.make_encoding_dict(ds_out)
    ds_out.to_netcdf(savename,encoding=edict)
