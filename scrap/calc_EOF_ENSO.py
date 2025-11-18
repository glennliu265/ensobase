#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Perform EOF Analysis and retrieve the principle component timeseries

Copied function from stochmod (2025.11.18)

Created on Mon Nov 17 23:13:16 2025

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

#%% User Edits (Loop for ERA5)
infile='/home/niu4/gliu8/projects/common_data/ERA5/anom_detrend1/sst_1979_2024.nc'
outpath='/home/niu4/gliu8/projects/ccfs/enso_eof/'
expname='ERA5_1979_2024'
savename = "%s%s_ENSO_EOF.nc" % (outpath,expname)

# Dataset Information
vname    = 'sst'
timename = 'valid_time'
latname  = 'latitude'
lonname  = 'longitude'

# Calculation Options
tstart = '1979-01-01'
tend   = '2024-12-31'


#lensflag = False


bbox_takahashi = [120,360-70,-10,10]
pcrem          = 5
sep_mon        = False

if sep_mon:
    savename_in = proc.addstrtoext(savename,"_sepmon",adjust=-1)
else:
    savename_in = savename

#%% User Edits (Loop for AWI-CM3 )

vname    = "sst"
timename = 'time_counter'
latname  = 'lat'
lonname  = 'lon'

tstart   = None
tend     = None
awipath  = "/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1/"
outpath  = '/home/niu4/gliu8/projects/ccfs/enso_eof/'

expnames = ["TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]#["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
nexps = len(expnames)
for ex in range(nexps):
    
    expname  = expnames[ex]
    infile   = "%s%s_%s.nc" % (awipath,expname,'sst',)
    savename = "%s%s_ENSO_EOF.nc" % (outpath,expname)   
    if sep_mon:
        savename_in = proc.addstrtoext(savename,"_sepmon",adjust=-1)
    else:
        savename_in = savename
 
#%% Unindent this section for ERA5 Processing
    st = time.time()
    ds = xr.open_dataset(infile)['sst']
    ds = proc.format_ds(ds,timename=timename,
                        lonname=lonname,latname=latname,lon180=False)
    ds_tpac = proc.sel_region_xr(ds,bbox_takahashi)
    if tstart is not None and tend is not None: # (Don't think I have to do this)
        ds_tpac = ds_tpac.sel(time=slice(tstart,tend)).load()
    print("Loaded Output in %.2fs" % (time.time()-st))
    
    #%% Calculate ENSO
    
    
    invar   = ds_tpac.transpose('time','lat','lon').data
    lat     = ds_tpac.lat
    lon     = ds_tpac.lon
    ensoout = calc_enso(invar,lon,lat,pcrem,bbox=bbox_takahashi,sep_mon=sep_mon)
    eofall,pcall,varexpall = ensoout
    
    #%%
    
    pcnums  = np.arange(1,pcrem+1)
    times   = ds_tpac.time
    
    if sep_mon:
        mons    = np.arange(1,13,1)
        years   = np.arange(int(len(times)/12))
        
        # Make Dictionary
        coords_eofs   = dict(lat=lat,lon=lon,month=mons,pc=pcnums) # 
        coords_pcs    = dict(year=years,month=mons,pc=pcnums)
        coords_varexp = dict(month=mons,pc=pcnums)
    else:
        # Make Dictionary
        coords_eofs   = dict(lat=lat,lon=lon,pc=pcnums) # 
        coords_pcs    = dict(time=times,pc=pcnums)
        coords_varexp = dict(pc=pcnums)
    
    
    # if lensflag:
    #     ens     = np.arange(1,nens+1,1)
    #     # Unpack and repack dict to append item to start # https://www.geeksforgeeks.org/python-append-items-at-beginning-of-dictionary/
    #     coords_eofs,coords_pcs,coords_varexp = [{**{'ens':ens},**dd} for dd in [coords_eofs,coords_pcs,coords_varexp]]
    
    
    da_eofs       = xr.DataArray(eofall,coords=coords_eofs,dims=coords_eofs,name='eofs')
    da_pcs        = xr.DataArray(pcall,coords=coords_pcs,dims=coords_pcs,name='pcs')
    da_varexp     = xr.DataArray(varexpall,coords=coords_varexp,dims=coords_varexp,name='varexp')
    
    # Merge everything
    da_out        = xr.merge([da_eofs,da_pcs,da_varexp])
    
    # Add Additional Variables
    if sep_mon:
        da_out['time']      = times
    da_out['enso_bbox']     = bbox_takahashi
    
    edict = proc.make_encoding_dict(da_out)
    da_out.to_netcdf(savename_in,encoding=edict)




 