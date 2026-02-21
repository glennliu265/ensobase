#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate the PC-based NAO Index by performing an EOF analysis separately for
anomalies of each month.

Created on Fri Feb 20 14:38:38 2026

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
import matplotlib as mpl

import xeofs as xe

#%% Helper Functions

def lon360to180_xr(ds,lonname='lon'):
    # Based on https://stackoverflow.com/questions/53345442/about-changing-longitude-array-from-0-360-to-180-to-180-with-python-xarray
    ds.coords[lonname] = (ds.coords[lonname] + 180) % 360 - 180
    ds = ds.sortby(ds[lonname])
    return ds

def swap_rename(ds,chkvar,newvar):
    if chkvar in list(ds.coords):
        print("Renaming [%s] to [%s]" % (chkvar,newvar))
        ds = ds.rename({chkvar:newvar})
    return ds

def standardize_names(ds):
    
    ds = swap_rename(ds,'valid_time','time')
    ds = swap_rename(ds,'time_counter','time')
    ds = swap_rename(ds,"TIME_COUNTER",'time')
    ds = swap_rename(ds,"LON","lon")
    ds = swap_rename(ds,"LAT","lat")
    ds = swap_rename(ds,"longitude","lon")
    ds = swap_rename(ds,"latitude","lat")
    ds = swap_rename (ds,"LAT232_409","lat")
    
    # Other preprocessing
    # drop LON_bnds, TIME_COUNTER_bnds
    dropvars = ["LON_bnds","TIME_COUNTER_bnds"]
    for dropvar in dropvars:
        if dropvar in ds:
            ds = ds.drop_vars(dropvar)
    return ds

def xrdeseason(ds,check_mon=True):
    """ Remove seasonal cycle, given an Dataarray with dimension 'time'"""
    if check_mon:
        try: 
            if ds.time[0].values.item().month != 1:
                print("Warning, first month is not Jan...")
        except:
            print("Warning, not checking for feb start")
    
    return ds.groupby('time.month') - ds.groupby('time.month').mean('time')

def xrdetrend(ds,timename='time',verbose=True):
    
    detrendlin = lambda x: sp.signal.detrend(x)
    result = xr.apply_ufunc(
        detrendlin,
        ds,
        input_core_dims=[["time"]],
        output_core_dims=[["time"]],
        vectorize=True,
        dataset_fill_value=np.nan
        )
    result['time'] = ds['time']
    return result

def selmon(ds,mon):
    return ds.sel(time=ds.time.dt.month.isin([mon]))


def xr_eof(ds,N_mode,lat_weight=True):
    
    # Latitude Weighting
    if lat_weight:
        wgt             = np.sqrt(np.cos(np.radians(ds.lat.values))) # [Lat]
        dswgt           = ds * wgt[:,None,None]
    else:
        dswgt           = ds
    
    # Reshape
    dswgt           = dswgt.transpose('lat','lon','time')
    nlat,nlon,ntime = dswgt.shape
    inarr           = dswgt.data.reshape(nlat*nlon,ntime)
    
    # Perform EOF
    eofs,pcs,varexp = eof_simple(inarr,N_mode,0)
    
    # Reshape
    eofs            = eofs.reshape(nlat,nlon,N_mode)
    if lat_weight:
        eofs        = eofs * 1/wgt[:,None,None]
    
    # Make Dictionaries
    coordseof = dict(lat=ds.lat,lon=ds.lon,mode=np.arange(1,N_mode+1))
    daeof     = xr.DataArray(eofs,coords=coordseof,dims=coordseof,name="eofs")
    
    coordspc  = dict(time=ds.time,mode=np.arange(1,N_mode+1))
    dapcs     = xr.DataArray(pcs,coords=coordspc,dims=coordspc,name="pcs")

    coordsvar = dict(mode=np.arange(1,N_mode+1),)
    davarexp  = xr.DataArray(varexp,coords=coordsvar,dims=coordsvar,name="varexp")
    
    ds_eof    = xr.merge([daeof,dapcs,davarexp])
    
    return ds_eof

def eof_simple(pattern,N_mode,remove_timemean):
    """
    Simple EOF function based on script by Yu-Chiao Liang
    
    Inputs:
        1) pattern: Array of Space x Time [MxN], no NaNs
        2) N_mode:  Number of Modes to output
        3) remove_timemean: Set 1 to remove mean along N
    
    Outputs:
        1) eof: EOF patterns   [M x N_mode]
        2) pcs: PC time series [N x N_mode]
        3) varexp: % Variance explained [N_mode]
    
    Dependencies:
        import numpy as np
    
    """
    pattern1 = pattern.copy()
    nt = pattern1.shape[1] # Get time dimension size
    ns = pattern1.shape[0] # Get space dimension size
    
    if N_mode > nt:
        print("Warning, number of requested modes greater than length of time dimension. Adjusting to size of time.")
        N_mode = nt
    
    # Preallocate
    eofs = np.zeros((ns,N_mode))
    pcs  = np.zeros((nt,N_mode))
    varexp = np.zeros((N_mode))
    
    # Remove time mean if option is set
    if remove_timemean == 1:
        pattern1 = pattern1 - pattern1.mean(axis=1)[:,None] # Note, the None adds another dimension and helps with broadcasting
    
    # Compute SVD
    [U, sigma, V] = np.linalg.svd(pattern1, full_matrices=False)
    
    # Compute variance (total?)
    norm_sq_S = (sigma**2).sum()
    
    for II in range(N_mode):
        
        # Calculate explained variance
        varexp[II] = sigma[II]**2/norm_sq_S
        
        # Calculate PCs
        pcs[:,II] = np.squeeze(V[II,:]*np.sqrt(nt-1))
        
        # Calculate EOFs and normalize
        eofs[:,II] = np.squeeze(U[:,II]*sigma[II]/np.sqrt(nt-1))
    return eofs, pcs, varexp    
    


# def sel_region_xr(ds,bbox):
#     """
#     Selects region given bbox = [West Bnd, East Bnd, South Bnd, North Bnd]
    
#     Parameters
#     ----------
#     ds : xr.DataArray or Dataset
#         Assumes "lat" and "lon" variables are [present]
#     bbox : LIST
#         Boundaries[West Bnd, East Bnd, South Bnd, North Bnd]
        
#     Returns
#     -------
#         Subsetted datasetor dataarray
#     """
#     return ds.sel(lon=slice(bbox[0],bbox[1]),lat=slice(bbox[2],bbox[3]))

# def calc_eof(daf,bboxeof,N_mode=None,concat_ens=True,mask=None):
    
#     """
#     daf : time ensemble lat lon
#     sel_region_xr
#     find_nan
#     eof_simple
#     regress_2d
    
#     """
    
#     if 'ens' not in daf.dims: # Add dummy ens variable
#         daf = daf.expand_dims(dim={'ens':[0,]},axis=1)
#         print("Adding ens dim")
#     daf      = daf.transpose('time','ens','lat','lon')
#     flxa     = daf # [Time x Ens x Lat x Lon] # Anomalize variabless
    
#     # Apply area weight
#     wgt    = np.sqrt(np.cos(np.radians(daf.lat.values))) # [Lat]
#     flxwgt = flxa * wgt[None,None,:,None]
    
#     # Apply Max if needed
#     if mask is not None:
#         print("Applying provided mask...")
#         flxwgt = flxwgt * mask
    
#     # Select Region
#     flxreg     = sel_region_xr(flxwgt,bboxeof)
    
#     flxout     = flxreg.values
#     ntime,nens,nlatr,nlonr = flxout.shape
    
#     if concat_ens:
#         # IMPORTANT NOTE (implement fix later)
#         # Variable must be stacked as [ens x time x otherdims]
#         if flxout.shape[0] != nens:
#             ens_reshape_flag = True
#             print("Warning, since ensemble dimension is NOT first, temporarily permuting array to ens x time")
#             flxout = flxout.transpose(1,0,2,3)
#         else:
#             ens_reshape_flag = False
#         print("Stacking Dimensions")
#         flxout = flxout.reshape(nens*ntime,1,nlatr,nlonr)
#         ntime,nens,nlatr,nlonr = flxout.shape
#     npts       = nlatr*nlonr
#     nyr        = int(ntime/12)
#     if N_mode is None: # Set EOFs to number of years
#         N_mode=nyr
    
#     # Repeat for full variable
#     flxout_full= flxa.values
#     _,_,nlat,nlon=flxout_full.shape
#     if ens_reshape_flag:
#         print("Permuting full variable")
#         print("\tOriginal Shape %s" % str(flxout_full.shape))
#         flxout_full = flxout_full.transpose(1,0,2,3)
#         print("\tNew Shape %s" % str(flxout_full.shape))
#     npts_full  = nlat*nlon
#     if concat_ens:
#         flxout_full = flxout_full.reshape(ntime,1,nlat,nlon)
#     print("\tFinal Shape %s" % str(flxout_full.shape))
    
#     # Check to see if N_mode exceeds nyrs
#     if N_mode > nyr:
#         print("Requested N_mode exists the maximum number of years, adjusting....")
#         N_mode=nyr
    
#     # Preallocate for EOF Analysis
#     eofall    = np.zeros((N_mode,nens,nlat*nlon)) * np.nan
#     pcall     = np.zeros((N_mode,nens,ntime)) * np.nan
#     varexpall = np.zeros((N_mode,nens)) * np.nan
        
#     # Loop for ensemble memmber
#     for e in tqdm.tqdm(range(nens)):
        
#         # Remove NaN Points
#         flxens            = flxout[:,e,:,:].reshape(ntime,npts) #  Time x Space
#         okdata,knan,okpts = find_nan(flxens,0)
#         _,npts_valid = okdata.shape
        
#         # Repeat for full data
#         flxens_full       = flxout_full[:,e,:,:].reshape(ntime,npts_full)
#         okdataf,knanf,okptsf = find_nan(flxens_full,0)
#         _,npts_validf = okdataf.shape
        
#         # Reshape to [yr x mon x pts]
#         okdatar  = okdata.reshape(ntime,npts_valid)
#         okdatarf = okdataf.reshape(ntime,npts_validf)
        
#         # Calculate EOF by month
#         #for im in range(12):
            
#         # Compute EOF
#         datain          = okdatar.T # --> [space x time]
#         eofs,pcs,varexp = eof_simple(datain,N_mode,1)
        
#         # Standardize PCs
#         pcstd = pcs / pcs.std(0)[None,:]
        
#         # Regress back to dataset
#         datainf = okdatarf[:,:].T
#         eof,b = regress_2d(pcstd.T,datainf.T) # [time x pts]
        
        
#         # Save the data
#         eofall[:,e,okptsf] = eof.copy()
#         pcall[:,e,:] = pcs.T.copy()
#         varexpall[:,e] = varexp.copy()
    
#     # Reshape the variable
#     eofall = eofall.reshape(N_mode,nens,nlat,nlon) # (86, 42, 96, 89)
    
    
#     # # Flip Signs
#     # if bbox_check is not None:
#     #     print("Flipping boxes based on [bbox_check]")
#     #     nmode_check = len(bbox_check)
#     #     for N in tqdm.tqdm(range(nmode_check)):
#     #         chkbox = bbox_check[N]
#     #         for e in range(nens):
#     #             for m in range(12):
                    
                    
#     #                 sumflx = sel_region(eofall[N,[m],e,:,:].transpose(2,1,0),flxa.lon.values,flxa.lat.values,chkbox,reg_avg=True)
#     #                 #sumslp = proc.sel_region(eofslp[:,:,[m],N],lon,lat,chkbox,reg_avg=True)
                    
#     #                 if sumflx > 0:
#     #                     print("Flipping sign for NHFLX, mode %i month %i" % (N+1,m+1))
#     #                     eofall[N,m,e,:,:]*=-1
#     #                     pcall[N,m,e,:] *= -1
#     # else:
#     #     print("Sign of EOF pattern will not be checked.")
    
#     startyr   = daf.time.data[0]
#     nyrs      = int(len(daf.time)/12)
#     if concat_ens:
#         tnew      = np.arange(0,int(ntime/12))
#     else:
#         tnew      = xr.cftime_range(start=startyr,periods=nyrs,freq="YS",calendar="noleap")

#     # Make Dictionaries
#     coordseof = dict(mode=np.arange(1,N_mode+1),ens=np.arange(1,nens+1,1),lat=flxa.lat,lon=flxa.lon)
#     daeof     = xr.DataArray(eofall,coords=coordseof,dims=coordseof,name="eofs")

#     coordspc  = dict(mode=np.arange(1,N_mode+1),ens=np.arange(1,nens+1,1),yr=tnew)
#     dapcs     = xr.DataArray(pcall,coords=coordspc,dims=coordspc,name="pcs")

#     coordsvar = dict(mode=np.arange(1,N_mode+1),)
#     davarexp  = xr.DataArray(varexpall,coords=coordsvar,dims=coordsvar,name="varexp")
    
#     ds_eof    = xr.merge([daeof,dapcs,davarexp])
#     return ds_eof.squeeze()

# def sel_region(var,lon,lat,bbox,reg_avg=0,reg_sum=0,warn=1,autoreshape=False,returnidx=False,awgt=None):
#     """
    
#     Select Region
    
#     Inputs
#         1) var: ARRAY, variable with dimensions [lon x lat x otherdims]
#         2) lon: ARRAY, Longitude values
#         3) lat: ARRAY, Latitude values
#         4) bbox: ARRAY, bounding coordinates [lonW lonE latS latN]
#         5) reg_avg: BOOL, set to 1 to return regional average
#         6) reg_sum: BOOL, set to 1 to return regional sum
#         7) warn: BOOL, set to 1 to print warning text for region selection
#         8) awgt: INT, type of area weighting to apply (default is None, 1=cos(lat),2=cos^2(lat))
#     Outputs:
#         1) varr: ARRAY: Output variable, cut to region
#         2+3), lonr, latr: ARRAYs, new cut lat/lon
    
#     Assume longitude is always searching eastward...
#     Assume var is of the form [lon x lat x otherdims]
#     bbox is [lonW lonE latS latN]
    
    
#     """    
#     # Reshape to combine dimensions
#     dimflag = False 
#     if autoreshape:
#         var,vshape,dimflag=combine_dims(var,2,debug=True)
    
#     # Find indices
#     klat = np.where((lat >= bbox[2]) & (lat <= bbox[3]))[0]
#     if bbox[0] < bbox[1]:
#         klon = np.where((lon >= bbox[0]) & (lon <= bbox[1]))[0]
#     elif bbox[0] > bbox[1]:
#         if warn == 1:
#             print("Warning, crossing the prime meridian!")
#         klon = np.where((lon <= bbox[1]) | (lon >= bbox[0]))[0]
    
#     if returnidx:
#         return klon,klat
    
#     lonr = lon[klon]
#     latr = lat[klat]
    
#     #print("Bounds from %.2f to %.2f Latitude and %.2f to %.2f Longitude" % (latr[0],latr[-1],lonr[0],lonr[-1]))
    
#     # Index variable
#     varr = var[klon[:,None],klat[None,:],...]
    
#     if reg_avg==1:
#         if awgt is not None:
#             varr = area_avg(varr,bbox,lonr,latr,awgt)
#         else:
#             varr = np.nanmean(varr,(0,1))
#         return varr
#     elif reg_sum == 1:
#         varr = np.nansum(varr,(0,1))
#         return varr
    
#     # Reshape variable automatically
#     if dimflag:
#         newshape = np.hstack([[len(lonr),len(latr)],vshape[2:]])
#         varr     = varr.reshape(newshape)
    
#     return varr,lonr,latr

# def area_avg(data,bbox,lon,lat,wgt=None):
#     """
#     Function to find the area average of [data] within bounding box [bbox], 
#     based on wgt type (see inputs)
#     Inputs:
#         1) data: target array [lon x lat x otherdims]
#         2) bbox: bounding box [lonW, lonE, latS, latN]
#         3) lon:  longitude coordinate
#         4) lat:  latitude coodinate
#         5) wgt:  number (or str) to indicate weight type
#                     0 or None     no weighting
#                     1 or 'cos'  = cos(lat)
#                     2 or 'cossq' = sqrt(cos(lat))
                
    
#     Output:
#         1) data_aa: Area-weighted array of size [otherdims]
        
#     Dependencies:
#         numpy as np
    

#     """
    
#     # Check order of longitude
#     # vshape = data.shape
#     #nlon = lon.shape[0]
#     #nlat = lat.shape[0]
    
#     # Find lat/lon indices 
#     kw = np.abs(lon - bbox[0]).argmin()
#     ke = np.abs(lon - bbox[1]).argmin()
#     ks = np.abs(lat - bbox[2]).argmin()
#     kn = np.abs(lat - bbox[3]).argmin()
    
#     # Select the region
#     sel_data = data[kw:ke+1,ks:kn+1,:]
    
#     # If wgt == 1, apply area-weighting 
#     if (wgt != 0) or (wgt is not None):
        
#         # Make Meshgrid
#         _,yy = np.meshgrid(lon[kw:ke+1],lat[ks:kn+1])
        
#         # Calculate Area Weights (cosine of latitude)
#         if (wgt == 1) or (wgt == "cos"):
#             wgta = np.cos(np.radians(yy)).T
#         elif (wgt == 2) or (wgt == "cossq"):
#             wgta = np.sqrt(np.cos(np.radians(yy))).T
        
#         # Remove nanpts from weight, ignoring any pt with nan in otherdims
#         nansearch = np.sum(sel_data,2) # Sum along otherdims
#         wgta[np.isnan(nansearch)] = 0
        
#         # Apply area weights
#         #data = data * wgtm[None,:,None]
#         sel_data  = sel_data * wgta[:,:,None]
    
#     # Take average over lon and lat
#     if (wgt != 0) or (wgt is not None):

#         # Sum weights to get total area
#         sel_lat  = np.sum(wgta,(0,1))
        
#         # Sum weighted values
#         data_aa = np.nansum(sel_data/sel_lat,axis=(0,1))
#     else:
#         # Take explicit average
#         data_aa = np.nanmean(sel_data,(0,1))
#     return data_aa

# def combine_dims(var,nkeep,debug=True):
#     """
#     Keep first n dimensions of var, and 
#     reshape to combine the rest

#     Parameters
#     ----------
#     var : ARRAY
#         Array to reshape
#     nkeep : INT
#         Number of dimensions to avoid reshape
#     debug : BOOL (optional)
#         Prints warning message if reshaping occurs
        
    
#     Returns
#     -------
#     var : ARRAY
#         Reshaped Variable
#     vshape : TYPE
#         Original shape of the variable
#     dimflag : BOOL
#         True if variable is reshaped

#     """
#     dimflag = False
#     vshape  = var.shape
#     if (len(var.shape) > nkeep+1):
#         dimflag=True
#         if debug:
#            print("Warning, variable has more than %i dimensions, combining last!"% (nkeep+1))
#         otherdims = np.prod(vshape[nkeep:])
#         newshape = np.hstack([vshape[:nkeep],[otherdims]])
#         var = var.reshape(vshape[0],vshape[1],otherdims)
#     return var,vshape,dimflag


# def find_nan(data,dim,val=None,return_dict=False,verbose=True):
#     """
#     For a 2D array, remove any point if there is a nan in dimension [dim].
    
#     Inputs:
#         1) data        : 2d array, which will be summed along last dimension
#         2) dim         : dimension to sum along. 0 or 1.
#         3) val         : value to search for (default is NaN)
#         4) return_dict : Set to True to return dictionary with clearer arguments...
#     Outputs:
#         1) okdata : data with nan points removed
#         2) knan   : boolean array with indices of nan points
#         3) okpts  : indices for non-nan points
#     """
    
#     # Sum along select dimension
#     if len(data.shape) > 1:
#         datasum = np.sum(data,axis=dim)
#     else:
#         datasum = data.copy()
    
#     # Find non nan pts
#     if val is None:
#         knan  = np.isnan(datasum)
#     else:
#         knan  = (datasum == val)
#     okpts = np.invert(knan)
    
#     if len(data.shape) > 1:
#         if dim == 0:
#             okdata = data[:,okpts]
#             clean_dim = 1
#         elif dim == 1:    
#             okdata = data[okpts,:]
#             clean_dim = 0
#     else:
#         okdata = data[okpts]
#     if verbose:
#         print("Found %i NaN Points along axis %i." % (data.shape[clean_dim] - okdata.shape[clean_dim],clean_dim))
#     if return_dict: # Return dictionary with clearer arguments
#         nandict = {"cleaned_data" : okdata,
#                    "nan_indices"  : knan,
#                    "ok_indices"   : okpts,
#                    }
#         return nandict
#     return okdata,knan,okpts

# def regress_2d(A,B,nanwarn=1,verbose=True):
#     """
#     Regresses A (independent variable) onto B (dependent variable), where
#     either A or B can be a timeseries [N-dimensions] or a space x time matrix 
#     [N x M]. Script automatically detects this and permutes to allow for matrix
#     multiplication.
#     Note that if A and B are of the same size, assumes axis 1 of A will be regressed to axis 0 of B
    
#     Returns the slope (beta) for each point, array of size [M]
    
    
#     """
    
#     # Determine if A or B is 2D and find anomalies
#     bothND = False # By default, assume both A and B are not 2-D.
#     # Note: need to rewrite function such that this wont be a concern...
    
#     # Accounting for the fact that I dont check for equal dimensions below..
#     #B = B.squeeze()
#     #A = A.squeeze() Commented out below because I still need to fix some things
#     # Compute using nan functions (slower)
#     if np.any(np.isnan(A)) or np.any(np.isnan(B)):
#         if nanwarn == 1:
#             print("NaN Values Detected...")
        
#         # 2D Matrix is in A [MxN]
#         if len(A.shape) > len(B.shape):
            
#             # Tranpose A so that A = [MxN]
#             if A.shape[1] != B.shape[0]:
#                 A = A.T
            
#             # Set axis for summing/averaging
#             a_axis = 1
#             b_axis = 0
            
#             # Compute anomalies along appropriate axis
#             Aanom = A - np.nanmean(A,axis=a_axis)[:,None]
#             Banom = B - np.nanmean(B,axis=b_axis)
            
#         # 2D matrix is B [N x M]
#         elif len(A.shape) < len(B.shape):
            
#             # Tranpose B so that it is [N x M]
#             if B.shape[0] != A.shape[0]:
#                 B = B.T
            
#             # Set axis for summing/averaging
#             a_axis = 0
#             b_axis = 0
            
#             # Compute anomalies along appropriate axis        
#             Aanom = A - np.nanmean(A,axis=a_axis)
#             Banom = B - np.nanmean(B,axis=b_axis)[None,:]
        

#         # A is [P x N], B is [N x M]
#         elif len(A.shape) == len(B.shape):
#             if verbose:
#                 print("Note, both A and B are 2-D...")
#             bothND = True
#             if A.shape[1] != B.shape[0]:
#                 print("WARNING, Dimensions not matching...")
#                 print("A is %s, B is %s" % (str(A.shape),str(B.shape)))
#                 print("Detecting common dimension")
#                 # Get intersecting indices 
#                 intersect, ind_a, ind_b = np.intersect1d(A.shape,B.shape, return_indices=True)
#                 if ind_a[0] == 0: # A is [N x P]
#                     A = A.T # Transpose to [P x N]
#                 if ind_b[0] == 1: # B is [M x N]
#                     B = B.T # Transpose to [N x M]
#                 print("New dims: A is %s, B is %s" % (str(A.shape),str(B.shape)))
                
#             # Set axis for summing/averaging
#             a_axis = 1 # Assumes dim 1 of A will be regressed to dim 0 of b
#             b_axis = 0
            
#             # Compute anomalies along appropriate axis        
#             Aanom = A - np.nanmean(A,axis=a_axis,keepdims=True)#[:,None] # Anomalize w.r.t. dim 1 of A
#             Banom = B - np.nanmean(B,axis=b_axis,keepdims=True)# # Anonalize w.r.t. dim 0 of B
            
#         # Calculate denominator, summing over N
#         Aanom2 = np.power(Aanom,2)
#         denom  = np.nansum(Aanom2,axis=a_axis,keepdims=True)     # Sum along dim 1 of A (lets say this is time)
        
#         # Calculate Beta
#         #if 
#         if len(denom.shape)==1 or not bothND: # same as both not ND
#             print("Adding singleton dimension to denom")
#             denom = denom[:,None]
#         beta = Aanom @ Banom / denom#[:,None] # Denom is [A[mode,time]@ B[time x space]], output is [mode x pts]
        
#         b = (np.nansum(B,axis=b_axis,keepdims=True) - beta * np.nansum(A,axis=a_axis,keepdims=True))/A.shape[a_axis]
#         # b is [mode x pts] [or P x M]
            
#     else:
#         # 2D Matrix is in A [MxN]
#         if len(A.shape) > len(B.shape):
            
#             # Tranpose A so that A = [MxN]
#             if A.shape[1] != B.shape[0]:
#                 A = A.T
                
#             a_axis = 1
#             b_axis = 0
            
#             # Compute anomalies along appropriate axis
#             Aanom = A - np.mean(A,axis=a_axis)[:,None]
#             Banom = B - np.mean(B,axis=b_axis)
            
#         # 2D matrix is B [N x M]
#         elif len(A.shape) < len(B.shape):
            
#             # Tranpose B so that it is [N x M]
#             if B.shape[0] != A.shape[0]:
#                 B = B.T
            
#             # Set axis for summing/averaging
#             a_axis = 0
#             b_axis = 0
            
#             # Compute anomalies along appropriate axis        
#             Aanom = A - np.mean(A,axis=a_axis)
#             Banom = B - np.mean(B,axis=b_axis)[None,:]
            
#         # A is [P x N], B is [N x M]
#         elif len(A.shape) == len(B.shape):
#             if verbose:
#                 print("Note, both A and B are 2-D...")
#             bothND = True
#             if A.shape[1] != B.shape[0]:
#                 print("WARNING, Dimensions not matching...")
#                 print("A is %s, B is %s" % (str(A.shape),str(B.shape)))
#                 print("Detecting common dimension")
#                 # Get intersecting indices 
#                 intersect, ind_a, ind_b = np.intersect1d(A.shape,B.shape, return_indices=True)
#                 if ind_a[0] == 0: # A is [N x P]
#                     A = A.T # Transpose to [P x N]
#                 if ind_b[0] == 1: # B is [M x N]
#                     B = B.T # Transpose to [N x M]
#                 print("New dims: A is %s, B is %s" % (str(A.shape),str(B.shape)))
            
#             # Set axis for summing/averaging
#             a_axis = 1
#             b_axis = 0
            
#             # Compute anomalies along appropriate axis        
#             Aanom = A - np.mean(A,axis=a_axis)[:,None]
#             Banom = B - np.mean(B,axis=b_axis)[None,:]

#         # Calculate denominator, summing over N
#         Aanom2 = np.power(Aanom,2)
#         denom  = np.sum(Aanom2,axis=a_axis,keepdims=True)
#         if not bothND:
            
#             denom = denom[:,None] # Broadcast
            
#         # Calculate Beta
#         beta = Aanom @ Banom / denom
            
#         if bothND:
#             b = (np.sum(B,axis=b_axis)[None,:] - beta * np.sum(A,axis=a_axis)[:,None])/A.shape[a_axis]
#         else:
#             b = (np.sum(B,axis=b_axis) - beta * np.sum(A,axis=a_axis))/A.shape[a_axis]
    
#     return beta,b

    
#%% Input Information

use_xeof = False # Currently doesn't work
bbox     = [-90,40,20,80] # Assumes with degrees west #[-90+360, 40, 20, 80]

# Information for SLP Data (Input)
vname    = "msl" # Name of the variable
datpath  = "/home/niu4/gliu8/projects/scrap/processed_global/" # Input Path
expname  = "TCo319-DART-ssp585d-gibbs-charn" # Name of the experiment (for naming)
expnames = ["TCo319_ctl1950d",
            "TCo319_ssp585",
            "TCo319-DART-ctl1950d-gibbs-charn",
            "TCo319-DART-ssp585d-gibbs-charn",
            "TCo95-hi1950d",
            "TCo95-ssp585d",
            "TCo1279-DART-2060",
            "TCo1279-DART-1950",
            "TCo1279-DART-2090"] # ,,,

for expname in expnames:
    
    # Additional Input Information (for expname loop)
    ncname   = "%s_%s.nc" % (expname,vname)# Name of Netcdf
    
    # Info for output
    outpath     = "/home/niu4/gliu8/projects/scrap/nao_indices/"
    ncname_out  = "%s%s_NAO_Indices.nc" % (outpath,expname,)
    ncname_out_allmon  = "%s%s_NAO_Indices_AllMonths.nc" % (outpath,expname,)
    
    
    
    #%% Load the data
    
    # Open a View
    ds      = xr.open_dataset(datpath+ncname)[vname]
    
    # Select the NAO Region
    if np.any(ds.lon.data > 180):
        print("Flipping Longitude...")
        ds      = lon360to180_xr(ds) # Correct Longitude (assuming it has values over 360)
    dsreg   = ds.sel(lon=slice(bbox[0],bbox[1]),lat=slice(bbox[2],bbox[3]))
    
    # Load the variable (took ~400 seconds)
    st      = time.time()
    dsreg   = dsreg.load()
    print("Loaded in %.2fs" % (time.time()-st))
    
    dsreg   = standardize_names(dsreg)
    dsreg   = dsreg / 100 # Convert to HectoPascals
    
    #%% Perform Preprocessing (Detrend, Deseason)
    
    # Remove seasonal cycle and detrend
    st      = time.time()
    dsa     = xrdeseason(dsreg) # Remove mean seasonal cycle
    dsa_dt  = xrdetrend(dsa) # Remove simple linear trend
    print("Deseason/Detrend in %.2fs" % (time.time()-st))
    
    #%% Do NAO Calculations
    
    N_mode    = 3
    months    = np.arange(1,13,1)
    
    eofmon    = []
    pcmon     = []
    varexpmon = []
    for mon in tqdm.tqdm(months):
        dsmon = selmon(dsa_dt,mon)
        
        
        # Use xEOFs to compute necessary information
        if use_xeof:
            model           = xe.single.EOF(use_coslat=True,n_modes=N_mode)
            st              = time.time()
            model.fit(dsmon,dim='time') # Convert to hectopascals
            
            sigma     = model.singular_values()
            norm_sq_S = (sigma**2).sum('mode')
            nt        = len(dsmon.time)
            
            eofcorr   = sigma/np.sqrt(nt)
            pccorr    = np.sqrt(nt-1)
            
            
            eofall          = model.components() # [N_mode x lat x lon]
            pcall           = model.scores()
            varexpall       = model.explained_variance_ratio()
            print("Computed EOF in %.2fs" % (time.time()-st))
        
            # Need to set this from {} --> 'none' for to_netcdf to work later
            eofall.attrs['solver_kwargs']='none'
            pcall.attrs['solver_kwargs']='none'
            varexpall.attrs['solver_kwargs']='none'
        else:
            ds_eof = xr_eof(dsmon,N_mode=N_mode,lat_weight=True)
            
            eofall = ds_eof.eofs.transpose('mode','lat','lon')
            pcall  = ds_eof.pcs.transpose('mode','time')
            varexpall = ds_eof.varexp
        
        # Flip Signs where necessary
        spgbox     = [-60,20,45,80]
        eapbox     = [-60,20,45,60] # Shift Box west for EAP
        bbox_check = [spgbox,eapbox,]    
        print("Flipping boxes based on [bbox_check]")
        nmode_check = len(bbox_check)
        for N in range(nmode_check):
            chkbox = bbox_check[N]
            
            sumflx = eofall.isel(mode=N).sel(lon=slice(chkbox[0],chkbox[1]),lat=slice(chkbox[2],chkbox[3])).mean().data.item()
            
            if sumflx > 0:
                print("Flipping sign for SLP, mode %i" % (N+1))
                eofall[N,:,:] *= -1
                pcall[N,:] *= -1
        
        pcall['month']     = mon
        eofall['month']    = mon
        varexpall['month'] = mon
        
        # Append Output
        eofmon.append(eofall)
        pcmon.append(pcall)
        varexpmon.append(varexpall)
    
    # Concatenate and Save Output
    eofmon    = xr.concat(eofmon,dim='month')
    pcmon     = xr.concat(pcmon,dim='month')
    varexpmon = xr.concat(varexpmon,dim='month')
    dsout     = xr.merge([eofmon.rename('eof'),
                      pcmon.rename('pc'),
                      varexpmon.rename('varexp')])
    
    dsout.to_netcdf(ncname_out)
        
                
    # Repeat Again for All Months Together ========================================
    
    # Use xEOFs to compute necessary information
    if use_xeof:
        model           = xe.single.EOF(use_coslat=True,n_modes=N_mode)
        st              = time.time()
        model.fit(dsa_dt,dim='time')
        eofall          = model.components()
        pcall           = model.scores()
        varexpall       = model.explained_variance_ratio()
        print("Computed EOF in %.2fs" % (time.time()-st))
    else:
        ds_eof = xr_eof(dsa_dt,N_mode=N_mode,lat_weight=True)
        eofall = ds_eof.eofs.transpose('mode','lat','lon')
        pcall  = ds_eof.pcs.transpose('mode','time')
        varexpall = ds_eof.varexp
    
    
    # Flip Signs where necessary
    spgbox     = [-60,20,45,80]
    eapbox     = [-60,20,45,60] # Shift Box west for EAP
    bbox_check = [spgbox,eapbox,]    
    print("Flipping boxes based on [bbox_check]")
    nmode_check = len(bbox_check)
    for N in range(nmode_check):
        chkbox = bbox_check[N]
        
        sumflx = eofall.isel(mode=N).sel(lon=slice(chkbox[0],chkbox[1]),lat=slice(chkbox[2],chkbox[3])).mean().data.item()
        
        if sumflx > 0:
            print("Flipping sign for SLP, mode %i" % (N+1))
            eofall[N,:,:] *= -1
            pcall[N,:] *= -1
    
    if use_xeof:
        #Need to set this from {} --> 'none' for to_netcdf to work later
        eofall.attrs['solver_kwargs']='none'
        pcall.attrs['solver_kwargs']='none'
        varexpall.attrs['solver_kwargs']='none'
    
    dsout     = xr.merge([eofall.rename('eof'),
                      pcall.rename('pc'),
                      varexpall.rename('varexp')])
    dsout.to_netcdf(ncname_out_allmon)
    
    # #%%
    # # # Apply area weight
    # dsa_dt = dsa_dt.transpose('time','lat','lon')
    # wgt    = np.sqrt(np.cos(np.radians(dsa_dt.lat.values))) # [Lat]
    # # dswgt  = dsa_dt.data * wgt[None,:,None]






        



