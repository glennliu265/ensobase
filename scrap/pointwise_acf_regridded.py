#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate the pointwise ACF for regridded files

Created on Wed Dec 10 14:29:39 2025

@author: gliu


"""

import numpy as np
import scipy as sp
import xarray as xr
import scipy as sp
from scipy import signal,stats


#%%
def calc_lag_covar_ann(var1,var2,lags,dim,detrendopt,verbose=True):
    """
    Calculate lag **correlation**
    var1 is lagged, var2 remains the same (base)
    
    """
    
    # Move time to the first dimension (assume var1.shape==var2.shape)
    invars      = [var1,var2]
    oldshape    = var1.shape
    reshapevars = []
    for v in invars:
        vreshape,neworder = dim2front(v,dim,combine=True,return_neworder=True)
        reshapevars.append(vreshape)
    
    # Remove Nan Points (if any are found)
    
    # Get total number of lags
    var1,var2=reshapevars
    lagdim   = len(lags)
    if lagdim > var1.shape[0]:
        if verbose:
            print("\tWarning! maximum lag  (%i) exceeds the timeseries length (%i)" % (lagdim,var1.shape[0]))
    
    # Get timeseries length # [yr x npts]
    ntime  = var1.shape[0]
    npts   = var1.shape[1]
    
    # Detrend variables if option is set
    if detrendopt == 1:
        if verbose:
            print("\tWarning! Variable will be detrended linearly.")
        var1 = signal.detrend(var1,0,type='linear')
        var2 = signal.detrend(var2,0,type='linear')
    
    # Preallocate
    corr_ts        = np.zeros((lagdim,npts)) * np.nan
    window_lengths = []
    for l,lag in enumerate(lags):
        varlag   = var1[lag:,:]
        varbase  = var2[:(ntime-lag),:]
        window_lengths.append(varlag.shape[0])
        
        # Calculate correlation
        corr_ts[l,:] = pearsonr_2d(varbase,varlag,0)    
    
    # Replace back into old shape
    size_combined_dims = tuple(np.array(oldshape)[neworder][1:]) # Get other dims
    reshape_corr       = (lagdim,) + size_combined_dims
    corr_ts            = corr_ts.reshape(reshape_corr)
    
    return corr_ts,window_lengths


def dim2front(x,dim,verbose=True,combine=False,flip=False,return_neworder=False):
    """
    Move dimension in position [dim] to the front
    
    Parameters
    ----------
    x    [NDARRAY] : Array to Reorganize
    dim      [INT] : Target Dimension
    combine [BOOL] : Set to True to combine all dims
    flip    [BOOL] : Reverse the dimensions
    
    Returns
    -------
    y [ND Array]   : Output

    """
    if dim <0:
        dim = len(x.shape) + dim
        
    neworder = np.concatenate([[dim,],
                         np.arange(0,dim),
                         np.arange(dim+1,len(x.shape))
                         ])
    y = x.transpose(neworder)
    if combine:
        y = y.reshape(y.shape[0],np.prod(y.shape[1:]))
    if flip:
        y = y.transpose(np.flip(np.arange(0,len(y.shape))))
    if verbose:
        print("New Order is : %s"%str(neworder))
        print("Dimensions are : %s"%str(y.shape))
    if return_neworder:
        return y,neworder
    return y


def pearsonr_2d(A,B,dim,returnsig=0,p=0.05,tails=2,dof='auto'):
    """
    Calculate Pearson's Correlation Coefficient for two 2-D Arrays
    along the specified dimension. Input arrays are anomalized.
    
    Option to perform students t-test.
    
    Inputs
    -------
    1) A : ARRAY
        First variable, 2D
    2) B : ARRAY
        Second variable, 2D (same axis arrangement as A)
    3) dim : INT
        Dimension to compute correlation along
    OPTIONAL ARGS
    4) returnsig : BOOL
        Return significance test result    
    5) p : FLOAT
        P-value for significance testing
    6) tails: INT
        Number of tails (1 or 2)
    7) dof: "auto" or INT
        Degress of freedom method. Set to "auto" for ndim - 2, or
        manually enter an integer value.

    Outputs
    --------
    1) rho : ARRAY
        Array of correlation coefficients
    OPTIONAL OUTPUTS (if returnsig=1)
    2) T : ARRAY
        T values from testing
    3) critval : FLOAT
        T - Critical value used as threshold
    4) sigtest : ARRAY of BOOLs
        Indicates points that passed significance threshold
    5) corrthres : FLOAT
        Correlation threshold corresponding to critval
    
    
    Calculates correlation between two matrices of equal size and also returns
    the significance test result
    
    Dependencies
        numpy as np
        scipy stats
    
    """
    
    # Find Anomaly
    Aanom = A - np.nanmean(A,dim)
    Banom = B - np.nanmean(B,dim)
    
    # Elementwise product of A and B
    AB = Aanom * Banom
    
    # Square A and B
    A2 = np.power(Aanom,2)
    B2 = np.power(Banom,2)
    
    # Compute Pearson's Correlation Coefficient
    rho = np.nansum(AB,dim) / np.sqrt(np.nansum(A2,dim)*np.nansum(B2,dim))
    
    if returnsig == 0:
        return rho
    else:
        
        # Determine DOF (more options to add later...)
        if dof == 'auto':
            # Assume N-2 dof
            n_eff = A.shape[dim]-2
        else:
            # Use manually supplied dof
            n_eff = dof
        
        # Compute p-value based on tails
        ptilde = p/tails
        
        # Compute T at each point
        T = rho * np.sqrt(n_eff / (1 - np.power(rho,2)))
        
        # Get threshold critical value
        critval = stats.t.ppf(1-ptilde,n_eff)
        
        # Perform test
        sigtest = np.where(np.abs(T) > critval)
        
        # Get critical correlation threshold
        corrthres = np.sqrt(1/ ((n_eff/np.power(critval,2))+1))
        
        return rho,T,critval,sigtest,corrthres

#%% # Load Anomalized and Detrended Datasets

import tqdm
expnames = ["TCo319_ctl1950d", "TCo319_ssp585",
            "TCo1279-DART-1950", "TCo1279-DART-2090", "TCo2559-DART-1950C"]
            
lagfit = np.arange(0,7)


for ex in range(5):
    
    datpath = "/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/"
    nc      = "%s_sst_regrid1x1.nc" % expnames[ex]
    ds      = xr.open_dataset(datpath+nc)['sst'].load()
    
    # calculate acf
    sst     = ds.data 
    lags    = np.arange(37)
    acf,window = calc_lag_covar_ann(sst,sst,lags,dim=0,detrendopt=0)
    
    # fit exponential
    def expfit(acf,lags,lagmax):
        # Copief from reemergence/estimate_damping_fit/12.S992 Final Project
        expf3      = lambda t,b: np.exp(b*t)         # No c and A
        funcin     = expf3
        x          = lags
        y          = acf
        popt, pcov = sp.optimize.curve_fit(funcin, x[:(lagmax+1)], y[:(lagmax+1)],maxfev=5000)
        
        tau_inv = popt[0] # 1/tau (tau = timescale),. np.exp(tau_inv*t)
        acf_fit = expf3(lags,tau_inv)
        outdict = {'tau_inv':tau_inv, 'acf_fit':acf_fit}
        return outdict
    
    _,nlat,nlon = acf.shape
    timescales  = np.zeros((nlat,nlon))
    for a in tqdm.tqdm(range(nlat)):
        for o in range(nlon):
            acfin = acf[:,a,o]
            if np.any(np.isnan(acfin)):
                continue
            outdict = expfit(acfin,lagfit,lagfit[-1])
            timescales[a,o] = 1/outdict['tau_inv'] * -1
    
    
    coords = dict(lat=ds.lat,lon=ds.lon)
    ds_tau = xr.DataArray(timescales,coords=coords,dims=coords,name='tau')
    outpath = "/home/niu4/gliu8/projects/awi_hackathon/acf_regrid/"
    outname = "%s%s_sst_acf.nc" % (outpath,expnames[ex])
    
    ds_tau.to_netcdf(outname)

# ds_lat = xr.DataArray(data.lat,coords=coords,dims=coords,name='lat')
# ds_lon = xr.DataArray(data.lat,coords=coords,dims=coords,name='lon')

# ds_out = xr.merge([ds_tau,data.lat,data.lon])
# ds
#ds_tau.to_netcdf()
#ds_tau.plot(vmin=0,vmax=6,cmap="RdBu_r"),plt.show()













#%%


