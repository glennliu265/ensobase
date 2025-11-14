#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate CCFs (in a rather coarse fashion...)
i.e. I dont have wind speed so using wind stress
and perform on 1x1 degree data that has been regridded...

Created on Fri Nov  7 14:55:27 2025

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


from sklearn.linear_model import LinearRegression
import sklearn


#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% Other Functions


def mlr(X,y):
    # MLR fit using scipy 
    
    # Initialize Model and Fit
    model             = LinearRegression()
    model.fit(X,y)
    pred              = model.predict(X)
    # Calculate Error and other variables
    mlr_out           = {}
    mlr_out['pred']   = pred
    mlr_out['err']    = y - pred # Model Error # []
    mlr_out['coeffs'] = model.coef_
    mlr_out['r2']     = sklearn.metrics.r2_score(y,pred)
    
    return mlr_out

def mlr_ccfs(ccfs,flx,standardize=True,fill_value=0,verbose=False):
    # Perform MLR
    #    ccfs: LIST of DataArrays [variable x time]
    #    flx:  DataArray [time x 1]
    
    # Set up Predictors and Target (convert to Numpy Arrays)
    predictors = np.array([ds for ds in ccfs]) # [variable x time]
    if standardize:
        if verbose:
            print("Standardizing each variable")
        predictors = np.array([ds/np.nanstd(ds) for ds in list(predictors)])
    X = predictors.T
    y = flxpt.data
    
    # Replace NaN Values in Predictors
    if verbose:
        if np.any(np.isnan(X)):
            print("NaN values detected! Replace with %f" % fill_value)
    X = np.where(np.isnan(X),fill_value,X) # Set NaN to zero
    
    # Use sklearn for now (can try LSE manual later...)
    mlr_out = mlr(X,y)
    return mlr_out

#%% Load Land Mask


landmask = ut.load_land_mask_awi("TCo319",regrid=True)
expnames = ["TCo1279-DART-1950","TCo2559-DART-1950C"]
datpath  = "/home/niu4/gliu8/projects/scrap/regrid_1x1/"

# ===================================
#%% Part (1)" Calculate Wind Speed and Tadv ===================================
# ===================================

"""

Note, skip this part if Tadv and WS are already calculated...

"""

#% Load input

vnames   = ["tx_sur","ty_sur","sst"] #Dont change order

dsbyexp = []
for ex in range(2):
    dsvars = []
    for v in range(3):
        ncname = "%s%s_%s_regrid1x1.nc" % (datpath,expnames[ex],vnames[v])
        ds = xr.open_dataset(ncname)[vnames[v]].load()
        ds = ut.standardize_names(ds)
        
            
        
        print(ds.shape)
        dsvars.append(ds)
    if ex == 1:
        dstx,dsst =  proc.match_time_month(dsvars[0],dsvars[2])
        dsty,dsst =  proc.match_time_month(dsvars[1],dsvars[2])
        dsvars[0] = dstx
        dsvars[1] = dsty
        dsvars[2] = dsst

    #dsvars = xr.merge(dsvars)
    
    
    dsbyexp.append(dsvars)
#%% Check the weird units for wind stresss
taux = dsvars[0].mean('time')
tauy = dsvars[1].mean('time')
    
#%% Calculation of Wind Speed

# Try a constant value of CD :( ....
CD      = 1.25e-3   # [dimensionless], from Kara et al. 2007
rho_air = 1.225     # [kg/m3], from Li et al. 2020 (Laifang's 2020 paper)

def calc_ws(ds):
    return ((ds[0])**2 + (ds[1])**2)**0.5

dsws = [calc_ws(ds) for ds in dsbyexp]

#%% Calculation of advection

"""



"""
ex        = 0
RE        = 6375e3 # Earth Radius, in meters

ds_Tadv = []
for ex in range(2):
    st        = time.time()
    
    U10       = dsbyexp[ex][0] / (rho_air*CD) # Using zonal Wind Stress for Now ...
    V10       = dsbyexp[ex][1] / (rho_air*CD) # Using meridional Wind Stress for Now ...
    phi       = np.radians(dsbyexp[ex][0].lat)
    lbd       = np.radians(dsbyexp[ex][0].lon)
    print("\tComputed Variables")
    
    dSST_dlbd = dsbyexp[ex][2].differentiate('lon')
    dSST_dphi = dsbyexp[ex][2].differentiate('lat')
    print("\tDifferentiated SST along Lon/Lat")
    
    # Term1 = U10/(RE*np.cos(phi)) * dSST_dlbd
    # #Term1 = U10 * dSST_dlbd
    # print("\tCaclulated 1st Term  (%.2fs)" % (time.time()-st))
    # Term2 = V10/RE               * dSST_dphi
    # print("\tCaclulated 2nd Term  (%.2fs)" % (time.time()-st))
    # Tadv  = - Term1 - Term2
    # print("\tCaclulated Tadv  (%.2fs)" % (time.time()-st))
    # ds_Tadv.append(Tadv)
    
    u10 = U10.data
    dT  = dSST_dlbd.data
    
    #Term1 = U10.data/(np.cos(phi.data))[None,:,None] * dSST_dlbd.data
    
    Term1 = U10.data / (np.cos(phi.data))[None,:,None] * dSST_dlbd.data
    
    print("\tCaclulated 1st Term  (%.2fs)" % (time.time()-st))
    Term2 = V10.data * dSST_dphi.data
    
    print("\tCaclulated 2nd Term  (%.2fs)" % (time.time()-st))
    Tadv  = (- Term1 - Term2)/RE
    print("\tCaclulated Tadv  (%.2fs)" % (time.time()-st))
    
    coords = dict(time=dsbyexp[ex][0].time,lat=dsbyexp[ex][0].lat,lon=dsbyexp[ex][0].lon,)
    Tadv = xr.DataArray(Tadv,coords=coords,dims=coords,name='Tadv')
    
    ds_Tadv.append(Tadv)

#%% Save each output

vnames_out = ["WS","Tadv"]
dsouts = [dsws,ds_Tadv]

for ex in range(2):
    
    for v in range (2):
        
        ncname = "%s%s_%s_regrid1x1.nc" % (datpath,expnames[ex],vnames_out[v])
        dsout  = dsouts[v][ex].rename(vnames_out[v])
        dsout.to_netcdf(ncname)
        
#%% Try Looking at time mean values

ds_Tadv[1].mean('time').plot(vmin=-2.5,vmax=2.5),plt.show()

tmean = ds_Tadv[1].mean('time')
dtday = 3600*24
(tmean*dtday).plot(vmin=-2.5,vmax=2.5,cmap='cmo.balance'),plt.show()

wsmean = dsws[0].mean('time')
wsmean.plot(vmin=.03,vmax=.13,cmap='cmo.balance'),plt.show()
(dsws[1].mean('time')).plot(vmin=3,vmax=13,cmap='cmo.balance'),plt.show()

# ===================================
# Part (2): Compute Radiative Kernels
# ===================================
#%% Now Load each variable to do CCFs calculation

ccf_vars = ["sst","eis","Tadv","r700","w700","WS","ucc"] 

def reduce_time(ds,dsst):
    dsnew,dsst = proc.match_time_month(ds,dsst)
    return dsnew

dsbyexp = []
for ex in range(2):
    dsvars = []
    for v in range(len(ccf_vars)):
        
        ncname = "%s%s_%s_regrid1x1.nc" % (datpath,expnames[ex],ccf_vars[v])
        ds     = xr.open_dataset(ncname)[ccf_vars[v]].load()
        ds     = ut.standardize_names(ds)
        
        if expnames[ex] == 'TCo2559-DART-1950C':
            if ccf_vars[v] == "w700": # Duplicate a month for the missing variables
                w700data = ds.data.squeeze()
                w700data_duplicate_jan1950 = np.concatenate([w700data[[0],...],w700data],axis=0)
                newtime = dsvars[v-1].time
                coords  = dict(time=newtime,lat=ds.lat,lon=ds.lon)
                w700new = xr.DataArray(w700data_duplicate_jan1950,coords=coords,dims=coords,name='w700')
                ds = w700new
        
        print("%s, %s" % (expnames[ex],ccf_vars[v]))
        print(ds.shape)
        print("")
        dsvars.append(ds.squeeze())
    
    # SST is only 8 years, so reduce the time...
    if expnames[ex] == 'TCo2559-DART-1950C': 
        dsvars = [reduce_time(ds,dsvars[0]) for ds in dsvars]
    
    dsbyexp.append(dsvars)




#%% Anomalize and detrend variables

dsbyexp_anoms = []
for ex in range(2):
    dsvars_anoms = []
    
    for v in tqdm.tqdm(range(len(ccf_vars))):
        
        dsin    = dsbyexp[ex][v]
        dsanoms = ut.preprocess_enso(dsin)
        dsvars_anoms.append(dsanoms)
    
    dsbyexp_anoms.append(dsvars_anoms)

#%% Looping for each radiative fluxes

flxnames = ['allsky','clearsky','cre'] # ['cre',]

# MLR Calculation Options
standardize = True
fill_value  = 0
add_ucc     = True



for flxname in flxnames:
    
    # Load Fluxes for each experiment
    dsflxs  = []
    for ex in range(2):
        dsflx   = xr.open_dataset("%s%s_%s_regrid1x1.nc" % (datpath,expnames[ex],flxname)).load()
        dsflx   = ut.preprocess_enso(ut.standardize_names(dsflx[flxname]))
        dsflx   = ut.varcheck(dsflx,flxname,expnames[ex])
        
        dsflxs.append(dsflx)
        
    # Perform MLR for each experiment
    for ex in range(2):
        st = time.time()
        
        dsexp_sel      = dsbyexp_anoms[ex]
        if add_ucc is False:
            dsexp_sel = dsexp_sel[:-1]
        dsexp_flx      = dsflxs[ex]
        
        # Check time dimension
        ntimes_predictors = [len(ds.time) for ds in dsexp_sel]
        ntimes_flux       = len(dsexp_flx.time)
        if expnames[ex] == 'TCo2559-DART-1950C':
            print("Adjusting flux length")
            dsexp_flx,_ = proc.match_time_month(dsexp_flx,dsexp_sel[0])
        
        # Pre-allocate
        lon             = dsexp_flx.lon.data
        lat             = dsexp_flx.lat.data
        nlat,nlon,ntime = dsexp_flx.shape
        nccfs           = len(dsexp_sel)
        coeffs          = np.zeros((nlat,nlon,nccfs)) * np.nan # [ Lat x Lon x CCFs ]
        ypred           = np.zeros(dsexp_flx.shape) * np.nan   # [ Lat x Lon x Time ]
        r2              = np.zeros((nlat,nlon)) * np.nan       # [ Lat x Lon ]
        
        # Do a silly loop (took 5 min 17 sec)
        for o in tqdm.tqdm(range(nlon)):
            lonf = lon[o]
            
            for a in range(nlat):
                latf = lat[a]
                
                chkland = proc.selpt_ds(landmask,lonf,latf).data
                if np.isnan(chkland):
                    continue
                
                
                # Check for NaN in predictor
                dspts  = [proc.selpt_ds(ds,lonf,latf) for ds in dsexp_sel]
                chknan = [np.any(np.isnan(ds.data)) for ds in dspts]
                if np.any(chknan):
                    iinan = np.where(chknan)[0][0]
                    #print("NaN detected for variables %s, lon (%.2f), lat (%.2f)... skipping." % (chknan,lonf,latf))
                    continue
                # Check for NaN in target
                flxpt  = proc.selpt_ds(dsexp_flx,lonf,latf)
                if np.any(np.isnan(flxpt.data)):
                    #print("NaN detected for Flux, lon (%.2f), lat (%.2f)... skipping." % (lonf,latf))
                    continue
                
                # Do calculations
                mlr_out = mlr_ccfs(dspts,flxpt,standardize=standardize,verbose=False)
                
                r2[a,o] = mlr_out['r2']
                ypred[a,o,:] = mlr_out['pred']
                coeffs[a,o,:] = mlr_out['coeffs']
        
        #%% Write the Output to DataArrays
         
        
        outpath         = "/home/niu4/gliu8/projects/scrap/regrid_1x1/ccfs_regression_global/"
        
        if add_ucc:
            ccfnames = ccf_vars
        else:
            ccfnames = ccf_vars[:-1]
        
        coords_r2       = dict(lat=lat,lon=lon)
        coords_coeffs   = dict(lat=lat,lon=lon,ccf=ccfnames)
        coords_pred     = dict(lat=lat,lon=lon,time=dsexp_flx.time)
        
        da_r2           = xr.DataArray(r2,coords=coords_r2,dims=coords_r2,name='r2')
        da_coeffs       = xr.DataArray(coeffs,coords=coords_coeffs,dims=coords_coeffs,name='coeffs')
        da_pred         = xr.DataArray(ypred,coords=coords_pred,dims=coords_pred,name='ypred')
        ds_out          = xr.merge([da_r2,da_coeffs,da_pred])
        edict           = proc.make_encoding_dict(ds_out)
        outname         = "%s%s_%s_CCFs_Regression_standardize%i_regrid1x1_adducc%i.nc" % (outpath,expnames[ex],flxname,standardize,add_ucc)
        ds_out.to_netcdf(outname,encoding=edict)
        
        print("Completed CCF kernel calculation for %s (%s) in %.2fs" % (flxname,expnames[ex],time.time()-st))


        }}
 

    


#%% Now Perform multiple Linear Regression (Test)

# ex      = 0
# dsflx   = dsflxs[ex]
# lonf    = 330
# latf    = 50

# flxpt = proc.selpt_ds(dsflx,lonf,latf)


# dssel = dsbyexp_anoms[ex]
# dspts = [proc.selpt_ds(ds,lonf,latf) for ds in dssel]
# #flxpt,_ = proc.match_time_month(flxpt,dspts[0])

# #% Do MLR


# std_predictors = True

# # Prepare Training Data
# # X[nsamples,nfeatures]
# # Y[nsamples]
# predictors = np.array([ds for ds in dspts]) # [variable x time]
# if std_predictors:
#     predictors = np.array([ds/np.nanstd(ds) for ds in list(predictors)])

# X          = predictors.T
# X          = np.where(np.isnan(X),0,X)
# y          = flxpt.data


# model      = LinearRegression()
# model.fit(X,y)
# #dstest = [xr.merge(ds) for ds in dsbyexp]
# pred       = model.predict(X)
# err        = y - pred # Model Error
# coeffs     = model.coef_

# r2         = sklearn.metrics.r2_score(y,pred)


#%% Now try a loop for every point for one experiment





ex             = 0


dsexp_sel      = dsbyexp_anoms[ex]
if add_ucc is False:
    dsexp_sel = dsexp_sel[:-1]
dsexp_flx      = dsflxs[ex]

# Check time dimension
ntimes_predictors = [len(ds.time) for ds in dsexp_sel]
ntimes_flux       = len(dsexp_flx.time)
if expnames[ex] == 'TCo2559-DART-1950C':
    print("Adjusting flux length")
    dsexp_flx,_ = proc.match_time_month(dsexp_flx,dsexp_sel[0])

# Pre-allocate
lon             = dsexp_flx.lon.data
lat             = dsexp_flx.lat.data
nlat,nlon,ntime = dsexp_flx.shape
nccfs           = len(dsexp_sel)
coeffs          = np.zeros((nlat,nlon,nccfs)) * np.nan # [ Lat x Lon x CCFs ]
ypred           = np.zeros(dsexp_flx.shape) * np.nan   # [ Lat x Lon x Time ]
r2              = np.zeros((nlat,nlon)) * np.nan       # [ Lat x Lon ]

# Do a silly loop (took 5 min 17 sec)
for o in tqdm.tqdm(range(nlon)):
    lonf = lon[o]
    
    for a in range(nlat):
        latf = lat[a]
        
        chkland = proc.selpt_ds(landmask,lonf,latf).data
        if np.isnan(chkland):
            continue
        
        
        # Check for NaN in predictor
        dspts  = [proc.selpt_ds(ds,lonf,latf) for ds in dsexp_sel]
        chknan = [np.any(np.isnan(ds.data)) for ds in dspts]
        if np.any(chknan):
            iinan = np.where(chknan)[0][0]
            #print("NaN detected for variables %s, lon (%.2f), lat (%.2f)... skipping." % (chknan,lonf,latf))
            continue
        # Check for NaN in target
        flxpt  = proc.selpt_ds(dsexp_flx,lonf,latf)
        if np.any(np.isnan(flxpt.data)):
            #print("NaN detected for Flux, lon (%.2f), lat (%.2f)... skipping." % (lonf,latf))
            continue
        
        # Do calculations
        mlr_out = mlr_ccfs(dspts,flxpt,standardize=standardize,verbose=False)
        
        r2[a,o] = mlr_out['r2']
        ypred[a,o,:] = mlr_out['pred']
        coeffs[a,o,:] = mlr_out['coeffs']

#%% Write the Output to DataArrays
 

outpath         = "/home/niu4/gliu8/projects/scrap/regrid_1x1/ccfs_regression_global/"

if add_ucc:
    ccfnames = ccf_vars
else:
    ccfnames = ccf_vars[:-1]

coords_r2       = dict(lat=lat,lon=lon)
coords_coeffs   = dict(lat=lat,lon=lon,ccf=ccfnames)
coords_pred     = dict(lat=lat,lon=lon,time=dsexp_flx.time)

da_r2           = xr.DataArray(r2,coords=coords_r2,dims=coords_r2,name='r2')
da_coeffs       = xr.DataArray(coeffs,coords=coords_coeffs,dims=coords_coeffs,name='coeffs')
da_pred         = xr.DataArray(ypred,coords=coords_pred,dims=coords_pred,name='ypred')
ds_out          = xr.merge([da_r2,da_coeffs,da_pred])
edict           = proc.make_encoding_dict(ds_out)
outname         = "%s%s_%s_CCFs_Regression_standardize%i_regrid1x1_adducc%i.nc" % (outpath,expnames[ex],flxname,standardize,add_ucc)
ds_out.to_netcdf(outname,encoding=edict)

#%% Check the Values


figpath = "/home/niu4/gliu8/figures/bydate/2025-11-12/"
import cartopy

def init_globalmap(nrow=1,ncol=1,figsize=(12,8)):
    proj            = ccrs.Robinson(central_longitude=-180)
    bbox            = [-180,180,-90,90]
    fig,ax          = plt.subplots(nrow,ncol,subplot_kw={'projection':proj},figsize=figsize,constrained_layout=True)
    
    multiax = True
    if (type(ax) == mpl.axes._axes.Axes) or (type(ax) == cartopy.mpl.geoaxes.GeoAxes):
        ax = [ax,]
        multiax = False
    
    for a in ax:
        a.coastlines(zorder=10,lw=0.75,transform=proj)
        a.gridlines(ls ='dotted',draw_labels=True)
        
    if multiax is False:
        ax = ax[0]
    return fig,ax

proj  = ccrs.PlateCarree()
for ii in range(len(ccfnames)):
    
    fig,ax = init_globalmap(figsize=(8,3.5))
    
    ds_out.coeffs.isel(ccf=ii).plot(ax=ax,vmin=-5,vmax=5,cmap='cmo.balance',transform=proj)
    ax.set_title("d%s / d%s" % (flxname,ccf_vars[ii]))
    savename = "%s%s_%s_CCFs_Regression_standardize%i_regrid1x1_%s_adduccc%i_check.png" % (figpath,expnames[ex],flxname,standardize,ccf_vars[ii],add_ucc)
    plt.savefig(savename,dpi=150,bbox_inches='tight')
    plt.show()



#%% Plot Mean State of Variables...

ex     = 1
varsin = dsbyexp[ex]

vlimsdict = dict(
    sst = [273,303],
    eis = [-8,8],
    Tadv = [-2.5,2.5],
    r700 = [20,70],
    WS   = [.03,.13],
    w700 = [-.1,.1],#w700 = [-65,65]   
    )
    
for vv in tqdm.tqdm(range(len(varsin))):
    
    fig,ax      = init_globalmap(figsize=(8,3.5))
    
    plotvar     = varsin[vv].mean('time')
    if plotvar.name == "Tadv":
        plotvar = plotvar * dtday
    if ex == 1 and plotvar.name == "sst":
        plotvar = plotvar + 273.15
    vlims       = vlimsdict[plotvar.name]
    pcm         = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,cmap='cmo.balance',vmin=vlims[0],vmax=vlims[1])
    cb          = fig.colorbar(pcm,ax=ax,pad=0.01,fraction=0.025)
    
    #ds_out.coeffs.isel(ccf=ii).plot(ax=ax,vmin=-5,vmax=5,cmap='cmo.balance',transform=proj)
    
    ax.set_title(plotvar.name)
    savename    = "%s%s_%s_CCFs_TimeMean_check.png" % (figpath,expnames[ex],plotvar.name)
    plt.savefig(savename,dpi=150,bbox_inches='tight')
# -----------------------------------------------------------------------------
#%% Scrap Below:
#%% Check Values Out

# SST
ds_out.coeffs.isel(ccf=0).plot(vmin=-5,vmax=5,cmap='cmo.balance'),plt.show()

# 
ds_out.coeffs.isel(ccf=1).plot(vmin=-5,vmax=5,cmap='cmo.balance'),plt.show()    

#%%

fig,ax = plt.subplots(1,1)
ax.plot(y,label='model')
ax.plot(pred,label="fit")
ax.legend()
plt.show()

#%% Compare to LSE Method

def LSE(E,Winv,y,Cnn): # Least Square Estimator from PSET1/2
    """
    Solve using Least Squares Estimator (simple inversion) given the following...
    
    inputs
    ------
        1) E : ndarray
            Observation matrix, [m x n]
        2) Winv: ndarray
            Inverse of the weights, [m x m]
        3) y : ndarray
            Observations to fit [m x 1]
        4) Cnn : ndarray
            Noise covariance of observations[m x m]
            
    outputs
    -------
        1) F@y : ndarray
            Least-squares estimated (F being the estimator) [n x 1]
        2) Cxx : ndarray
            Solution covariance [m x m]
    
    """
    
    F = np.linalg.inv(E.T@Winv@E)@E.T@Winv
    Cxx = F@Cnn@F.T
    
    return F@y,Cxx

E       = X.copy()   # [t x nvars]
y2      = y[:,None] # [t x 1]
ntime   = len(y)
Winv    = np.diag(1/ntime*np.ones((ntime)))
Cnn     = np.zeros((ntime,ntime))

coeffs_LSE,cxx = LSE(E,Winv,y2,Cnn)

print("Coefficients from Scipy")
print(coeffs)

print("Coefficients from tbx.LSE")
print(coeffs_LSE.squeeze())

# Note that they are the same for the test case
"""

Coefficients from Scipy
[-0.47701344  0.5051356   0.23814469  0.71822259  4.13316803]
Coefficients from tbx.LSE
[-0.47701344  0.5051356   0.23814469  0.71822259  4.13316803]

"""


    
    
#%% Debugging

vv          = 2

fig,ax      = init_globalmap(figsize=(8,3.5))
vlims       = [-.10,.10]

plotvar     = varsin[vv].mean('time')
#vlims       = vlimsdict[plotvar.name]
pcm         = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,cmap='cmo.balance')#,vmin=vlims[0],vmax=vlims[1])
pcm         = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,cmap='cmo.balance',vmin=vlims[0],vmax=vlims[1])
cb          = fig.colorbar(pcm,ax=ax,pad=0.01,fraction=0.025)

#ds_out.coeffs.isel(ccf=ii).plot(ax=ax,vmin=-5,vmax=5,cmap='cmo.balance',transform=proj)
ax.set_title(plotvar.name)
plt.show()
    
    


#%% Standardize Output and perform multiple linear regression




#%%

