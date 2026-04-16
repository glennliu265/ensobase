#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

# Get ENSO-related component by first regressing the CCF
based on 'get_enso_component'

Created on Wed Apr 15 15:43:29 2026

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

# =============
#%% User Edits
# =============

expname = "CERES_EBAF_ERA5_2001_2024"
flxname = "cre"
figpath = "/home/niu4/gliu8/figures/bydate/2026-04-21/"
proc.makedir(figpath)

# ===============================
#%% Additional Functions
# ===============================

# def extract_linear_component(x,y):
#     if np.any(np.isnan(x)) or np.any(np.isnan(y)): # Skip if NaN is found
#         return np.nan * np.ones(len(x))
#     # Get Fit and Calculate model (copied from detrend_poly)
#     fit     = np.polyfit(x,y,1)
#     inputs  = np.array([np.power(x,d) for d in reversed(range(len(fit)))])
#     model   = fit.T.dot(inputs)
#     return model

# def pointwise_linear_fit(ds_index,ds_target):
#     st = time.time()
#     ds_lincomp = xr.apply_ufunc(
#         extract_linear_component,
#         ds_index,
#         ds_target,
#         input_core_dims=[['time'],['time']],
#         output_core_dims=[['time']],
#         vectorize=True,
#         )
#     ds_lincomp['time'] = ds_index['time']
#     print("Extracted Component in %.2fs" % (time.time()-st))
#     return ds_lincomp

# def pointwise_r2(ds_index,ds_target):
#     st = time.time()
#     # Do separately for each point
#     calc_r2_point = lambda x,y: np.corrcoef(x,y)[0,1]**2
    
#     ds_r2   = xr.apply_ufunc(
#         calc_r2_point,
#         ds_index,
#         ds_target,
#         input_core_dims=[['time'],['time']],
#         output_core_dims=[[]],
#         vectorize=True,
#         )
    
#     print("Calculated point-wise r2 in %.2fs" % (time.time()-st))
#     return ds_r2


def get_nino_r2(nino_indices,flx_target,calc_gmean=True):
    flx_gmean    = proc.area_avg_cosweight(flx_target) 
    
    r2_flx_enso  = []
    for nn in range(2):
        nino_in                 = nino_indices[nn]
        flx_in                  = flx_target
        
        # Match Time
        flx_in,nino_in          = proc.match_time_month(flx_in,nino_in)
        
        # Get ENSO Component
        flx_ensocomp            = proc.pointwise_linear_fit(nino_in,flx_in)
        
        # Now compute local r2 for each point (to compare to ENSO component...)
        flx_raw_in,flx_ensocomp = proc.match_time_month(flx_in,flx_ensocomp)
        pointwise_r2_flx        = proc.pointwise_r2(flx_ensocomp,flx_raw_in)
        print(proc.area_avg_cosweight(pointwise_r2_flx))
        
        r2_flx_enso.append(pointwise_r2_flx)
        
    # Also compute the local r2 between flx and gmean flx
    if calc_gmean:
        r2_flx_gmean = pointwise_r2(flx_gmean,flx_target)
        return r2_flx_enso,r2_flx_gmean
    else:
        return r2_flx_enso

# ===============================
#%% Load ENSO Modes (rotated EOF)
# ===============================

st           = time.time()
if "ERA5" in expname:
    expname_enso = "ERA5_1979_2024"
else:
    expname_enso = expname


ep,cp        = ut.load_enso_eof(expname_enso,apply_smoothing=False)

proc.printtime(st)
nino_indices = [ep,cp]
nino_names   = ["EP","CP"]

# =============================
#%% Load Raw Flux
# =============================

anompath  = "/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/%s/anom_detrend1/" % expname
ncname    = "%s%s.nc" % (anompath,flxname)
ds        = xr.open_dataset(ncname)[flxname].load()
flx_raw   = ut.standardize_names(ds)

flx_gmean = proc.area_avg_cosweight(flx_raw)


# =============================
#%% Load CCF variables
# =============================

# Note: For globally averaged analysis, need to carefully consider whether or not to do land/ice masking 

ccf_vars      = ["sst","eis","Tadv","r700","w700","ws10"]
nccfs         = len(ccf_vars)

# Load CCFs
ds_ccfs       = []
for cc in tqdm(range(nccfs)):
    ccfname = ccf_vars[cc]
    ncname  = "%s%s.nc" % (anompath,ccfname)
    ds      = xr.open_dataset(ncname)[ccfname].load().squeeze()
    ds      = ut.standardize_names(ds)
    ds_ccfs.append(ds)


#ds_ccfs       = xr.merge(ds_ccfs)

#%% Compute ENSO-related component of each

nn     = 0
ninoin = nino_indices[nn]


ds_ccfs_ensocomp = []
for cc in tqdm(range(nccfs)):
    
    ccfin        = ds_ccfs[cc]
    ccfin,ninoin = proc.match_time_month(ccfin,ninoin)
    
    # Get ENSO Component
    ccf_ensocomp = proc.pointwise_linear_fit(ninoin,ccfin)
    ds_ccfs_ensocomp.append(ccf_ensocomp)
    

# --------------------------------
#%% Convolve with Radiative Kernel
# --------------------------------

# Copied loading sequence from calc_ccf_radiation_byexp.py
dskernel             = ut.load_ccf_kernel(expname,flxname,seasonal=False)
dskernel_seasonal    = ut.load_ccf_kernel(expname,flxname,seasonal=True)

ccfrad_enso          = []
ccfrad_enso_seasonal = []
for cc in tqdm(range(nccfs)):
    
    ccfname = ccf_vars[cc]
    ccf_in  = ds_ccfs_ensocomp[cc]  # ENSO Component 
    
    # Do All months
    kernel_in               = dskernel
    ccf_rad_ensocomp        = ut.convolve_kernel_ccf(ccf_in,kernel_in,ccfname,seasonal=False)
    ccfrad_enso.append(ccf_rad_ensocomp)
    
    # Do Seasonal
    kernel_in               = dskernel_seasonal
    ccf_rad_ensocomp_seas   = ut.convolve_kernel_ccf(ccf_in,kernel_in,ccfname,seasonal=True)
    ccfrad_enso_seasonal.append(ccf_rad_ensocomp_seas)
    
#%% Now calculate r2 with Global Mean cre

r2_byccf          = []
r2_byccf_seasonal = []
for cc in tqdm(range(nccfs)):
    
    # All Months
    flx_in          = ccfrad_enso[cc]
    flx_in,gmean_in = proc.match_time_month(flx_in,flx_gmean)
    r2all           = pointwise_r2(gmean_in,flx_in)
    r2_byccf.append(r2all)
    
    # Seasonal
    flx_in          = ccfrad_enso_seasonal[cc]
    flx_in,gmean_in = proc.match_time_month(flx_in,flx_gmean)
    r2all_seasonal  = pointwise_r2(gmean_in,flx_in)
    r2_byccf_seasonal.append(r2all_seasonal)
    

#%% Also Calculate for the summed values

dssum           = xr.concat(ccfrad_enso,dim='ccf').sum('ccf')
flx_in,gmean_in = proc.match_time_month(dssum,flx_gmean)
r2_ccfsum       = pointwise_r2(gmean_in,flx_in)

dssum_seasonal  = xr.concat(ccfrad_enso_seasonal,dim='ccf').sum('ccf')
flx_in,gmean_in = proc.match_time_month(dssum,flx_gmean)
r2_ccfsum       = pointwise_r2(gmean_in,flx_in)

# Note: I Hit this point and realized that the values are the same everywhere, since each component
# has correlation with ENSO index of 1, the value is just correlation between ENSO index and Gmean Flux

#%%

# =============================================================================
# =============================================================================
#%% Scrap Below

#%%

            
        

    nccf  = len(ccf_vars)
    dtday = 3600*60
    for vv in tqdm(range(nccf)):
        
        # Get Variable
        ccfname        = ccf_vars[vv]
        varanom        = dsvars_anom[vv]
        
        varanom_std     = varanom / varanom.std('time') # Standardize Variable
        
        # Multiple by the Coefficient
        coeff_allmons  = dsall.coeffs.sel(ccf=ccfname) #/dtday
        R_component    = coeff_allmons * varanom_std
        vname_new      = ccfname
        R_component    = R_component.rename(vname_new)
        if customname is not None: # Add after Flux
            rname_out      = "%s%s_%s_%s_component.nc" % (outpath,flxname,customname,vname_new,)
        else:
            rname_out      = "%s%s_%s_component.nc" % (outpath,flxname,vname_new,)
        R_component.to_netcdf(rname_out)
        
        if calc_seasonal:
            coeff_seasonal = [ds.coeffs.sel(ccf=ccfname) for ds in ds_byseason] # /dtday
            R_component_seasonal = []
            for ss in range(4):
                
                selmons    = selmons_loop[ss]
                
                varmon     = proc.selmon_ds(varanom_std,selmons)
                varmon_out = varmon * coeff_seasonal[ss]
                R_component_seasonal.append(varmon_out)
                
            R_component_seasonal = xr.concat(R_component_seasonal,dim='time')
            R_component_seasonal = R_component_seasonal.sortby('time')
            R_component_seasonal = R_component_seasonal.rename(vname_new)
            R_component_seasonal = R_component_seasonal#.rename({'time':'valid_time'})
            if customname is not None:
                rname_out      = "%s%s_%s_%s_component_seasonal.nc" % (outpath,flxname,customname,vname_new,)
            else:
                rname_out      = "%s%s_%s_component_seasonal.nc" % (outpath,flxname,vname_new,)
            R_component_seasonal.to_netcdf(rname_out)
    



#%% Load CCF components (copied from check_ccf_seasonality_sep_nep.ipynb)

# Load CCF Radiation Components
st = time.time()
ccf_fluxes  = ut.load_ccf_radiation(expname,flxname)
ccf_fluxes_seas   = ut.load_ccf_radiation(expname,flxname,seasonal=True)
proc.printtime(st,"Loaded CCF Components")

#%% Load Raw Flux Anomaly

anompath  = "/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/%s/anom_detrend1/" % expname
ncname    = "%s%s.nc" % (anompath,flxname)
ds        = xr.open_dataset(ncname)[flxname].load()
flx_raw   = ut.standardize_names(ds)

flx_gmean = proc.area_avg_cosweight(flx_raw)

# ====================================================
#%% Check Seasonal versus All Month reconstructed flux 
# ====================================================

ccf_fluxes          = [ccf.squeeze() for ccf in ccf_fluxes]
ccf_fluxes_sum      = xr.concat(ccf_fluxes,dim='ccf')
ccf_fluxes_sum      = ccf_fluxes_sum.sum('ccf')

ccf_fluxes_seas     = [ccf.squeeze() for ccf in ccf_fluxes_seas]
ccf_fluxes_seas_sum = xr.concat(ccf_fluxes_seas,dim='ccf')
ccf_fluxes_sum_seas = ccf_fluxes_seas_sum.sum('ccf')

residual_allmon = flx_raw - ccf_fluxes_sum
residual_seas   = flx_raw - ccf_fluxes_sum_seas
    
# ================================================
#%% Part (1) Get ENSO-related component (Raw Flux)
# ================================================

nn                      = 1

for nn in range(2):
    
    # Loopy Stuff
    nino_in                 = nino_indices[nn]
    flx_in                  = flx_raw
    
    # ----- This should be within loop
    
    # Match Time
    flx_in,nino_in          = proc.match_time_month(flx_in,nino_in)
    
    # Get ENSO Component
    flx_ensocomp            = pointwise_linear_fit(nino_in,flx_in)
    
    # Now compute local r2 for each point (to compare to ENSO component...)
    flx_raw_in,flx_ensocomp = proc.match_time_month(flx_raw,flx_ensocomp)
    pointwise_r2_flx        = pointwise_r2(flx_ensocomp,flx_raw_in)
    print(proc.area_avg_cosweight(pointwise_r2_flx))
    
    # Also compute the local r2 between flx and gmean flx
    pointwise_r2_global_mean = pointwise_r2(flx_gmean,flx_raw_in)
    
    #%% Make Some Plots
    
    #%% R2 Between Flx and ENSO-related component at each point
    
    proj    = ccrs.PlateCarree()
    #cints   = np.arange(-1,1.1,0.1)
    cints   = np.arange(-0.40,0.45,0.05)
    fig,ax  = ut.init_globalmap(1,1)
    
    plotvar = pointwise_r2_flx
    pcm     = ax.contourf(plotvar.lon,plotvar.lat,plotvar,transform=proj,
                          levels=cints,cmap='cmo.balance')
    cb      = viz.hcbar(pcm)
    ax.set_title("Pointwise $r^2$, %s-ENSO and %s" % (nino_names[nn],flxname))
    #ax.coastlines()
    
    figname = "%sR2_Pointwise_%s_%s_vs_%s-ENSO.png" % (figpath,expname,flxname,nino_names[nn])
    plt.savefig(figname,dpi=150,bbox_inches='tight')
    
    #%% R2 between Flx and Global Mean at each point
    fig,ax  = ut.init_globalmap(1,1)
    cints   = np.arange(-0.30,0.33,0.03)
    plotvar = pointwise_r2_global_mean
    pcm     = ax.contourf(plotvar.lon,plotvar.lat,plotvar,transform=proj,
                          levels=cints,cmap='cmo.balance')
    cb      = viz.hcbar(pcm)
    ax.set_title("Pointwise $r^2$, Global mean and local %s" % (flxname))
    
    figname = "%sR2_Pointwise_%s_%s_GlobalMean_vs_Local.png" % (figpath,expname,flxname)
    plt.savefig(figname,dpi=150,bbox_inches='tight')
    
    plt.show()

#%% Check Both

    
r2_flx_ccfsum,r2_gmean_ccfsum           = get_nino_r2(nino_indices,ccf_fluxes_sum)
r2_flx_ccfsum_seas,r2_gmean_ccfsum_seas = get_nino_r2(nino_indices,ccf_fluxes_sum_seas)


#%% Plot

r2_in = [r2_flx_ccfsum,r2_flx_ccfsum_seas]
r2names = ["All Months","Seasonal"]

fig,axs  = ut.init_globalmap(2,2,figsize=(20,14))
axs = axs.reshape((2,2))
for nn in range(2):
    for ss in range(2):
        
        ax = axs[nn,ss]
        if nn == 0:
            ax.set_title(r2names[ss])
        elif ss == 0:
            viz.add_ylabel(nino_names[nn],ax=ax)
        
        plotvar=r2_in[ss][nn]
        
        pcm     = ax.contourf(plotvar.lon,plotvar.lat,plotvar,transform=proj,
                              levels=cints,cmap='cmo.balance')


cb      = viz.hcbar(pcm,ax=axs.flatten())

figname = "%sR2_Pointwise_%s_%s_CCFReconstruction_Seasonality_Effect.png" % (figpath,expname,flxname)
plt.savefig(figname,dpi=150,bbox_inches='tight')

plt.show()

# ===========================
#%% Now Compute for each CCF
# ===========================
# i.e. how much does EP/CP Index explain for local radiation variability (by CCF)

r2_byccf      = []
r2_byccf_seas = []
for cc in tqdm(range(6)):
    
    r2out,_      = get_nino_r2(nino_indices,ccf_fluxes[cc])
    r2_byccf.append(r2out)
    
    r2out_seas,_ = get_nino_r2(nino_indices,ccf_fluxes_seas[cc])
    r2_byccf_seas.append(r2out_seas)


#%% Visualize by CCF

cints         = np.arange(-0.40,0.45,0.05)
proj          = ccrs.PlateCarree()
ccf_vars      = ["sst","eis","Tadv","r700","w700","ws10"]

ccfcalc_names = ["Allmon","Seasonal"]
r2_byccf_all  = [r2_byccf,r2_byccf_seas]

for ic in range(2):
    r2_byccf_in =r2_byccf_all[ic]
    
    fig,axs  = ut.init_globalmap(2,6,figsize=(28,8))
    axs      = axs.reshape(2,6)
    
    for nn in range(2):
        
        for cc in range(6):
            
            
            ax = axs[nn,cc]
            if nn == 0:
                ax.set_title(ccf_vars[cc])
            if cc == 0:
                viz.add_ylabel(nino_names[nn],ax=ax)
            
            
            plotvar = r2_byccf_in[cc][nn]
            pcm     = ax.contourf(plotvar.lon,plotvar.lat,plotvar,transform=proj,
                                  levels=cints,cmap='cmo.balance')
    
    #cb = fig.colorbar(pcm,ax=axs.flatten(),pad=0.02,fraction=0.01)
    
    figname = "%sR2_Pointwise_%s_%s_vs_ENSO_ccfcalc_%s.png" % (figpath,expname,flxname,ccfcalc_names[ic])
    plt.savefig(figname,dpi=150,bbox_inches='tight')

# ========================================================================================
#%% Do Some Seasonal Analysis (1): EP and CP ENSO Correlations with global mean components
# ========================================================================================
season_indices = [[12,1,2],[3,4,5],[6,7,8],[9,10,11]]
flxcalcnames   = ["FlxRaw","AllMonRecon","SeasRecon","AllMonRes","SeasRes"]
flxes_in       = [flx_raw,ccf_fluxes_sum,ccf_fluxes_sum_seas,residual_allmon,residual_seas]
season_names   = ['All','DJF','MAM','JJA','SON']

corr_byflx = []
for ff in range(5):
    
    flx_in       = flxes_in[ff]
    flxcalcname  = flxcalcnames[ff]#"FlxRaw"
    flx_in_gmean = proc.area_avg_cosweight(flx_in)
    
    corr_bynino  = []
    for nn in range(2):
        
        corr_byseason        = []
        nino_in              = nino_indices[nn]
        
        nino_in,flx_in_gmean = proc.match_time_month(nino_in,flx_in_gmean)
        corr_byseason.append(np.corrcoef(nino_in,flx_in_gmean)[0,1]**2)
        
        for ss in range(4):
            
            sid       = season_indices[ss]
            nino_seas = proc.selmon_ds(nino_in,sid)
            flx_seas  = proc.selmon_ds(flx_in_gmean,sid)
            
            corrseas  = np.corrcoef(nino_seas,flx_seas)[0,1]**2
            
            corr_byseason.append(corrseas)
        
        corr_bynino.append(corr_byseason)
    
    corr_bynino = np.array(corr_bynino)
    corr_byflx.append(corr_bynino)
    
    # # Make Plot
    # fig,ax = plt.subplots(1,1,figsize=(10,4.5),constrained_layout=True)
    # pcm    = ax.pcolormesh(season_names,nino_names,corr_bynino)
    # ax.set_aspect('equal')
    # #ax.grid(True,ls='solid',lw=0.55,c="k")
    # fig.colorbar(pcm,ax=ax,fraction=0.01,pad=0.05)
    # plt.show()
    
    # # Make Plot
    # fig,ax = plt.subplots(1,1,figsize=(10,4.5),constrained_layout=True)
    # pcm    = ax.pcolormesh(season_names,nino_names,corr_bynino)
    # ax.set_aspect('equal')
    # #ax.grid(True,ls='solid',lw=0.55,c="k")
    # fig.colorbar(pcm,ax=ax,fraction=0.01,pad=0.05)
    # plt.show()
    
    #%%
    
    # Bar PLot in the style of Moon et al. 2025...
    
    # Plotting Options
    x            = np.arange(5)  # the label locations
    width        = 0.25           # the width of the bars
    label_bars   = True
    
    # Other Edits
    nino_colors  = ['salmon','cornflowerblue']
    
    nset   = 1
    # Initialize Figure
    fig,ax    = plt.subplots(1,1,layout='constrained',figsize=(8,4.5))
    
    multiplier = 0
    for ss in range(5):
        offset = width * multiplier
    
        # Loop for each Experiment
        for ii in range(2):
                
            offset      = width * ii
            color_in    = nino_colors[ii]
            
            if ss == 0: # Label onlyon first iteration
                label=nino_names[ii]
            else:
                label=""
            
            # Get Data and PLot
            measurement = corr_bynino[ii,:]
            rects       = ax.bar(x + offset, measurement, width, label=label,
                            color=color_in,linewidth=0.5)
            
            # Set Axis Ticks and option to label bars
            ax.set_xticks(x,season_names)
            if label_bars:
                ax.bar_label(rects, padding=3,fmt="%.02f")
            multiplier += 1
    
        # Axis Labeling and Legend
        ax.set_ylabel(r"$r^2$(ENSO Index,Gmean TOA %s)" % (flxname))
        ax.set_xlabel("Months used in Calculation")
        ax.set_ylim([0,0.60])
        ax.legend()
    
    figname = "%s%s_R2_TOA_%s_Globalmean_vs_ENSO_Indices_%s.png" % (figpath,expname,flxname,flxcalcname)
    plt.savefig(figname,dpi=150,bbox_inches='tight')
    plt.close()
    #plt.show()

# ====================================================
#%% Do Some Seasonal Analysis (2): Pointwise Analysis
# ====================================================

flx=flx_raw

r2_byseason = [] # [Season][Nino][Lat x Lon]
for ss in tqdm(range(4)):
    
    sid          = season_indices[ss]
    flx_in_seas  = proc.selmon_ds(flx_raw,sid)
    nino_in_seas = [proc.selmon_ds(nn,sid) for nn in nino_indices]
    r2out,_      = get_nino_r2(nino_in_seas,flx_in_seas)
    
    r2_byseason.append(r2out)


#%%

cints   = np.arange(-0.40,0.45,0.05)
proj    = ccrs.PlateCarree()


fig,axs = ut.init_globalmap(2,4,figsize=(28,8))
axs     = axs.reshape(2,4)

for nn in range(2):
    for ss in range(4):
        
        ax = axs[nn,ss]
        if nn == 0:
            ax.set_title(season_names[ss+1])
        if ss == 0:
            viz.add_ylabel(nino_names[nn],ax=ax)
        
        plotvar = r2_byseason[ss][nn]
        pcm     = ax.contourf(plotvar.lon,plotvar.lat,plotvar,transform=proj,
                              levels=cints,cmap='cmo.balance')
        
        
            

#cb = fig.colorbar(pcm,ax=axs.flatten(),pad=0.02,fraction=0.01)

figname = "%sR2_Pointwise_%s_%s_vs_ENSO_seasonal.png" % (figpath,expname,flxname)
plt.savefig(figname,dpi=150,bbox_inches='tight')


plt.show()


# =============================================================================
#%% Update 2026.04.14 : Try Getting ENSO-related component of a CCF, then
# Convolving with the kernel
# =============================================================================





#%%



#%% Scrap/Old Script Below

#fig,ax = plt.subplots(2,1)


#%%


# ---------


## Calculate Global Mean (note this just returns the correlation between flx_gmean and ep/cp enso)
gmean_in,flx_ensocomp   = proc.match_time_month(flx_gmean,flx_ensocomp)
r2_flx                  = pointwise_r2(flx_gmean,flx_ensocomp,)


# -----------------------------------------------------------------------------
# Debug Block for Linear Component Extraction
x = nino_in
y = proc.selpt_ds(flx_in,240,0)
fit     = np.polyfit(x,y,1)
inputs  = np.array([np.power(x,d) for d in reversed(range(len(fit)))])
model   = fit.T.dot(inputs)
plt.scatter(x,y),plt.plot(x,model),plt.show()
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Debug block for R2 calculation
x = gmean_in
y = proc.selpt_ds(flx_ensocomp,240,0)
calc_r2_point = lambda x,y: np.corrcoef(x,y)[0,1]**2
print(calc_r2_point(x,y))
# -----------------------------------------------------------------------------



    
    




    
#%% First Part: Compute for the full fluxes












