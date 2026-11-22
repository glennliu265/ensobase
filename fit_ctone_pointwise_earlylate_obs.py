#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Script version of find_earlylate_split_ssp585.ipynb
Runs on niu

Performs for early period (2 full awi-cm3 ssp585 runs)
and for late period (2 full and 3 ensemble spinoff)


Created on Fri Nov 20 13:48:13 2026

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
import scipy as sp
import os

from tqdm import tqdm

#%% Import Custom Modules
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

tbxpath = "/home/niu4/gliu8/scripts/commons/tbx"
sys.path.append(tbxpath)
import tbx as tbx

import xskillscore as xs
import cartopy.feature as cfeature

#%% Helper Functions

# Anomalize and Detrend
def preproc(ds):
    dsanom = proc.xrdeseason(ds)
    dsanom = proc.xrdetrend_dim(dsanom,dim='time',deg=2)
    return dsanom

def preproc_dataset(ds,varlist):
    dsanoms = []
    for vname in varlist:
        dsa = preproc(ds[vname])
        dsanoms.append(dsa.rename(vname))
    return xr.merge(dsanoms)

def cmode_fit_loop_lpf(varanom_in,ninoid_in,lpfmon,periodname,
                       poster_path,expname,vname,
                       ):
    """
    Loop Version of below. Fit Cmode Model Hierarchy and Save Output
    
    Inputs
        varanom_in  (xr.DataArray) --> Variable to Fit
        ninoid_in   (xr.DataArray) --> Nino Index
        lpfmon      (INT) --> Number of months to low pass filter
        periodname  (STR) --> Name of the period
        poster_path (STR) --> Output Path
        expname     (STR) --> Name of Experiment
        vname       (STR) --> Variable Name
    
    """
    fitnames         = ["NinoOnly","CmodeOnly","CmodeNino"]
    
    # Apply Low Pass Filter
    st               = time.time()
    varanom_in_lpf   = proc.pointwise_lp(varanom_in,lpfmon)
    print("Applied LPF in %.fs" % (time.time()-st))
    
    # Fit Models
    models           = fit_cmode_enso_hierarchy(varanom_in,ninoid_in)
    models_lpf       = fit_cmode_enso_hierarchy(varanom_in_lpf,ninoid_in)
    
    # Save Output 
    st = time.time()
    for ff in tqdm(range(3)):
        fitname      = fitnames[ff]
        
        # Save Regular
        ds_out       = models[fitname]
        outname      = "%s%s_fit_%s_%s_%s.nc" % (poster_path,vname,expname,periodname,fitname) # "cre_fit_NinoOnly
        ds_out.to_netcdf(outname)
        
        # Save LPF
        ds_out_lpf   = models_lpf[fitname]
        outname_lpf  = proc.addstrtoext(outname,"_LPF%02i" % lpfmon)
        ds_out_lpf.to_netcdf(outname_lpf)
        
    # Also Save Variables (Original and LPF), and Nino3.4
    # Full Anomaly
    outname      = "%s%s_anomaly_%s_%s.nc" % (poster_path,vname,expname,periodname)
    varanom_in.to_netcdf(outname)
    # LPF Anomaly
    outname_lpf  = proc.addstrtoext(outname,"_LPF%02i" % lpfmon)
    varanom_in_lpf.to_netcdf(outname_lpf)
    
    # Nino3.4 Index
    
    outname_nino =  "%sNino34_%s_%s.nc" % (poster_path,expname,periodname)
    if not os.path.exists(outname_nino):
        ninoid_in.to_netcdf(outname_nino)
    else:
        print("Not Saving Nino3.4 as it has already exists.")
    print("Period Output for %s saved in %.fs" % (expname,time.time()-st))
    
    return varanom_in_lpf,models,models_lpf

def fit_cmode_enso_hierarchy(v_in,ninoid_in):
    # note: Time Dimensions for varanom_in and ninoid_in should match
    # outputs dict with
    # "NinoOnly" Fit to Linear ENSO
    # "CmodeOnly" Fit to Cmode Only
    # "CmodeNino" Fit to Cmode + Linear ENSO
    
    st               = time.time()

    # Fit for ENSO Only
    ninofit_var      = proc.xr_linfit(v_in,ninoid_in,deg=1) 
    # Fit for Cmode + ENSO Only
    ctonefit_var     = ut.fit_ctone_pointwise(v_in,ninoid_in,linear_enso=True,fix_frequency_offset=True)
    # Fit for Cmode Only
    ctoneonlyfit_var = ut.fit_ctone_pointwise(v_in,ninoid_in,linear_enso=False,fix_frequency_offset=True)
    print("Completed All Fits in %.2fs" % (time.time()-st))
    
    models = {
        "NinoOnly"  : ninofit_var,
        "CmodeOnly" : ctoneonlyfit_var,
        "CmodeNino" : ctonefit_var,
        }
    return models

def fit_cmode_cp_ep_hierarchy(v_in,cpep_in):
    # EP.CP Version of [fit_cmode_enso_hierarchy]
    # note: Time Dimensions for varanom_in and ninoid_in should match
    # outputs dict with
    # "NinoOnly" Fit to Linear ENSO
    # "CmodeOnly" Fit to Cmode Only
    # "CmodeNino" Fit to Cmode + Linear ENSO
    
    st               = time.time()
    
    # Set Up
    cp_in           = cpep_in.cp_index
    ep_in           = cpep_in.ep_index
    
    # Fit for CP/EP Indices (MLR)
    predictors       = xr.concat([cp_in,ep_in],dim='predictors').transpose('time','predictors')
    ninofit_var      = proc.pointwise_mlr(predictors,v_in,standardize=False,predictor_names=["cp","ep"],return_pred=True)
    
    # Fit for Cmode + ENSO 
    ctonefit_var     = ut.fit_ctone_cp_ep_pointwise(v_in,cp_in,ep_in,return_ds=True,
                        linear_enso=True)
    
    # Fit for Cmode Only
    ctoneonlyfit_var = ut.fit_ctone_cp_ep_pointwise(v_in,cp_in,ep_in,return_ds=True,
                        linear_enso=False)
    print("Completed All Fits in %.2fs" % (time.time()-st))
    
    models = {
        "NinoOnly"  : ninofit_var,
        "CmodeOnly" : ctoneonlyfit_var,
        "CmodeNino" : ctonefit_var,
        }
    return models

def cmode_fit_loop_lpf_cpep(varanom_in,cpep_in,lpfmon,periodname,
                       poster_path,expname,vname,
                       ):
    """
    CP/EP Version of [cmode_fit_loop_lpf]
    Loop Version of below. Fit Cmode Model Hierarchy and Save Output
    
    Inputs
        varanom_in  (xr.DataArray) --> Variable to Fit
        cpep_in   (xr.DataArray) --> CP/EP Index
        lpfmon      (INT) --> Number of months to low pass filter
        periodname  (STR) --> Name of the period
        poster_path (STR) --> Output Path
        expname     (STR) --> Name of Experiment
        vname       (STR) --> Variable Name
    
    """
    fitnames         = ["NinoOnly","CmodeOnly","CmodeNino"]
    
    # Apply Low Pass Filter
    st               = time.time()
    varanom_in_lpf   = proc.pointwise_lp(varanom_in,lpfmon)
    print("Applied LPF in %.fs" % (time.time()-st))
    
    # Fit Models
    models           = fit_cmode_cp_ep_hierarchy(varanom_in,cpep_in)
    models_lpf       = fit_cmode_cp_ep_hierarchy(varanom_in_lpf,cpep_in)
    
    # Save Output 
    st = time.time()
    for ff in tqdm(range(3)):
        fitname      = fitnames[ff]
        
        # Save Regular
        ds_out       = models[fitname]
        outname      = "%s%s_fit_CPEP_%s_%s_%s.nc" % (poster_path,vname,expname,periodname,fitname) # "cre_fit_NinoOnly
        ds_out.to_netcdf(outname)
        
        # Save LPF
        ds_out_lpf   = models_lpf[fitname]
        outname_lpf  = proc.addstrtoext(outname,"_LPF%02i" % lpfmon)
        ds_out_lpf.to_netcdf(outname_lpf)
        
    # Also Save Variables (Original and LPF), and Nino3.4
    # Full Anomaly
    outname      = "%s%s_anomaly_%s_%s.nc" % (poster_path,vname,expname,periodname)
    if not os.path.exists(outname):
        varanom_in.to_netcdf(outname)
    else:
        print("Not Saving Full Anomaly as it has already exists.")
    # LPF Anomaly
    outname_lpf  = proc.addstrtoext(outname,"_LPF%02i" % lpfmon)
    if not os.path.exists(outname_lpf):  
        varanom_in_lpf.to_netcdf(outname_lpf)
        print("Not Saving LPF Anomaly as it has already exists.")
    
    # CP/EP Index
    outname_nino =  "%sCPEP_Index_%s_%s.nc" % (poster_path,expname,periodname)
    if not os.path.exists(outname_nino):
        cpep_in.to_netcdf(outname_nino)
    else:
        print("Not Saving CP.EP Indices as it has already exists.")
    print("Period Output for %s saved in %.fs" % (expname,time.time()-st))
    
    return varanom_in_lpf,models,models_lpf

#%% Set some Paths and Variables

figpath = "/home/niu4/gliu8/figures/bydate/2026-11-23/"
proc.makedir(figpath)
proj    = ccrs.PlateCarree()

# Output Path will Be Here
poster_path = "/home/niu4/gliu8/projects/ccfs/enso_ctone_fits/fit_ctone_pointwise/"

# Set Variable Properties
vnames    = ["cre","tscre","ttcre","sst"]
expnames  = ["TCo319-DART-ssp585d-gibbs-charn",
             "TCo319_ssp585",
             "TCo319_ssp585_ens01",
             "TCo319_ssp585_ens02",
             "TCo319_ssp585_ens03",
             "CERES",]

expcolors = [
    "blue",
    "lightsalmon",
    "orange",
    "darkviolet",
    "gold",
    "k"
    ]

flxnames = ["cre","tscre","ttcre"]

lpfmon         = 6  

"""

--------------------
Two Different Modes
--------------------
nino34 : Fit to Niño3.4 Index
cpep   : Fit to CP/EP Index

"""
mode  = "cpep" # see above

nexps = len(expnames)
nvars = len(vnames)

#%%  Load Variables (AWI-CM3)

# Load Variables
varsbyexp = []
for ex in tqdm(range(nexps-1)):
    expname = expnames[ex]

    vbv = []
    for vv in range(nvars):
        vname = vnames[vv]
        ds    = ut.loadregrid(expname,vname,reformat=True)
        vbv.append(ds)
    varsbyexp.append(vbv)
    
# Merge Dataset for AWI
varsawi = [xr.merge(ds) for ds in varsbyexp]

#%% Load for Observations

varsobs = []
for vv in range(nvars):
    vname = vnames[vv]
    if vname == "sst":
        expin = "ERA5"
    else:
        expin = "CERES_EBAF"
        print("Loading %s for %s" % (expin,vname))
    ds    = ut.loadregrid(expin,vname,reformat=True)

    #if expin == "ERA5":
    ds = ds.sel(time=slice('2000-03-01','2024-12-31'))
    
    varsobs.append(ds)
    print(ds.shape)

varsobs        = xr.merge(varsobs)

#%% Determine Early_Late Split

# Indicate Periods
trange_early   = ["2015-01-01","2054-12-31"]
trange_late    = ["2060-01-01","2099-12-31"]

# Select Periods
vars_awi_early = [ds.sel(time=slice(*trange_early)) for ds in varsawi[:2]] # Only the first 2 runs
vars_awi_late  = [ds.sel(time=slice(*trange_late)) for ds in varsawi]

# Make String for Saving Later
early_str      = "early_%sto%s" % (trange_early[0][:4],trange_early[1][:4])
late_str       = "late_%sto%s" % (trange_late[0][:4],trange_late[1][:4])

# ================================
#%% Preprocess (AWI-CM3 Runs)
# ================================

awi_early_anom = [preproc_dataset(ds,vnames) for ds in vars_awi_early]
awi_late_anom  = [preproc_dataset(ds,vnames) for ds in vars_awi_late]

#%% Calculate ENSO ID for each case

if mode == "nino34":
    
    print("Calculating and fitting to the Niño3.4 Index")
    
    bbox_nino34       = [-170+360,-120+360,-5,5]     # Nino3.4 Box
    ninoid_early      = [proc.aavg(ds.sst,bbox_nino34) for ds in awi_early_anom]
    ninoid_late       = [proc.aavg(ds.sst,bbox_nino34) for ds in awi_late_anom]
    
    
    # ===============================
    # %% Perform Fitting (AWI-CM3)
    # ===============================
    
    """
    
    Looping for each variable 
    Looping for each experiment
    Looping for Early/Late
    
    - Get Anomalies
    - Get Nino3.4
    - LPF Anomalies
    
    - Loop for LPF and Full Anomalies...
        - Do Pointwise Fit (ENSO Only)
        - Do Pointwise Fit (ENSO + Cmode)
        - Do Pointwise Fit (Cmode Only)
    
    """
    
    
    
    #%% Loop Through Early/Late Periods for AWI-CM3
    
    #fitnames       = ["NinoOnly","CmodeOnly","CmodeNino"]
    
    
    for vname in flxnames: # "tscre",
        #vname          = "cre"
        for ex in tqdm(np.arange(5)): # range(5)
            print("Starting fit for: %s" % expnames[ex])
            
            # Process Late Period First (Valid for all experiments) --------
            varanom_in        = awi_late_anom[ex][vname]
            ninoid_in         = ninoid_late[ex]
            periodname        = late_str
            
            # Fit Models and Save
            varanom_in_lfp,models,models_lpf = cmode_fit_loop_lpf(varanom_in,ninoid_in,lpfmon,periodname,
                                   poster_path,expnames[ex],vname,)
            
            if ex < 2:
                print("\tAlso Performing Early Calculations for %s" % (expnames[ex]))
                # Set Parameters
                varanom_in       = awi_early_anom[ex][vname]
                ninoid_in        = ninoid_early[ex]
                periodname       = early_str
                
                # Fit Models and Save
                varanom_in_lfp,models,models_lpf = cmode_fit_loop_lpf(varanom_in,ninoid_in,lpfmon,periodname,
                                       poster_path,expnames[ex],vname,)
                # End Early Period Calculation
            # End Experiment Loop
                
        #%% Now Repeat for Observations
        
        # Take Anomalies
        varsobs_anom = preproc_dataset(varsobs,vnames)
        
        # Compute Nino3.4
        ninoid_obs   = proc.aavg(varsobs_anom.sst,bbox_nino34)
        
        
        # Perform Analysis
        varanom_in        = varsobs_anom[vname]
        ninoid_in         = ninoid_obs
        expname           = "Obs"
        periodname        = "%sto%s" % (
            str(varsobs_anom.time[0].data)[:4],
            str(varsobs_anom.time[-1].data)[:4],
            )
        
        # Fit Models and Save
        _,_,_ = cmode_fit_loop_lpf(varanom_in,ninoid_in,lpfmon,periodname,
                               poster_path,expname,vname,)
     
elif mode == "cpep":
    print("Calculating and fitting to the CP/EP Index")
    
    # First, Calculate CP and EP Index
    cpep_early       = [ut.calc_cp_ep_index(ds) for ds in awi_early_anom]
    cpep_late        = [ut.calc_cp_ep_index(ds) for ds in awi_late_anom]
    
    
    
    # Loop by Variable
    for vname in flxnames:  # ["cre",]:#
        for ex in tqdm(np.arange(5)): # range(5)
            print("Starting fit (CP/EP) for: %s" % expnames[ex])
            
            # Process Late Period First (Valid for all experiments) --------
            varanom_in        = awi_late_anom[ex][vname]
            cpep_in           = cpep_late[ex]
            periodname        = late_str
            
            # Fit Models and Save
            varanom_in_lfp,models,models_lpf = cmode_fit_loop_lpf_cpep(varanom_in,cpep_in,lpfmon,periodname,
                                   poster_path,expnames[ex],vname,)
            
            if ex < 2:
                print("\tAlso Performing Early Calculations for %s" % (expnames[ex]))
                # Set Parameters
                varanom_in       = awi_early_anom[ex][vname]
                cpep_in        = cpep_early[ex]
                periodname       = early_str
                
                # Fit Models and Save
                varanom_in_lfp,models,models_lpf = cmode_fit_loop_lpf_cpep(varanom_in,cpep_in,lpfmon,periodname,
                                       poster_path,expnames[ex],vname,)
                # End Early Period Calculation
            # End Experiment Loop
                
        #%% Now Repeat for Observations
        
        # Take Anomalies
        varsobs_anom = preproc_dataset(varsobs,vnames)
        
        # Compute CP/EP Index
        cpep_obs     = ut.calc_cp_ep_index(varsobs_anom)
        
        
        # Perform Analysis
        varanom_in        = varsobs_anom[vname]
        cpep_in         = cpep_obs
        expname           = "Obs"
        periodname        = "%sto%s" % (
            str(varsobs_anom.time[0].data)[:4],
            str(varsobs_anom.time[-1].data)[:4],
            )
        
        # Fit Models and Save
        _,_,_ = cmode_fit_loop_lpf_cpep(varanom_in,cpep_in,lpfmon,periodname,
                               poster_path,expname,vname,)
     
        
    
    # # Troubleshooting, test for pointwise fit
    
    # target = awi_early_anom[0].sst
    # cp_in      = cpep_early[0].cp_index
    # ep_in      = cpep_early[0].ep_index
    
    # predictors = xr.concat([cp_in,ep_in],dim='predictors').transpose('time','predictors')
    # mlr_out    = proc.pointwise_mlr(predictors,target,standardize=False,predictor_names=["cp","ep"],return_pred=True)
    
   
    



#%%
# #%% Scrap Below...
    
    
    


# #%%






# #ex               = 0
# for ex in np.arange(5): # range(5)
    
#     # Process Late Period First (Valid for all experiments) --------
#     varanom_in       = awi_late_anom[ex][vname]
#     ninoid_in        = ninoid_late[ex]
#     periodname       = late_str
    
#     # Apply Low Pass Filter
#     st               = time.time()
#     varanom_in_lpf   = proc.pointwise_lp(varanom_in,lpfmon)
#     print("Applied LPF in %.fs" % (time.time()-st))
    
#     models           = fit_cmode_enso_hierarchy(varanom_in,ninoid_in)
#     models_lpf       = fit_cmode_enso_hierarchy(varanom_in_lpf,ninoid_in)
        
#     # Save Output 
#     st = time.time()
#     for ff in tqdm(range(3)):
#         fitname      = fitnames[ff]
        
#         # Save Regular
#         ds_out       = models[fitname]
#         outname      = "%s%s_fit_%s_%s_%s.nc" % (poster_path,vname,expnames[ex],periodname,fitname) # "cre_fit_NinoOnly
#         ds_out.to_netcdf(outname)
        
#         # Save LPF
#         ds_out_lpf   = models_lpf[fitname]
#         outname_lpf  = proc.addstrtoext(outname,"_LPF%02i" % lpfmon)
#         ds_out_lpf.to_netcdf(outname_lpf)
        
#         # End Fit Loop
    
    
#     # Also Save Variables (Original and LPF), and Nino3.4
#     # Full Anomaly
#     outname      = "%s%s_anomaly_%s_%s.nc" % (poster_path,vname,expnames[ex],periodname)
#     varanom_in.to_netcdf(outname)
#     # LPF Anomaly
#     outname_lpf  = proc.addstrtoext(outname,"_LPF%02i" % lpfmon)
#     varanom_in_lpf.to_netcdf(outname_lpf)
#     # Nino3.4 Index
#     outname_nino =  "%sNino34_%s_%s.nc" % (poster_path,expnames[ex],periodname)
#     ninoid_in.to_netcdf(outname_nino)
#     print("Late Period Output for %s saved in %.fs" % (expnames[ex],time.time()-st))
    
#     #  -------- -------- -------- -------- -------- -------- -------- --------
    
#     # Repeat Above for Early if Needed
#     if ex < 2:
#         print("Also Performing Early Calculations for %s" % (expnames[ex]))
        
#         # <This is the only part that changes, really>
#         varanom_in       = awi_early_anom[ex][vname]
#         ninoid_in        = ninoid_early[ex]
#         periodname       = early_str
        
#         # ---- Pretty much same as above (make into function?)
#         # Apply Low Pass Filter
#         st               = time.time()
#         varanom_in_lpf   = proc.pointwise_lp(varanom_in,lpfmon)
#         print("Applied LPF in %.fs" % (time.time()-st))
    
#         models           = fit_cmode_enso_hierarchy(varanom_in,ninoid_in)
#         models_lpf       = fit_cmode_enso_hierarchy(varanom_in_lpf,ninoid_in)
        
#         # Save Output 
#         st = time.time()
#         for ff in tqdm(range(3)):
#             fitname      = fitnames[ff]
            
#             # Save Regular
#             ds_out       = models[fitname]
#             outname      = "%s%s_fit_%s_%s_%s.nc" % (poster_path,vname,expnames[ex],periodname,fitname) # "cre_fit_NinoOnly
#             ds_out.to_netcdf(outname)
            
#             # Save LPF
#             ds_out_lpf   = models_lpf[fitname]
#             outname_lpf  = proc.addstrtoext(outname,"_LPF%02i" % lpfmon)
#             ds_out_lpf.to_netcdf(outname_lpf)
            
#             # End Fit Loop
        
#         # Also Save Variables (Original and LPF), and Nino3.4
#         # Full Anomaly
#         outname      = "%s%s_anomaly_%s_%s.nc" % (poster_path,vname,expnames[ex],periodname)
#         varanom_in.to_netcdf(outname)
#         # LPF Anomaly
#         outname_lpf  = proc.addstrtoext(outname,"_LPF%02i" % lpfmon)
#         varanom_in_lpf.to_netcdf(outname_lpf)
#         # Nino3.4 Index
#         outname_nino =  "%sNino34_%s_%s.nc" % (poster_path,expnames[ex],periodname)
#         ninoid_in.to_netcdf(outname_nino)
#         # --- Same as for late loop Above ---------------
        
#         print("Early Period Output for %s saved in %.fs" % (expnames[ex],time.time()-st))
#         # End Early Period Calculation
#     # End Experiment Loop
        
        
        
        
        
    
    


    



# # # Loop by variable
# # varanom_loop     = [varanom_in,varanom_in_lpf]
# # #
# # vv = 0 
# # #for vv in range(2)
# # st               = time.time()
# # v_in             = varanom_loop[vv]
# # # Fit for ENSO Only
# # ninofit_var      = proc.xr_linfit(v_in,ninoid_in,deg=1) 
# # # Fit for Cmode + ENSO Only
# # ctonefit_var     = ut.fit_ctone_pointwise(v_in,ninoid_in,linear_enso=True,fix_frequency_offset=True)
# # # Fit for Cmode Only
# # ctoneonlyfit_var = ut.fit_ctone_pointwise(v_in,ninoid_in,linear_enso=False,fix_frequency_offset=True)
# # print("Completed All Fits in %.2fs" % (time.time()-st))





#     # "%s%s_fit_%s" % (outpath,vname,periodname,)

# # #%% Optionally Save Nino3.4 Output...

# # for ex in tqdm(range(5)):
    
# #     ncname_out = "%sNino34_%s_%s.nc" % (poster_path,)


