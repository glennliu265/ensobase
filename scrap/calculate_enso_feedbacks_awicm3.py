#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Estimate ENSO feedbacks on Regridded AWI-CM3 Output

Created on Mon Jan 12 14:41:57 2026

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

#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% User Edits

# Simulation Names and Variables-----
expname         = "TCo319_ctl1950d"
vnames          = ["sst","temp75","ssh","uvel_ML50","vvel_ML50","wvel50","tx_sur","qnet"]
tstart          = '1950-01-01'
tend            = '2100-12-31'


# Load UNANOMALIZED Data
datpath         = "/home/niu4/gliu8/projects/scrap/regrid_1x1/"





#%% Plotting Variables
proj            = ccrs.PlateCarree()
figpath         = "/home/niu4/gliu8/figures/bydate/2026-01-13/"
proc.makedir(figpath)


#%% Load Variables (~14 sec for regridded output)

dsall = []
for vname in tqdm.tqdm(vnames):
    ncname = "%s%s_%s_regrid1x1.nc" % (datpath,expname,vname)
    ds     = xr.open_dataset(ncname)[vname].load()
    dsall.append(ds)
    

    
#%% Crop to time, deseason, detrend

lon   = ds[0].lon
lat   = ds[0].lat
times = ds.time_counter.data

def preprocess(ds,tstart,tend,lon,lat,times):
    ds  = ut.standardize_names(ds)
    ds  = ds.sel(time=slice(tstart,tend))
    dsa = proc.xrdeseason(ds)
    dsa = proc.xrdetrend(dsa)
    # Replace Lat Lon Time to match dimensions (they should match...)
    dsa['lon'] = lon
    dsa['lat'] = lat
    dsa['time'] = times # Some are in middle (12-16) others are at end (12-31)
    return dsa.squeeze()

st      = time.time()
dsanoms = [preprocess(ds,tstart,tend,lon,lat,times) for ds in dsall]
print("Preprocessed data in %.2fs" % (time.time()-st))


dsanoms = xr.merge(dsanoms)
#[print(ds.time) for ds in dsanoms]
#[print(ds.shape) for ds in dsanoms]
print("Preprocessed data in %.2fs" % (time.time()-st))

#%% Take Necessary Area-Averages

# Note, Assumes vnames are in the correct order
def aavg_region(ds,bbox):
    dsreg = proc.sel_region_xr(ds,bbox)
    return proc.area_avg_cosweight(dsreg)


# Bounding Boxes from Jin et al. 2020 Eqn. 6.6  ----- Taken from commons.py
bbox_cep        = [150      , -130+360 , -5, 5]  # Central Equatorial Pacific, for [tau_x], 
bbox_nino3      = [-150+360 , -90+360  , -5, 5]  # Nino 3 Box: For SST, <tau_x>
bbox_nino34     = [-170+360 , -120+360 , -5, 5]
bbox_epac       = [-155+360 , -80+360  , -5, 5]  # Eastern Pacific (for h_e calculation)
bbox_wpac       = [120      , -155+360 , -5, 5]  # Western Pacific (for h_w calculation)

nino3       = aavg_region(dsanoms.sst,bbox_nino3)
nino34      = aavg_region(dsanoms.sst,bbox_nino34)

taux_cep    = aavg_region(dsanoms.tx_sur,bbox_cep)
taux_nino3  = aavg_region(dsanoms.tx_sur,bbox_nino3)

he          = aavg_region(dsanoms.ssh,bbox_epac)
hw          = aavg_region(dsanoms.ssh,bbox_wpac)
hdiff       = he-hw

w_nino3     = aavg_region(dsanoms.wvel50,bbox_nino3)
u_nino3     = aavg_region(dsanoms.uvel_ML50,bbox_nino3)
v_nino3     = aavg_region(dsanoms.uvel_ML50,bbox_nino3)

Q_nino3     = aavg_region(dsanoms.qnet,bbox_nino3)
Tsub_nino3  = aavg_region(dsanoms.temp75,bbox_nino3)

Q_nino3_conv = ut.varcheck(Q_nino3,'slhf',expname) # Using 'slhf' to trigger flux conversion... manual update later to accept 'qnet'


rho   = 1025 # [kg m^3]
Cp    = 3994 # [J / (kg K)]
H     = 50   # [meters]
dtmon = 3600*24*30


Q_nino3_degC = Q_nino3_conv/(rho*Cp*H) * dtmon

#%% Estimate Coefficients Using (Multiple) Linear Regression

mu_a        = ut.mlr_ccfs([nino3,],taux_cep,standardize=False,verbose=True)   # Eqn 6.6
mu_a_star   = ut.mlr_ccfs([nino3,],taux_nino3,standardize=False,verbose=True) # Eqn 6.7
beta_h      = ut.mlr_ccfs([taux_cep,],hdiff,standardize=False,verbose=True)   # Eqn 6.8

w_coeffs    = ut.mlr_ccfs([-taux_cep,-taux_nino3,hw],w_nino3,standardize=False,verbose=True) # Eqn 6.9
u_coeffs    = ut.mlr_ccfs([taux_cep,taux_nino3,hw],w_nino3,standardize=False,verbose=True) # Eqn 6.10
v_coeffs    = ut.mlr_ccfs([taux_cep,taux_nino3,hw],w_nino3,standardize=False,verbose=True) # Eqn 6.11

a_h         = ut.mlr_ccfs([he,],Tsub_nino3,standardize=False,verbose=True)                      # Eqn 6.12
alpha       = ut.mlr_ccfs([-1*nino3,],Q_nino3_degC,standardize=False,verbose=True)   # Eqn 6.13







#%% Perform the parameterizations (using MLR)

def init_fitplot():
    fig             = plt.figure(figsize=(12,3))
    gs              = gridspec.GridSpec(4,10)

    # --------------------------------- Timeseries
    ax11            = fig.add_subplot(gs[:,:8],)

    # --------------------------------- # Scatterplot
    ax22            = fig.add_subplot(gs[:,8:],)
    ax22            = viz.add_axlines(ax22)
    #ax22.set_aspect('equal')
    
    return fig,gs,ax11,ax22

# =====================
#%% Plot for Mu_a
# =====================

target      = taux_cep
target_name = r"[$\tau_x$]" # r"$ \langle \tau_x \rangle>" # r"[$\tau_x$]"
target_unit = r"N m$^{-2}$"
model_params = mu_a
model_name  = r"$\mu_a T_E$" 
coeff_name  = r"$\mu_a$"
coeff_fn    = "mu_a"

predictor      = nino3
predictor_name = r"$T_E$"
predictor_unit = r"$\degree C$"



fig,gs,ax11,ax22 = init_fitplot()

# Axis 1: Timeseries ----------------------------------------------------------
ax               = ax11

# Plot Target
plotvar          = target
ax.plot(plotvar.time,plotvar,c ='gray',lw=1.5,label=target_name)

# Plot Prediction
plotvar                 = model_params['pred']
ax.plot(times,plotvar,c ='red',label=model_name,ls='dashed',lw=0.5)

# Labeling
coeff_val        = model_params['coeffs'][0]
ax.set_title(r"%s = $%.2e$" % (coeff_name,coeff_val))
ax.set_ylabel(r"%s (%s)" % (target_name,target_unit))
#ax.set_xlim([times[0],times[12*45]])
ax.legend()


# Axis 2: Scatterplot  --------------------------------------------------------
ax = ax22

# Plot Scatter
# r2 = model_params['r2']
# ax.scatter(predictor,target,c=times,label=r"$R^2$=%.2e" % (r2),alpha=0.75,s=2.5)
r = np.corrcoef(predictor,target)[0,1]
ax.scatter(predictor,target,c=times,label=r"$R$=%.2f" % (r),alpha=0.75,s=2.5)

# Laeling
ax.set_xlabel(r"%s [%s]" % (predictor_name,predictor_unit))
ax.legend(frameon=False)
#ax.set_aspect(1, adjustable='box')
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()

# Save Figure and Show
savename = "%sBWJ_Feedback_Estimate_%s_%s_regrid1x1.png" % (figpath,coeff_fn,expname)
plt.savefig(savename,dpi=150,transparent=True,bbox_inches='tight')
plt.show()

# ==================
#%% mu_a_star
# ==================

target      = taux_nino3
target_name = r"$ \langle \tau_x \rangle $"
target_unit = r"N m$^{-2}$"
model_params = mu_a_star
model_name  = r"$\mu_a^* T_E$" 
coeff_name  = r"$\mu_a^*$"
coeff_fn    = "mu_a_star"

predictor      = nino3
predictor_name = r"$T_E$"
predictor_unit = r"$\degree C$"


fig,gs,ax11,ax22 = init_fitplot()

# Axis 1: Timeseries ----------------------------------------------------------
ax               = ax11

# Plot Target
plotvar          = target
ax.plot(plotvar.time,plotvar,c ='gray',lw=1.5,label=target_name)

# Plot Prediction
plotvar                 = model_params['pred']
ax.plot(times,plotvar,c ='red',label=model_name,ls='dashed',lw=0.5)

# Labeling
coeff_val        = model_params['coeffs'][0]
ax.set_title(r"%s = $%.2e$" % (coeff_name,coeff_val))
ax.set_ylabel(r"%s (%s)" % (target_name,target_unit))
#ax.set_xlim([times[0],times[12*45]])
ax.legend()


# Axis 2: Scatterplot  --------------------------------------------------------
ax = ax22

# Plot Scatter
# r2 = model_params['r2']
# ax.scatter(predictor,target,c=times,label=r"$R^2$=%.2e" % (r2),alpha=0.75,s=2.5)
r = np.corrcoef(predictor,target)[0,1]
ax.scatter(predictor,target,c=times,label=r"$R$=%.2f" % (r),alpha=0.75,s=2.5)

# Laeling
ax.set_xlabel(r"%s [%s]" % (predictor_name,predictor_unit))
ax.legend(frameon=False)
#ax.set_aspect(1, adjustable='box')
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()

# Save Figure and Show
savename = "%sBWJ_Feedback_Estimate_%s_%s_regrid1x1.png" % (figpath,coeff_fn,expname)
plt.savefig(savename,dpi=150,transparent=True,bbox_inches='tight')
plt.show()


# ======================
#%% beta_h
# ======================

target      = hdiff
target_name = r"$ ssh_e - ssh_w $"
target_unit = r"$m$"
model_params = beta_h
model_name  = r"$\beta_h [\tau_x]$" 
coeff_name  = r"$\beta_h$"
coeff_fn    = "beta_h"

predictor      = taux_cep
predictor_name = r"[$\tau_x$]"
predictor_unit = r"N m$^{-2}$"


fig,gs,ax11,ax22 = init_fitplot()

# Axis 1: Timeseries ----------------------------------------------------------
ax               = ax11

# Plot Target
plotvar          = target
ax.plot(plotvar.time,plotvar,c ='gray',lw=1.5,label=target_name)

# Plot Prediction
plotvar                 = model_params['pred']
ax.plot(times,plotvar,c ='red',label=model_name,ls='dashed',lw=0.5)

# Labeling
coeff_val        = model_params['coeffs'][0]
ax.set_title(r"%s = $%.2e$" % (coeff_name,coeff_val))
ax.set_ylabel(r"%s (%s)" % (target_name,target_unit))
#ax.set_xlim([times[0],times[12*45]])
ax.legend()


# Axis 2: Scatterplot  --------------------------------------------------------
ax = ax22

# Plot Scatter
# r2 = model_params['r2']
# ax.scatter(predictor,target,c=times,label=r"$R^2$=%.2e" % (r2),alpha=0.75,s=2.5)
r = np.corrcoef(predictor,target)[0,1]
ax.scatter(predictor,target,c=times,label=r"$R$=%.2f" % (r),alpha=0.75,s=2.5)

# Laeling
ax.set_xlabel(r"%s [%s]" % (predictor_name,predictor_unit))
ax.legend(frameon=False)
#ax.set_aspect(1, adjustable='box')
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()

# Save Figure and Show
savename = "%sBWJ_Feedback_Estimate_%s_%s_regrid1x1.png" % (figpath,coeff_fn,expname)
plt.savefig(savename,dpi=150,transparent=True,bbox_inches='tight')
plt.show()


# ======================
#%% a_h
# ======================



target      = Tsub_nino3
target_name = r"$\rangle T_{sub} \rangle$"
target_unit = r"$\degree C$"
model_params = a_h
model_name  = r"$a_h h_e$" 
coeff_name  = r"$a_h$"
coeff_fn    = "a_h"

predictor      = he

predictor_name = r"[$ssh_e$]"
predictor_unit = r"m"


fig,gs,ax11,ax22 = init_fitplot()

# Axis 1: Timeseries ----------------------------------------------------------
ax               = ax11

# Plot Target
plotvar          = target
ax.plot(plotvar.time,plotvar,c ='gray',lw=1.5,label=target_name)

# Plot Prediction
plotvar                 = model_params['pred']
ax.plot(times,plotvar,c ='red',label=model_name,ls='dashed',lw=0.5)

# Labeling
coeff_val        = model_params['coeffs'][0]
ax.set_title(r"%s = $%.2e$" % (coeff_name,coeff_val))
ax.set_ylabel(r"%s (%s)" % (target_name,target_unit))
#ax.set_xlim([times[0],times[12*45]])
ax.legend()


# Axis 2: Scatterplot  --------------------------------------------------------
ax = ax22

# Plot Scatter
# r2 = model_params['r2']
# ax.scatter(predictor,target,c=times,label=r"$R^2$=%.2e" % (r2),alpha=0.75,s=2.5)
r = np.corrcoef(predictor,target)[0,1]
ax.scatter(predictor,target,c=times,label=r"$R$=%.2f" % (r),alpha=0.75,s=2.5)

# Laeling
ax.set_xlabel(r"%s [%s]" % (predictor_name,predictor_unit))
ax.legend(frameon=False)
#ax.set_aspect(1, adjustable='box')
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()

# Save Figure and Show
savename = "%sBWJ_Feedback_Estimate_%s_%s_regrid1x1.png" % (figpath,coeff_fn,expname)
plt.savefig(savename,dpi=150,transparent=True,bbox_inches='tight')
plt.show()


# ======================
#%% alpha
# ======================

target       = Q_nino3_degC
target_name  = r"$\langle Q \rangle$"
target_unit  = r"$\degree C$ $month^{-1}$"
model_params = alpha
model_name   = r"$-\alpha T_E$" 
coeff_name   = r"$\alpha$"
coeff_fn     = "alpha"


predictor      = -1*nino3
predictor_name = r"$T_E$"
predictor_unit = r"$\degree C$"


fig,gs,ax11,ax22 = init_fitplot()

# Axis 1: Timeseries ----------------------------------------------------------
ax             = ax11

# Plot Target
plotvar          = target
ax.plot(plotvar.time,plotvar,c ='gray',lw=1.5,label=target_name)

# Plot Prediction
plotvar                 = model_params['pred']
ax.plot(times,plotvar,c ='red',label=model_name,ls='dashed',lw=0.5)

# Labeling
coeff_val        = model_params['coeffs'][0]
ax.set_title(r"%s = $%.2e$" % (coeff_name,coeff_val))
ax.set_ylabel(r"%s (%s)" % (target_name,target_unit))
#ax.set_xlim([times[0],times[12*45]])
ax.legend()


# Axis 2: Scatterplot  --------------------------------------------------------
ax = ax22

# Plot Scatter
# r2 = model_params['r2']
# ax.scatter(predictor,target,c=times,label=r"$R^2$=%.2e" % (r2),alpha=0.75,s=2.5)
r = np.corrcoef(predictor,target)[0,1]
ax.scatter(predictor,target,c=times,label=r"$R$=%.2f" % (r),alpha=0.75,s=2.5)

# Laeling
ax.set_xlabel(r"%s [%s]" % (predictor_name,predictor_unit))
ax.legend(frameon=False)
#ax.set_aspect(1, adjustable='box')
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()

# Save Figure and Show
savename = "%sBWJ_Feedback_Estimate_%s_%s_regrid1x1.png" % (figpath,coeff_fn,expname)
plt.savefig(savename,dpi=150,transparent=True,bbox_inches='tight')
plt.show()

# ======================
#%% u
# ======================


target       = u_nino3
target_name  = r"$\langle u \rangle$"
target_unit  = r"$m s^{-1}$"
model_params = u_coeffs

prednames    = [r"\langle \tau_x \rangle", r"[\tau_x]","h_w"]
coeff_name   = [r"\beta_{ul}",r"\beta_{ur}",r"\beta_{uh}"]
model_name   = r"Fit $\langle u \rangle$"
coeff_fn     = "u_coeffs"

predictor      = model_params['pred']
predictor_name = []
for ii in range(3):
    predictor_name.append(coeff_name[ii] + prednames[ii])
predictor_name = " + ".join(predictor_name)
predictor_unit = r"$m s^{-1}$"


fig,gs,ax11,ax22 = init_fitplot()

# Axis 1: Timeseries ----------------------------------------------------------

ax               = ax11

# Plot Target
plotvar = target
ax.plot(plotvar.time,plotvar,c ='gray',lw=1.5,label=target_name)

# Plot Prediction
plotvar                 = model_params['pred']
ax.plot(times,plotvar,c = 'red',label=model_name,ls='dashed',lw=0.5)

# Labeling
coeff_val               = model_params['coeffs']
title = []
for ii in range(3):
    title.append(r"$%s=%.2e$" % (coeff_name[ii],coeff_val[ii]))
title = ", ".join(title)
ax.set_title(title)
ax.set_ylabel(r"%s (%s)" % (target_name,target_unit))
#ax.set_xlim([times[0],times[12*45]])
ax.legend()


# Axis 2: Scatterplot  --------------------------------------------------------
ax = ax22

# Plot Scatter
# r2 = model_params['r2']
# ax.scatter(predictor,target,c=times,label=r"$R^2$=%.2e" % (r2),alpha=0.75,s=2.5)
r = np.corrcoef(predictor,target)[0,1]
ax.scatter(predictor,target,c=times,label=r"$R$=%.2f" % (r),alpha=0.75,s=2.5)

# Laeling
ax.set_xlabel(r"$%s$ [%s]" % (predictor_name,predictor_unit))
ax.legend(frameon=False)
#ax.set_aspect(1, adjustable='box')
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()

# Save Figure and Show
savename = "%sBWJ_Feedback_Estimate_%s_%s_regrid1x1.png" % (figpath,coeff_fn,expname)
plt.savefig(savename,dpi=150,transparent=True,bbox_inches='tight')
plt.show()

# ======================
#%% v
# ======================


target       = v_nino3
target_name  = r"$\langle v \rangle_A$"
target_unit  = r"$m s^{-1}$"
model_params = v_coeffs

#prednames    = [r"\langle \tau_x \rangle", r"[\tau_x]","h_w"]
coeff_name   = [r"\beta_{vl}",r"\beta_{vr}",r"\beta_{vh}"]
model_name   = r"Fit $\langle v \rangle_A$"
coeff_fn     = "v_coeffs"

predictor      = model_params['pred']
predictor_name = []
for ii in range(3):
    predictor_name.append(coeff_name[ii] + prednames[ii])
predictor_name = " + ".join(predictor_name)
predictor_unit = r"$m s^{-1}$"


fig,gs,ax11,ax22 = init_fitplot()

# Axis 1: Timeseries ----------------------------------------------------------

ax               = ax11

# Plot Target
plotvar = target
ax.plot(plotvar.time,plotvar,c ='gray',lw=1.5,label=target_name)

# Plot Prediction
plotvar                 = model_params['pred']
ax.plot(times,plotvar,c = 'red',label=model_name,ls='dashed',lw=0.5)

# Labeling
coeff_val               = model_params['coeffs']
title = []
for ii in range(3):
    title.append(r"$%s=%.2e$" % (coeff_name[ii],coeff_val[ii]))
title = ", ".join(title)
ax.set_title(title)
ax.set_ylabel(r"%s (%s)" % (target_name,target_unit))
#ax.set_xlim([times[0],times[12*45]])
ax.legend()

# Axis 2: Scatterplot  --------------------------------------------------------
ax = ax22

# Plot Scatter
# r2 = model_params['r2']
# ax.scatter(predictor,target,c=times,label=r"$R^2$=%.2e" % (r2),alpha=0.75,s=2.5)
r = np.corrcoef(predictor,target)[0,1]
ax.scatter(predictor,target,c=times,label=r"$R$=%.2f" % (r),alpha=0.75,s=2.5)

# Laeling
ax.set_xlabel(r"$%s$ [%s]" % (predictor_name,predictor_unit))
ax.legend(frameon=False)
#ax.set_aspect(1, adjustable='box')
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()

# Save Figure and Show
savename = "%sBWJ_Feedback_Estimate_%s_%s_regrid1x1.png" % (figpath,coeff_fn,expname)
plt.savefig(savename,dpi=150,transparent=True,bbox_inches='tight')
plt.show()


# ======================
#%% w
# ======================


target       = w_nino3
target_name  = r"$\langle w \rangle$"
target_unit  = r"$m s^{-1}$"
model_params = w_coeffs

#prednames    = [r"\langle \tau_x \rangle", r"[\tau_x]","h_w"]
coeff_name   = [r" - \beta_{wl}",r" - \beta_{wr}",r" + \beta_{wh}"]
model_name   = r"Fit $\langle w \rangle$"
coeff_fn     = "w_coeffs"

predictor      = model_params['pred']
predictor_name = []
for ii in range(3):
    predictor_name.append(coeff_name[ii] + prednames[ii])
predictor_name = "".join(predictor_name)
predictor_unit = r"$m s^{-1}$"


fig,gs,ax11,ax22 = init_fitplot()

# Axis 1: Timeseries ----------------------------------------------------------

ax               = ax11

# Plot Target
plotvar = target
ax.plot(plotvar.time,plotvar,c ='gray',lw=1.5,label=target_name)

# Plot Prediction
plotvar                 = model_params['pred']
ax.plot(times,plotvar,c = 'red',label=model_name,ls='dashed',lw=0.5)

# Labeling
coeff_val               = model_params['coeffs']
title = []
for ii in range(3):
    title.append(r"$%s=%.2e$" % (coeff_name[ii],coeff_val[ii]))
title = ", ".join(title)
ax.set_title(title)
ax.set_ylabel(r"%s (%s)" % (target_name,target_unit))
#ax.set_xlim([times[0],times[12*45]])
ax.legend()

# Axis 2: Scatterplot  --------------------------------------------------------
ax = ax22

# Plot Scatter
# r2 = model_params['r2']
# ax.scatter(predictor,target,c=times,label=r"$R^2$=%.2e" % (r2),alpha=0.75,s=2.5)

r = np.corrcoef(predictor,target)[0,1]
ax.scatter(predictor,target,c=times,label=r"$R$=%.2f" % (r),alpha=0.75,s=2.5)

# Laeling
ax.set_xlabel(r"$%s$ [%s]" % (predictor_name,predictor_unit))
ax.legend(frameon=False)
#ax.set_aspect(1, adjustable='box')
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()

# Save Figure and Show
savename = "%sBWJ_Feedback_Estimate_%s_%s_regrid1x1.png" % (figpath,coeff_fn,expname)
plt.savefig(savename,dpi=150,transparent=True,bbox_inches='tight')
plt.show()

#%%

# Scrap Below



#%%Old 

fig,gs,ax11,ax22 = init_fitplot()

# Axis 1: Timeseries ----------------------------------------------------------
ax               = ax11

# Plot Target
plotvar          = taux_cep
ax.plot(plotvar.time,plotvar,c ='gray',lw=1.5,label=r"[$\tau_x$]")

# Plot Prediction
plotvar                 = mu_a['pred']
ax.plot(times,plotvar,c ='red',label=r"[$\mu_a T_E$]",ls='dashed',lw=0.5)

# Labeling
coeff_val        = mu_a['coeffs'][0]
ax.set_title(r"$\mu_a$ = $%.2f$" % (coeff_val))
ax.set_ylabel(r"[$\tau_x$] (N m$^{-2}$)")
ax.set_xlim([times[0],times[12*45]])
ax.legend()

# Axis 2: Scatterplot  --------------------------------------------------------
ax = ax22

# Plot Scatter
r2 = mu_a['r2']
ax.scatter(nino3,taux_cep,c=times,label=r"$R^2$=%.2e" % (r2),alpha=0.75,s=2.5)

# Laeling
ax.set_xlabel(r"$T_E$ [$\degree C$]")
ax.legend(frameon=False)
#ax.set_aspect(1, adjustable='box')
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()

# Save Figure and Show
savename = "%sBWJ_Feedback_Estimate_mu_a_%s_regrid1x1.png" % (figpath,expname)
plt.savefig(savename,dpi=150,transparent=True)
plt.show()

#%%




















