#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Low Cloud Unobscured Fluxes from CERES_FBCT
based on procedure outlined in Scott et al. 2020

   - Uses output preprocessed by [preprocess_CERES_FBCT.py]

Created on Wed Nov 26 13:37:51 2025

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
import importlib

from sklearn.linear_model import LinearRegression
import sklearn

#%% Import Custom Modules

amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% Load Output

figpath  = "/home/niu4/gliu8/figures/bydate/2025-12-02/"
proc.makedir(figpath)

datpath  = "/home/niu4/gliu8/share/CERES/processed/"
vnames   = ["tcc_cldtyp","ttr_cldtyp","tsr_cldtyp",
           "ttrc","tsrc","ttr","tsr"
           ]

ncsearch = datpath + "CERES_FBCT_%s_2002-07_to_2023-02.nc"

dsall = []
for vname in tqdm.tqdm(vnames):
    ncname = ncsearch % vname
    ds     = xr.open_dataset(ncname)[vname].load()
    dsall.append(ds)
    
    
dsall  = xr.merge(dsall)

idlow  = [0,1]
idhi   = [2,3,4,5,6]

#%% Calculation Options

recalc_frac = False # Recalculate Low/High Cloud Fractions

# =======================================================
#%% Part 1: Calculate Low/High/Unobscured Cloud Fractions
# =======================================================

if recalc_frac:
    print("Recalculating Cloud Fractions")
    L      = dsall.tcc_cldtyp.isel(press=idlow).sum(('press','opt'))/100
    U      = dsall.tcc_cldtyp.isel(press=idhi).sum(('press','opt'))/100
    Ln     = L/(1-U)
    
    #%% Visualize Cloud Fractions
    cloudfracs = [U,L,Ln]
    plotfracs  = cloudfracs
    plotmeans  = [ds.mean('time') for ds in plotfracs]
    plotnames  = ["Upper Cloud","Low Cloud","Unobscured Low Cloud"]
    
    proj      = ccrs.PlateCarree()
    
    fig,axs   = ut.init_globalmap(1,3,figsize=(24,4.5))
    
    for ii in range(3):
        ax = axs[ii]
        plotvar = plotmeans[ii]
        pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,
                                transform=proj,vmin=.05,vmax=.75,cmap='Blues_r')
        ax.set_title(plotnames[ii])
    cb         = viz.hcbar(pcm,ax=axs.flatten(),fraction=0.035)
    cb.set_label("Cloud Area Fraction")
    
    figname = "%sMeanCloudFraction_CERES_FBCT.png" % (figpath)
    plt.savefig(figname)
    plt.show()
    
    #%% Save variables
    
    vnames_frac = ["ucc","lcc","lncc"]
    for vv in range(3):
        fracname = vnames_frac[vv]
        dsout    = plotfracs[vv]
        ncout    = ncsearch % fracname
        print("Saving file to: %s" % ncout)
        dsout    = dsout.rename(fracname)
        dsout.to_netcdf(ncout)
else:
    print("Loading Cloud Fractions")
    cloudfracs  = []
    vnames_frac = ["ucc","lcc","lncc"]
    for vv in range(3):
        fracname = vnames_frac[vv]
        ncout    = ncsearch % fracname
        ds       = xr.open_dataset(ncout)[fracname].load()
        cloudfracs.append(ds)
    

#%%

# ========================================
#%% Part 2: Calculate the All Sky and Clear Sky (net) fluxes
# ========================================


#%% debugging

tsrmean = dsall.tsr_cldtyp.sum(('press','opt')).mean('time')
ttrmean = dsall.ttr_cldtyp.sum(('press','opt')).mean('time')

tsrnet  = (dsall.tsr_cldtyp * (dsall.tcc_cldtyp/100)).sum(('press','opt'))

tsrmean = dsall.tsr_cldtyp.isel(press=0,opt=0,time=0).plot()

tsrnet.mean('time').plot(vmin=0,vmax=160,cmap='cmo.solar'),plt.show()
(dsall.tsr-dsall.tsrc).mean('time').plot(vmin=0,vmax=160,cmap='cmo.solar'),plt.show()

#%% Try some other approaches

allsky_cldtyp = dsall.tsr_cldtyp - dsall.ttr_cldtyp

allsky_total  = (allsky_cldtyp * (dsall.tcc_cldtyp/100)).sum(('press','opt'))


dsall.tsr.mean('time').plot(vmin=50,vmax=350,cmap='cmo.solar'),plt.show()

#allsky_net = dsall.tsr =
#allsky_total.mean('time').plot(vmin=-100,vmax=100,cmap='cmo.balance'),plt.show()

#%% Lets try to see how the data works
lonf = 330 #240
latf = 50#20

dspt      = proc.selpt_ds(dsall,lonf,latf)

tcc_total = (dspt.tcc_cldtyp/100).sum(('opt','press'))
    
# Try with TSR ()
sumtsr    = (dspt.tsr_cldtyp * dspt.tcc_cldtyp/100).sum(('opt','press')) + dspt.tsrc * (1-tcc_total)
fig,ax    = plt.subplots(1,1,figsize=(12.5,4),constrained_layout=True)
ax.plot(dspt.tsr.time,dspt.tsr,lw=2.5,color='k',label='tsr')
ax.plot(sumtsr.time,sumtsr,lw=2,color='r',ls='dashed',label='tsr*frac')
ax.legend()
plt.show()

# Try with TTR
sumttr = (dspt.ttr_cldtyp * dspt.tcc_cldtyp/100).sum(('opt','press')) + dspt.ttrc * (1-tcc_total)
#sumttr = (dspt.ttr_cldtyp).sum(('opt','press'))
fig,ax = plt.subplots(1,1,figsize=(12.5,4),constrained_layout=True)
ax.plot(dspt.ttr.time,dspt.ttr,lw=2.5,color='k',label='ttr')
ax.plot(sumttr.time,sumttr,lw=2,color='r',ls='dashed',label='ttr*frac')
ax.legend()
plt.show()

#%% Try to calculate allsky

tcc_total        = (dsall.tcc_cldtyp/100).sum(('opt','press'))
allsky_total     = dsall.tsr - dsall.ttr


cre_cldtyp       = (dsall.tsr_cldtyp - dsall.ttr_cldtyp) #* (dspt.tcc_cldtyp/100)
clearsky         = (dsall.tsrc - dsall.ttrc) * (1-tcc_total)

# Calculate CRE by summing and multiplying by cloud fraction
cre_sum          = (cre_cldtyp * (dspt.tcc_cldtyp/100) ).sum(('opt','press'))

cre_sum.mean('time').plot(vmin=-100,vmax=100,cmap='cmo.balance'),plt.show()
cre_sum.mean('time').plot(vmin=-150,vmax=150,cmap='cmo.balance'),plt.show()




clearsky.mean('time').plot(vmin=-100,vmax=100,cmap='cmo.balance'),plt.show()
dsall.tsrc.mean('time').plot(vmin=50,vmax=350,cmap='cmo.solar'),plt.show()
(-1*dsall.ttrc.mean('time')).plot(vmin=-300,vmax=-140,cmap='cmo.dense_r'),plt.show()

dsall.tsr.mean('time').plot(vmin=50,vmax=350,cmap='cmo.solar'),plt.show()
(-1*dsall.ttr).mean('time').plot(vmin=-300,vmax=-140,cmap='cmo.dense_r'),plt.show()
allsky_total.mean('time').plot(vmin=-100,vmax=100,cmap='cmo.balance'),plt.show()


#%% Compute Fluxes by summing

st      = time.time()
tsr_sum = (dsall.tsr_cldtyp * (dspt.tcc_cldtyp/100) ).sum(('opt','press'))
ttr_sum = (dsall.ttr_cldtyp * (dspt.tcc_cldtyp/100) ).sum(('opt','press'))
print("Summed in %.2fs" % (time.time()-st))

tsrall  = tsr_sum + dsall.tsrc * (1-tcc_total)
ttrall  = ttr_sum + dsall.ttrc * (1-tcc_total)

tsrall.mean('time').plot(vmin=50,vmax=350,cmap='cmo.solar'),plt.show()
(-1*ttrall.mean('time')).plot(vmin=-300,vmax=-140,cmap='cmo.dense_r'),plt.show()



#%%

# Compute Net TOA flux variables
allsky_cldtyp = dsall.tsr_cldtyp - dsall.ttr_cldtyp
clearsky      = dsall.tsrc - dsall.ttrc
cre_cldtyp    = allsky_cldtyp - clearsky

# Compute cloud radiative effects
cre_sumlater  = cre_cldtyp.sum(('press','opt'))
cre_sumfirst  = allsky_cldtyp.sum(('press','opt')) - clearsky


# =============================================================================
# ========================================
#%% Part 3: Partition the Radiative Fluxes
# ========================================

# Compute Monthly Mean and Anomalies
#cloudfracs    = [U,L,Ln] # same as plotfracs
climmean_frac = [ds.groupby('time.month').mean('time') for ds in cloudfracs]
anom_frac     = [cloudfracs[vv].groupby('time.month') - climmean_frac[vv] for vv in range(3)]

Ubar,Lbar,Lnbar         = climmean_frac
Uprime,Lprime,Lnprime   = anom_frac

# Similarly, determine the anomalous and climatological values for other parameters
f       = dsall.tcc_cldtyp # tcc
fbar    = f.groupby('time.month').mean('time')
fprime  = f.groupby('time.month') - fbar

# Second order correction considering redistribution across optical and depth bins
fprime2_Lprime  = fprime - fbar * (Lprime/Lbar)
fprime2_Uprime  = fprime - fbar * (Uprime/Ubar)
fprime2_Lnprime = fprime - fbar * (Lnprime/Lnbar)

#%%

RL_term1 = (((allsky_bar*fbar)/Lbar).sum(('press','opt')) - clearsky_bar) * Lprime
RL_term2 = ((allsky_bar-clearsky_bar)*fprime2)


#%% Try to calculate RLn

# Fraction Variables
climmean_frac           = [ds.groupby('time.month').mean('time') for ds in cloudfracs]
anom_frac               = [cloudfracs[vv].groupby('time.month') - climmean_frac[vv] for vv in range(3)]
Ubar,Lbar,Lnbar         = climmean_frac
Uprime,Lprime,Lnprime   = anom_frac

# Flux Variables by level
Rpt     = dsall.tsr_cldtyp - dsall.ttr_cldtyp
Rptbar  = Rpt.groupby('time.month').mean('time')
Rclr    = dsall.tsrc - dsall.ttrc
Rclrbar = Rclr.groupby('time.month').mean('time')

# Fraction Variables by level
# Similarly, determine the anomalous and climatological values for other parameters
f       = dsall.tcc_cldtyp/100 # tcc
fbar    = f.groupby('time.month').mean('time')
fprime  = f.groupby('time.month') - fbar

# Correction Variables
# Second order correction considering redistribution across optical and depth bins
Lnratio         = Lnprime.groupby('time.month') / Lnbar
fprime2_Lnprime = fprime.isel(press=idlow) - fbar.isel(press=idlow) * (Lnratio).groupby('time.month') 

# Calculate Fractional Term
frac_mult = ( Lprime - Uprime.groupby('time.month')*Lnbar).groupby('time.month') / (1-Ubar)

# Calculate First Term (appears to be something wrong here)
#term1 =  ((Rptbar.isel(press=idlow) * fbar.isel(press=idlow)) / Lbar - Rclrbar).sum(('press','opt'))
term1 =  ((Rptbar.isel(press=idlow) * fbar.isel(press=idlow)) / Lbar).sum(('press','opt')) - Rclrbar

# Calcualte Second Term
Rdiff  = (Rptbar.isel(press=idlow) - Rclrbar)
term2  = (Rdiff * fprime2_Lnprime.groupby('time.month')).sum(('press','opt'))

# Calculate Full Term
RLnprime = term1*frac_mult.groupby('time.month') + term2#*fprime2_Lnprime


RLnprime.mean('time').plot(vmin=-100,vmax=100),plt.show()


RLnprime.isel(time=2).plot(vmin=-100,vmax=100,cmap='cmo.balance'),plt.show()

#%% Save the output

RLnprime = RLnprime.rename("creln")
outname  = ncsearch % "creln_anom"
RLnprime.to_netcdf(outname)


#%% Do a simple calculation




#%% Compute the 



#Lnbar.plot(vmin=5,vmax=75,cmap='Blues_r'),plt.show()






#press = dsall.press

