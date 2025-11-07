#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Visualize output of lag_regression enso

based on structure of visualize_enso_regressions

Created on Fri Oct 10 09:46:49 2025

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

amvpath  = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut

#%% User Edits

expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"] #[glorys]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"] #["GLORYS Reanalysis",]

#ensoid_name     = "nino34"
standardize     = False
ninoid_name     = "nino34"

monmode = "All_Months"
datpath = "/home/niu4/gliu8/projects/scrap/TP_crop/lag_regressions/" + monmode + "/"
figpath = "/home/niu4/gliu8/figures/bydate/2025-10-14/"
proc.makedir(figpath)

# Other Parameters

vname       = "lcc"
leadlags    = np.arange(-12,13)


#%% Variable Dictionary (copied form visualize_mean_state)

sstdict = dict(
    vname      = "sst",
    longname   = "sea surface temperature",
    units      = "$\degree C$",
    cints_mean = np.arange(18,35.5,.5),
    cmap_mean  = "cmo.balance",
    cints_enso = np.arange(-2.5,2.6,0.1),
    
    )

sshdict = dict(
    vname      = "ssh",
    longname   = "sea surface height",
    units      = "$meters$",
    cints_mean = np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.deep_r",
    )

strdict = dict(
    vname      = "str",
    longname   = "Surface Thermal Radiation",
    units      = r"$\frac{J}{m^2}$",
    cints_mean = np.arange(-100,-22.5,2.5),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "plasma_r",
    cints_enso = np.arange(-5,5.5,0.5),
    )

ssrdict = dict(
    vname      = "ssr",
    longname   = "Surface Solar Radiation",
    units      = r"$\frac{J}{m^2}$",
    cints_mean = np.arange(120,290,10),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.solar",
    cints_enso = np.arange(-24,26,2),
    )

txsurdict = dict(
    vname      = "tx_sur",
    longname   = "Zonal Wind Stress to Ocean",
    units      = r"$\frac{m}{s^2}$",
    cints_mean = np.arange(-0.14,0.15,0.01),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.curl",
    cints_enso = np.arange(-0.012,0.013,0.001),
    )

tysurdict = dict(
    vname      = "ty_sur",
    longname   = "Meridional Wind Stress to Ocean",
    units      = r"$\frac{m}{s^2}$",
    cints_mean = np.arange(-0.08,0.085,0.005),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.curl",
    )


lccdict = dict(
    vname       = "lcc",
    longname    = "low cloud cover",
    units       ="fraction",
    cints_mean  = np.arange(0,1.1,0.1),
    cmap_mean   = 'cmo.ice',
    )

tccdict = dict(
    vname       = "tcc",
    longname    = "total cloud cover",
    units       ="fraction",
    cints_mean  = np.arange(0,1.1,0.1),
    cmap_mean   = 'cmo.ice',
    )

ttrdict = dict(
    vname      = "ttr",
    longname   = "TOA Net Thermal Radiation",
    units      = r"$\frac{J}{m^2}$",
    cints_mean = np.arange(-300,-195,5),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "plasma_r",
    )

# Just copy to make the clear sky version
ttrcdict = ttrdict.copy()
ttrcdict['vname'] = 'ttrc'
ttrcdict['longname'] = "TOA Clear-Sky Thermal Radiation",

tsrdict = dict(
    vname      = "tsr",
    longname   = "TOA Net Shortwave Radiation",
    units      = r"$\frac{J}{m^2}$",
    cints_mean = np.arange(220,370,10),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.solar",
    )

tsrcdict = dict(
    vname      = "tsrc",
    longname   = "TOA Clear-Sky Shortwave Radiation",
    units      = r"$\frac{J}{m^2}$",
    cints_mean = np.arange(350,384,2),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.solar",
    )

cpdict = dict(
    vname      = "cp",
    longname   = "Convective Precipitation",
    units      = r"$\frac{mm}{day}$",
    cints_mean = np.arange(0,11.5,0.5),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.rain",
    )

lspdict = dict(
    vname      = "lsp",
    longname   = "Stratiform Precipitation",
    units      = r"$\frac{mm}{day}$",
    cints_mean = np.arange(0,11.5,0.5),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.rain",
    )

sshfdict = dict(
    vname      = "sshf",
    longname   = "Surface sensible heat flux",
    units      = r"$\frac{J}{m^2}$",
    cints_mean = np.arange(-100,-22.5,2.5),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.thermal",
    )


slhfdict = dict(
    vname      = "slhf",
    longname   = "Surface latent heat flux",
    units      = r"$\frac{J}{m^2}$",
    cints_mean = np.arange(-200,-45,5),#np.arange(0,1.05,0.05),#np.arange(18,35.5,.5),
    cmap_mean  = "cmo.thermal",
    )


eisdict = dict(
    vname     = "eis",
    longname  = "Estimated inversion strength",
    units     = "$EIS$",
    cints_mean = np.arange(-5,5.5,0.5),
    cmap_mean = 'cmo.balance',
    cints_enso = np.arange(-1,1.1,0.1),
    
    )

dicts_byvar = [sstdict,
               sshdict,
               strdict,
               ssrdict,
               txsurdict,
               tysurdict,
               lccdict,
               tccdict,
               ttrdict,
               ttrcdict,
               tsrdict,
               cpdict,
               lspdict,
               slhfdict,
               sshfdict,
               eisdict,
               ]

vnames_all       = [d["vname"] for d in dicts_byvar]
vdicts           = dict(zip(vnames_all,dicts_byvar))


# Bounding Boxes from Jin et al. 2020 Eqn. 6.6  -----
bbox_cep        = [150      , -130+360 , -5, 5]   # Central Equatorial Pacific, for [tau_x], 
bbox_nino3      = [-150+360 , -90+360  , -5, 5]  # Nino 3 Box: For SST, <tau_x>
bbox_nino34     = [-170+360 , -120+360 , -5, 5]
bbox_epac       = [-155+360 , -80+360  , -5, 5]  # Eastern Pacific (for h_e calculation)
bbox_wpac       = [120      , -155+360 , -5, 5]  # Western Pacific (for h_w calculation)

bboxes       = [bbox_cep,bbox_nino3,bbox_nino34,bbox_epac,bbox_wpac]
bbnames_long = ["Central Equatorial Pacific","$Ni\tilde{n}o3$","$Ni\tilde{n}o3.4$","Tropical Eastern Pacific","Tropical Western Pacific"]
bbnames      = ["CEP","nino3","nino34","EPac","WPac"]
bbcols       = ["blue","red","goldenrod","hotpink","magenta"]
bbls         = ["solid","dashed","dotted","dashed","dotted"]

#%% Other Plotting Variables
projd   = ccrs.PlateCarree()
proj    =  ccrs.PlateCarree(central_longitude=180)

#%% Load ENSO ID

ensoids = [ut.load_ensoid(expname,ninoid_name,standardize=standardize) for expname in expnames]



#%% Load the Lag Regressions
# LagRegression_AllMonths_TCo2559-DART-1950C_ttrc_nino34_standardize0_lag-12to12.nc
nexps  = len(expnames)
ds_all = []
sigs   = []
for ex in tqdm.tqdm(range(nexps)):
    ncname = "%sLagRegression_AllMonths_%s_%s_%s_standardize%i_lag%02ito%02i.nc" % (datpath,expnames[ex],vname,ninoid_name,
                                                                                    standardize,leadlags[0],leadlags[-1])
    try:
        ds     = xr.open_dataset(ncname).load()
    except:
        print("Could not find %s, %s" % (vname,expnames[ex]))
        ds_all.append(None)
        sigs.append(None)
        continue
    
    sigs.append(ds.sig) # Note since varcheck not converts ds --> da, should grab sig first
    
    ds     = ut.varcheck(ds,vname,expnames[ex]) 
    
    ds_all.append(ds)
    
#%% Check Bounding Boxes
ex      = 4
il      = 12
vmax    = 1#0.012

bbox_test=False

plot_bbox = [2,3,4]#np.arange(len(bboxes))#bboxes None#[0,1]


# ---------
indict  = vdicts[vname]

cmap    = 'cmo.tarn'#'cmo.curl'
cints   = np.arange(-0.1,0.12,0.02)#np.arange(-0.012,0.013,0.001)#np.arange(-2.5,2.6,0.1)
# ---------

color = "k"
linestyle='solid',
linewidth = 1#TCo1279


fig,ax  = ut.init_tp_map()

if ex == 0:
    sigint = 1
elif ex == 1:
    sigint = 5
elif ex < 4:
    sigint = 20
else:
    if vname in ['sst','ssh']:
        sigint = 12
    else:
        sigint = 40
            

plotvar = ds_all[ex].isel(lag=il)

# Plot the variable
if not bbox_test:
    pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,
                            transform=projd,cmap=cmap,vmin=-vmax,vmax=vmax)
    
    # Plot Significance
    plotmask = sigs[ex].isel(lag=il)#ds_all[ex].sig.isel(lag=il)
    lon      = plotmask.lon
    lat      = plotmask.lat
    
    if len(lon) < 500: # Adjust stippling based on length of longitude
        sigint = 1
    elif len(lon) < 1000:
        sigint = 5
    elif len(lon) < 4000:
        sigint = 20
    else:
        sigint = 40
    viz.plot_mask(lon[::sigint],
                  lat[::sigint],
                  plotmask.T[::sigint,::sigint],
                  reverse=False,ax=ax,proj=projd,geoaxes=True,markersize=0.65,color='gray')
else:
    it = 12
    itlag = 12    


for bb in range(len(bboxes)):
    if bb not in plot_bbox:
        continue
    
    viz.plot_box(bboxes[bb],leglab=bbnames_long[bb],
                 color=bbcols[bb],linestyle=bbls[bb],proj=projd,ax=ax,linewidth=1.25)
    #ax.set_extent([-180,180,-90,90])
    ax.set_title(bbnames_long[bb])
    

if bbox_test:
    figname = "%sBBox_Test_ENSO_Lag_Regression_AllMonths_%s_%s_iter%02i_lag%02i.png" % (figpath,vname,expnames[ex],it,itlag)
    clab    = ""
    title   = ""
else:
    ax.legend()
    figname = "%sTest_ENSO_Lag_Regression_AllMonths_%s_%s_iter%02i_lag%02i.png" % (figpath,vname,expnames[ex],it,itlag)
    clab = "%s  [%s]" % (vname,indict['units'])
    title = "Lag %i %s-ENSO Regression, %s" % (plotvar.lag.data.item(),vname,expnames_long[ex])
cb = fig.colorbar(pcm,ax=ax,pad=0.01,fraction=0.01)
cb.set_label(clab)
ax.set_title(title)

plt.savefig(figname,dpi=150,bbox_inches='tight',transparent=True)

plt.show()

#%% Plot SST Patterns

ex          = 1
lag         = 0
cmap        = 'cmo.tarn'
cints       = np.arange(-0.1,0.15,0.05)#[0,] # #np.arange(-0.1,0.1,0.5)

#cints   = np.arange(0,1.05,.05)
#cints   = np.arange(-1.5,1.6,0.1) # SST
#cints   = np.arange(-0.05,0.06,0.01) #np.arange(-1,1.1,0.1)#np.arange(-5,5.5,0.5)#np.arange(-24,26,2)
vmax        = cints[0]
indict      = vdicts[vname]

for ex in range(nexps):
    it = 0
    

    
    for il in tqdm.tqdm(range(len(leadlags))):
        fig,ax  = ut.init_tp_map()
        
        
        #plotvar = ds_all[ex][vname].sel(lag=lag)
        try:
            plotvar = ds_all[ex].isel(lag=il)
            sigin     = sigs[ex]
        except:
            continue
        
        itlag    = plotvar.lag.data.item()
        
        # Plot the variable
        pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,
                                transform=projd,cmap=cmap,vmin=-vmax,vmax=vmax)
        
        cl      = ax.contour(plotvar.lon,plotvar.lat,plotvar,
                             linewidths=0.55,levels=cints,
                                transform=projd,colors='k')
        ax.clabel(cl)
        
        
        # Plot Significance
        plotmask = sigin.isel(lag=il)
        lon      = plotmask.lon
        lat      = plotmask.lat
        if len(lon) < 500: # Adjust stippling based on length of longitude
            sigint = 1
        elif len(lon) < 1000:
            sigint = 5
        elif len(lon) < 4000:
            sigint = 20
        else:
            sigint = 40
        # elif ex < 4:
        #     sigint = 20
        # else:
        #     if vname in ['sst','ssh']:
        #         sigint = 12
        #     else:
        #         sigint = 40
        
        viz.plot_mask(lon[::sigint],
                      lat[::sigint],
                      plotmask.T[::sigint,::sigint],
                      reverse=False,ax=ax,proj=projd,geoaxes=True,markersize=0.65,color='gray')
        
        cb      = fig.colorbar(pcm,ax=ax,fraction=0.01,pad=0.01)
        cb.set_label("%s [%s per $1\sigma$ %s]" % (vname,indict['units'],ninoid_name))
        ax.set_title("AWI-CM3 (%s) %s Regression, Lag %02i" % (expnames_long[ex],indict['longname'],itlag))
        
        figname = "%sENSO_Lag_Regression_AllMonths_%s_%s_iter%02i_lag%02i.png" % (figpath,vname,expnames[ex],it,itlag)
        plt.savefig(figname,dpi=150,bbox_inches='tight')
        it+=1

        plt.close()

#%%
for ex in range(nexps):
    plotmask = ds_all[ex].sig.isel(lag=il)
    lon      = plotmask.lon
    lat      = plotmask.lat
    print(expnames_long[ex])
    print(lon.shape)
    print(lat.shape)
    print("\n")
