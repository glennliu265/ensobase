#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 11:23:36 2025

@author: gliu

"""

# Niu Paths

scmpath = "/home/niu4/gliu8/scripts/commons/stochmod/model"
amvpath = "/home/niu4/gliu8/scripts/commons"

sys.path.append(amvpath)
from amv import proc,viz


sys.path.append(scmpath)
import scm

ensopath = "/home/niu4/gliu8/scripts/ensobase"
sys.path.append(ensopath)
import utils as ut


#%% Shared variable names and experiment
# Copied from visualize composites

#datpath         = "/Users/gliu/Downloads/02_Research/01_Projects/07_ENSO/01_Data/TP_Crop/composites/"


# Simulation Names -----
expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090","TCo2559-DART-1950C"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090","5km 1950"]

timecrops       = [[1950,2100],None,None,None,None]


expcols = ["cornflowerblue",'lightcoral',
           "slateblue","firebrick",
           "midnightblue","k"] # Includes Glorys and different shade based on resolution
 

## Add Glorys


# Initial Variable Analysis -----
vnames          = ["sst","ssr","str","tx_sur","D20","Dmaxgrad"]
vunits          = [r"$\degree C$",r"$\frac{W}{m^2}$",r"$\frac{W}{m^2}$",r"$\frac{m}{s^2}$","m","m"]
vnames_long     = ["SST","Surface Shortwave","Surface Longwave","Zonal Wind Stress","Thermocline (20$\degree$ Isotherm)","Thermocline (Max Vertical Gradient)"]
vmaxes          = [2,40,20,0.02,20,20]


# ENSO Names -----
ninoname        = [r"$El$ $Ni\tilde{n}o$",r"$La$ $Ni\tilde{n}a$"]
ninocol         = ["cornflowerblue","firebrick"]
ninoshort       = ['nino','nina']

# Conversion for STR and SSR considering 3h Accumulation -----
conversion      = 1/(3*3600) # 3 h accumulation time...? #1/(24*30*3600)
# https://forum.ecmwf.int/t/surface-radiation-parameters-joule-m-2-to-watt-m-2/1588

# Bounding Boxes from Jin et al. 2020 Eqn. 6.6  -----
bbox_cep        = [150      , -130+360 , -5, 5]   # Central Equatorial Pacific, for [tau_x], 
bbox_nino3      = [-150+360 , -90+360  , -5, 5]  # Nino 3 Box: For SST, <tau_x>
bbox_nino34     = [-170+360 , -120+360 , -5, 5]
bbox_epac       = [-155+360 , -80+360  , -5, 5]  # Eastern Pacific (for h_e calculation)
bbox_wpac       = [120      , -155+360 , -5, 5]  # Western Pacific (for h_w calculation)
bbox_tropics    = [0        , 360      , -30,30] # Tropics (from Ceppi and Fueglistaler 2021)

bboxes      = [bbox_cep,bbox_nino3,bbox_nino34,bbox_epac,bbox_wpac]
bbnames_long = ["Central Equatorial Pacific","$Ni\tilde{n}o3$","$Ni\tilde{n}o3.4$","Tropical Eastern Pacific","Tropical Western Pacific"]
bbnames      = ["CEO","nino3","nino34","EPac","WPac"]


#%% CCF Variables

ccf_vars     = ["sst","eis","Tadv","r700","w700","ws10"]



