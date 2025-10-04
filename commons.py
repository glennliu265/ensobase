#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 11:23:36 2025

@author: gliu
"""

#%% Shared variable names and experiment
# Copied from visualize composites

#datpath         = "/Users/gliu/Downloads/02_Research/01_Projects/07_ENSO/01_Data/TP_Crop/composites/"


# Simulation Names -----
expnames        = ["TCo319_ctl1950d","TCo319_ssp585","TCo1279-DART-1950","TCo1279-DART-2090"]
expnames_long   = ["31km Control","31km SSP585","9km 1950","9km 2090"]


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
bbox_cep        = [-150    , -130+360, -5, 5]  # Central Equatorial Pacific, for [tau_x], 
bbox_nino3      = [-150    , -90+360 , -5, 5]  # Nino 3 Box: For SST, <tau_x>
bbox_nino34     = [-170+360,-120+360,-5,5]
bbox_epac       = [-155    , -80+360 , -5, 5]  # Eastern Pacific (for h_e calculation)
bbox_wpac       = [120     , -155+360     , -5, 5]  # Western Pacific (for h_w calculation)