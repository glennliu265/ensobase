#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Ensemble Average for MESACLIP Output


For Each File...

(1) Load in Corresponding File for each Ensemble member
(2) Remove the mean seasonal cycle from each member (uses output from mesacliip_calc_scycle.py)
(3) Average across the members

Created on Tue Jul 21 12:26:49 2026

@author: gliu
"""


import numpy as np
import xarray as xr
import glob
import pandas as pd
import tqdm
import time

#%% Helper Functions

def cftime2str(times): # Copied from amv.proc
    "Convert array of cftime objects to string (YYYY-MM-DD)"
    newtimes = []
    for t in range(len(times)):
        newstr = "%04i-%02i-%02i" % (times[t].year,times[t].month,times[t].day)
        newtimes.append(newstr)
    return np.array(newtimes)


def generate_file_datestr(ystart=1920,yend=2006):
    # Works for the hi res CESM1 MESACLIP Simulations
    ii        = 0
    yend_file = 0
    datestring_byfile = []
    while yend_file < yend:
        
        ystart_file = ystart+5*ii
        yend_file   = ystart_file + 4
        
        
        tstart      = "%04i0101" % (ystart_file)
        if yend_file > yend: # [End on yend-01-01]
            tend        = "%04i0101" % (yend)
        else:
            tend        = "%04i1231" % (yend_file)
        datestr_loop = "%s-%s" % (tstart,tend)
        datestring_byfile.append(datestr_loop)
        ii += 1
    return datestring_byfile

    

#%% User Edits

# Number of Ensembles To Loop
vname   = "TS"
expname = 'hires'
if expname == 'hires':
    ensall  = np.arange(1,11,1)
else:
    ensall  = np.arange(1,23,1)
nens = len(ensall)
    
# Get List of Date Strings
datestrings = generate_file_datestr(ystart=1920,yend=2006)
nfiles      = len(datestrings)

# Get Paths
rawpath   = "/home/niu4/gliu8/share/CESM1/MESACLIP/RDA/hires/BHIST/day_1/regrid_1x1"

#%% Part (1), Get File Lists for each Ensemble Member

nclist_byens = []
for e in range(nens):
    ensnum   = ensall[e]
    ncsearch = "%s/ens%03i/%s_*_regrid1x1.nc" % (rawpath,ensnum,vname)
    nclist   = glob.glob(ncsearch)
    nclist.sort()
    nclist_byens.append(nclist)

# Debugging...
[print(nc[0]) for nc in nclist_byens]
[print(nc[-1]) for nc in nclist_byens]
[print(len(nc)) for nc in nclist_byens]

nfiles = len(nclist_byens[0])



#%% Load Seaspal Cycle for each File

scycle_byens = []
for e in tqdm.tqdm(range(nens)):

    ensnum   = ensall[e]
    ncsearch = "%s/ens%03i/%s_scycle_*.nc" % (rawpath,ensnum,vname)
    nclist   = glob.glob(ncsearch)
    nclist.sort()
    if len(nclist) > 1:
        print("Warning, more than 1 seasonal cycle was found for ens %03i" % (ensnum))
        print("\t Taking First One: %s" % nclist[0])
    
    ds = xr.open_dataset(nclist[0])[vname].load()
    scycle_byens.append(ds)
    

#%% Loop for each file (took 6 min 25 sec)


ens_avg_byfile = []
for ff in tqdm.tqdm(range(nfiles)):
    
    for e in tqdm.tqdm(range(nens)):
        
        # Load the File
        ensnum = ensall[e]
        ncfile = nclist_byens[e][ff] #"%s/ens%03i/%s_%s_regrid1x1.nc" % (rawpath,ensnum,vname,datestr)
        ds     = xr.open_dataset(ncfile)[vname].load()
        
        # Remove the seasonal cycle
        scycle_in  = scycle_byens[e]
        dsanom     = ds.groupby('time.dayofyear') - scycle_in
        
        # Accumulate by Ensemble Member
        if e == 0:
            ens_sum = dsanom.copy()
        else:
            ens_sum = ens_sum + dsanom
            
        del dsanom
        del ds
    
    # Divide by Number of Ensemble Members
    ens_sum = ens_sum / nens
    
    # Save the Output
    times           = cftime2str(ens_sum.time.data)
    tstart          = times[0].replace("-","")
    tend            = times[-1].replace("-","")
    outname         = "%s/ensavg/%s_Anomaly_EnsAvg_nens%0i_%s-%s.nc" % (rawpath,vname,nens,tstart,tend)
    ens_sum.to_netcdf(outname)
    
    
    
    # Append Output
    #ens_avg_byfile.append(ens_sum)

# #%% Next, concat by time

# ens_avg_cattime = xr.concat(ens_avg_byfile,dim='time')
# times           = cftime2str(ens_avg_cattime.time.data)
# tstart          = times[0].replace("-","")
# tend            = times[-1].replace("-","")
# st = time.time()
# outname         = "%s%s_Anomaly_EnsAvg_nens%0i_%s-%s.nc" % (rawpath,vname,nens,tstart,tend)
# ens_avg_cattime.to_netcdf(outname)
# print("Saved data in %.2fs" % (time.time()-st))
    






