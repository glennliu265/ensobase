#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Seasonal Cycle for MESACLIP Output

Saves in the same ens folder with the filename <vname>_scycle_<tstart>_<tend>.nc

Created on Mon Jul 20 15:15:11 2026

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

#%% User Edits

# Number of Ensembles To Loop
ensall  = np.arange(1,23,1)

# Input File, variable, and Search String
vname     = "TS"
rawpath   = "/home/niu4/gliu8/share/CESM1/MESACLIP/RDA/lores/BHIST/day_1/regrid_1x1"
searchstr = "%s_*_regrid1x1.nc" % (vname)

#%% Start Loop
st = time.time()
nens  = len(ensall)
for e in tqdm.tqdm((range(nens))): # Loop by Ensemble
    
    ens       = ensall[e]
    searchstr = "%s/ens%03i/%s" % (rawpath,ens,searchstr)
    nclist    = glob.glob(searchstr)
    nclist.sort()
    nfiles    = len(nclist)
    print("Found %i files for ens%03i" % (nfiles,ens))
    
    times_byfile = [] # Loop and sum by File
    for ff in tqdm.tqdm(range(nfiles)): # Takes 26 Seconds (sum method), 2:24 (load all),
        # Load File
        nc = nclist[ff]
        ds = xr.open_dataset(nc)[vname].load()
        
        # Group by Day of year and sum
        dsday = ds.groupby('time.dayofyear')
        dssum = dsday.sum('time')
        
        # Delete to save memory
        times_byfile.append(ds.time)
        del ds
        

        # Get Number of Years (and check)
        nyrs         = np.array([len(dsday[np.int64(dd)].time)  for dd in np.arange(1,366,1)])
        equal_unique = (nyrs == np.unique(nyrs)[0]).sum().item()
        ngroups      = len(dsday.groups) # Number of Groups
        if ngroups != equal_unique:
            print("Warning, different number of years for certain days: %s" % (np.unique(nyrs)))
        del dsday
        
        # Append to running sum
        if ff == 0:
            ens_scycle = dssum.copy()
            ens_nyr    = np.array(nyrs).copy()
        else:
            ens_scycle = ens_scycle + dssum
            ens_nyr    = ens_nyr    + np.array(nyrs)
        # Dete to Save Memory
        del dssum
        # End File Loop
        
        
    # Method (2) (see mesaclip_calc_scycle_debug.py)
    scycle2 = ens_scycle / ens_nyr[:,None,None]
    
    # Save the Output
    tstart  = cftime2str(times_byfile[0].data)[0].replace("-","")#times_byfile[0][0].data.item()
    tend    = cftime2str(times_byfile[-1].data)[-1].replace("-","")#times_byfile[0][0].data.item()
    outname = "%s/ens%03i/%s_scycle_%s_%s.nc" % (rawpath,ens,vname,tstart,tend)
    scycle2.to_netcdf(outname)
   # End Ens loop
   
print("Script ran in %.2fs" % (time.time()-st))
        
        
        
    