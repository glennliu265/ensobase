#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Seasonal Cycle for MESACLIP Output (Debug Version)

Test 4 methods and calculate the differences

Created on Mon Jul 20 15:15:11 2026

@author: gliu

"""


import numpy as np
import xarray as xr
import glob
import pandas as pd
import tqdm
ensall  = np.arange(1,23,1)
rawpath = "/home/niu4/gliu8/share/CESM1/MESACLIP/RDA/lores/BHIST/day_1/regrid_1x1"
vname   = "TS"

def cftime2str(times):
    "Convert array of cftime objects to string (YYYY-MM-DD)"
    newtimes = []
    for t in range(len(times)):
        newstr = "%04i-%02i-%02i" % (times[t].year,times[t].month,times[t].day)
        newtimes.append(newstr)
    return np.array(newtimes)



debug = True

nens  = len(ensall)

for e in range(nens):
    
    ens       = ensall[e]
    searchstr = "%s/ens%03i/%s_*_regrid1x1.nc" % (rawpath,ens,vname)
    nclist    = glob.glob(searchstr)
    nclist.sort()
    nfiles    = len(nclist)
    print("Found %i files for ens%03i" % (nfiles,ens))
    
    #ens_nyr    = []
    if debug:
        dsall       = []
    for ff in tqdm.tqdm(range(nfiles)): # Takes 26 Seconds (sum method), 2:24 (load all),
        # Load File
        nc = nclist[ff]
        ds = xr.open_dataset(nc)[vname].load()
        
        # Group by Day of year and sum
        dsday = ds.groupby('time.dayofyear')
        dssum = dsday.sum('time')
        
        # Delete to save memory
        if debug:
            dsall.append(ds) # Check Memory Intensive CAse
        else:
            del ds
        
        # Fix Time from cftime.datetimeNoLeap to String
        #timesnew = cftime2str(ds.time.data)
        #timesnew = pd.to_datetime(timesnew)
        #ds['time'] = timesnew
        
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
    
#%% Some Debug Checks

"""
Evaluate some different methods

1. Concat all and take seasonal cycle
2. Sum all then divide by nyrs
3. Take file-wise average, then do weighted sum


"""

# Method (1)
dsallcat   = xr.concat(dsall,dim='time')
scycle1 = dsallcat.groupby('time.dayofyear').mean('time')

# Method (2)
scycle2 = ens_scycle / ens_nyr[:,None,None]

# Method (3)
def get_nyr(ds):
    ds = ds.groupby('time.dayofyear')
    nyr = np.array([len(ds[np.int64(dd)].time)  for dd in np.arange(1,366,1)])
    return nyr
scyclebyfile  = [ds.groupby('time.dayofyear').mean('time') for ds in dsall]
nyr_perfile   = np.array([get_nyr(ds) for ds in dsall])
nyr_byday     = nyr_perfile.sum(0) #SUm along File Dimension, Total Count per day
wgts_byfile   = nyr_perfile / nyr_byday
scyclebyfile  = xr.concat(scyclebyfile,dim='file')
scycle3       = scyclebyfile * wgts_byfile[:,:,None,None]
scycle3       = scycle3.sum('file')

#filewise_scycle = [ds.groupby('time.dayofyear').mean('time') for ds in dsall]

# Method (4) Load CDO Version
scycle4 = xr.open_dataset("/home/niu4/gliu8/share/CESM1/MESACLIP/RDA/lores/BHIST/day_1/regrid_1x1/ens001/ens001_merged_scycle.nc").load()




#%% 

# Just Checked all of these methods, they seem equivalent!


lonf = 330
latf = 50
import matplotlib.pyplot as plt

fig,ax = plt.subplots(1,1)

plotvar = scycle1
ax.plot(plotvar.dayofyear,plotvar.sel(lon=lonf,lat=latf,method='nearest'),label="Method 1")

plotvar = scycle2
ax.plot(plotvar.dayofyear,plotvar.sel(lon=lonf,lat=latf,method='nearest'),label="Method 2",ls='dashed')

plotvar = scycle3
ax.plot(plotvar.dayofyear,plotvar.sel(lon=lonf,lat=latf,method='nearest'),label="Method 3",ls='dotted',c='k')
doy = plotvar.dayofyear

plotvar = scycle4.TS
ax.plot(doy,plotvar.sel(lon=lonf,lat=latf,method='nearest'),label="Method 4 (CDO)",ls='dashdot',c='red')


ax.legend()
            
            
        
        
        
    