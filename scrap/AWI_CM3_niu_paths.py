#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 14:19:57 2025

@author: gliu
"""



import glob

expname = "TCo2559-DART-1950C"
vname   = "tx_sur"
ensnum  = None

get_rawpath_awi(expname,vname,ensnum)


def get_rawpath_awi(expname,vname,ensnum=None):
    # TCo319_ctl1950d
    ctlpath0_31= "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ctl1950d/monthly/"
    ctlpath1_31= "//export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_control/"

    # TCo319_ssp585
    ssppath0_31="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/"
    ssppath0_31_tropics="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/tropics_only/"
    if ensnum is not None:
        ssppath0_31_ens = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585_ens%02i/"
    
    # TCo1279_DART-1950
    ctlpath0_09_atm = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/"
    ctlpath0_09_ocn = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/ocn/mon/"

    # TCo1279_DART-2090
    ssppath0_09     = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/"

    # TCo2559-DART-1950C
    ctlpath0_05_atm = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo2559-DART-1950C_5km/atm/"
    ctlpath0_05_ocn = "/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo2559-DART-1950C_5km/ocn/"


    
    notfound=False
    # ========================================================
    if expname == "TCo319_ctl1950d": # 31 km
        # First Try the Control Path (post 1950 crops)
        searchstr = "%s/*_%s_*.nc" % (ctlpath0_31,vname)
        nclist    = glob.glob(searchstr)
        if len(nclist) > 0:
            print("\tFound in [TCo319_ctl1950d] path!")
            #continue
        else:
            print("\tNot found in [TCo319_ctl1950d] path...")
            # Next Try the older Control Path
            searchstr = "%s/*_%s_*.nc" % (ctlpath1_31,vname)
            nclist    = glob.glob(searchstr)
            if len(nclist) > 0:
                print("\tFound in [TCo319_control] path!")
            else:
                notfound=True
    # ========================================================
    elif expname == "TCo319_ssp585":
        if ensnum is not None:
            searchstr = ssppath0_31_ens % (ensnum) + "*_%s_*.nc" % (vname)
            nclist    = glob.glob(searchstr)
            print("\tFound in [TCo319_ssp585_ens%02i]!" % ensnum)
            #continue
        else:
            print("\tVariable not available for ens members... looking in original folder")
        
        # First Try the Control Path (post 1950 crops)
        searchstr = "%s/*_%s_*.nc" % (ssppath0_31,vname)
        nclist    = glob.glob(searchstr)
        if len(nclist) > 0:
            print("\tFound in [TCo319_ssp585] path!")
            #continue
        else: # Otherwise Try Tropics Only
            print("\tNot found in [TCo319_ssp585] path...")
            # First Try the Control Path (post 1950 crops)
            searchstr = "%s/*_%s_*.nc" % (ssppath0_31_tropics,vname)
            nclist    = glob.glob(searchstr)
            if len(nclist) > 0:
                print("\tFound in [TCo319_ssp585_tropics_only] path!")
                #continue
            else:
                notfound=True
    # ========================================================
    elif expname == "TCo1279-DART-1950":
        # First Try the ATM Path (post 1950 crops)
        searchstr = "%s/*_%s_*.nc" % (ctlpath0_09_atm,vname)
        nclist    = glob.glob(searchstr)
        if len(nclist) > 0:
            print("\tFound in [TCo1279-DART-1950_atm] path!")
            #continue
        else:
            print("\tNot found in [TCo1279-DART-1950_atm] path...")
            # First Try the ATM Path (post 1950 crops)
            searchstr = "%s/*_%s_*.nc" % (ctlpath0_09_ocn,vname)
            nclist    = glob.glob(searchstr)
            if len(nclist) > 0:
                print("\tFound in [TCo1279-DART-1950_ocn] path!")
                #continue
            else:
                notfound=True
    # ========================================================
    elif expname == "TCo1279-DART-2090":
        # First Try the Control Path (post 1950 crops)
        searchstr = "%s/*_%s_*.nc" % (ssppath0_09,vname)
        nclist    = glob.glob(searchstr)
        if len(nclist) > 0:
            print("\tFound in [TCo319_ctl1950d] path!")
        else:
            notfound=True
    # ========================================================
    elif expname == "TCo2559-DART-1950C":
        # First Try the ATM Path (post 1950 crops)
        searchstr = "%s/*_%s_*.nc" % (ctlpath0_05_atm,vname)
        nclist    = glob.glob(searchstr)
        if len(nclist) > 0:
            print("\tFound in [TCo2559-DART-1950C_atm] path!")
            #continue
        else:
            print("\tNot found in [TCo2559-DART-1950C_atm] path...")
            # First Try the ATM Path (post 1950 crops)
            searchstr = "%s/*_%s_*.nc" % (ctlpath0_05_ocn,vname)
            nclist    = glob.glob(searchstr)
            if len(nclist) > 0:
                print("\tFound in [TCo2559-DART-1950C_ocn] path!")
                #continue
            else:
                notfound=True
    # ========================================================
    else:
        print("%s not recognized..." % expname)
        notfound=True
    if notfound:
        print("%s Not found for %s" % (vname,expname))
        return None
    print(nclist)
    return None
    
    
            
        

            
        
        
    # ========================================================

    
        
        
            
        
            
            
        

        
            
            
            
    




