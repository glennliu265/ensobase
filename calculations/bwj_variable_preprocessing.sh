


# 31 km Control

# ```
# TCo319_ctl1950d
# TCo319_ssp585
# tx_sur
# /home/niu4/gliu8/projects/scrap/processed_global/
# /home/niu4/gliu8/projects/scrap/regrid_1x1/
# ```

# ---------------
# Zonal Wind Stress tx_sur
# ---------------

# Copy Zonal Wind Stress Over to folder (31km Future)
cp /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/TCo319_ssp585_tx_sur_1m_2015-2100_1x1regrid.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo319_ssp585_tx_sur_regrid1x1.nc

# Copy Zonal Wind Stress Over to folder (31km Control)
cp /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ctl1950d/monthly/awicm3_tco319_ctl1950d_tx_sur_fesom_r360x180_1950-2100.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo319_ctl1950d_tx_sur_regrid1x1.nc

# ---------------
# Tsub
# Slice at 75 meters
# Rename to temp75
# ---------------

# Try to select 75 meters (31km Control)
cdo sellevel,75 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ctl1950d/monthly/awicm3_tco319_ctl1950d_temp_fesom_r360x180_1950-2100.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel.nc
cdo chname,temp,temp75 /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo319_ctl1950d_temp75_regrid1x1.nc
rm templevel.nc

# Try to select 75 meters (31km SSP5856)
cdo sellevel,75 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/TCo319_ssp585_temp_1m_2015-2100_1x1regrid.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel.nc
cdo chname,temp,temp75 /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo319_ssp585_temp75_regrid1x1.nc
rm templevel.nc

# ---------------
# U, V 
# (select to 50 meters, then take depth weighted mean) ----
# Rename to uvel_ML50, vvel_ML50
# Note that depths are 15,25,35,45, midpoints of first 5 layers. I thus just took the *unweighted* mean...
# ---------------

# U 31km Control
cdo sellevel,0/50 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ctl1950d/monthly/awicm3_tco319_ctl1950d_u_fesom_r360x180_1950-2100.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel.nc
cdo vertmean /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel_mean.nc
cdo chname,u,uvel_ML50 /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel_mean.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo319_ctl1950d_uvel_ML50_regrid1x1.nc
#cdo chname,u,uML50 /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel_mean.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo319_ctl1950d_uML50.nc
#mv TCo319_ctl1950d_uML50.nc TCo319_ctl1950d_uML50_regrid1x1.nc
# Change the Name
#cdo chname,uML50,uvel_ML50 TCo319_ctl1950d_uML50_regrid1x1.nc TCo319_ctl1950d_uvel_ML50_regrid1x1.nc
#rm TCo319_ctl1950d_uML50_regrid1x1.nc

# ...

# V 31km Control
cdo sellevel,0/50 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ctl1950d/monthly/awicm3_tco319_ctl1950d_v_fesom_r360x180_1950-2100.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel.nc
cdo vertmean /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel_mean.nc
#cdo chname,v,vML50 /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel_mean.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo319_ctl1950d_vML50.nc
# Change the Name
#mv TCo319_ctl1950d_vML50.nc TCo319_ctl1950d_vML50_regrid1x1.nc
#cdo chname,vML50,vvel_ML50 TCo319_ctl1950d_vML50_regrid1x1.nc TCo319_ctl1950d_vvel_ML50_regrid1x1.nc

# ...

# U 31km Future
cdo sellevel,0/50 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/TCo319_ssp585_u_1m_2015-2100_1x1regrid.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel.nc
cdo vertmean /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel_mean.nc
#cdo chname,u,uML50 /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel_mean.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo319_ssp585_uML50.nc
#rm "templevel*.nc"
# Change the Name
#cdo chname,uML50,uvel_ML50 TCo319_ssp585_uML50_regrid1x1.nc TCo319_ssp585_uvel_ML50_regrid1x1.nc

# V 31km Future
cdo sellevel,0/50 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/TCo319_ssp585_v_1m_2015-2100_1x1regrid.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel.nc
cdo vertmean /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel_mean.nc
#cdo chname,v,vML50 /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel_mean.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo319_ssp585_vML50.nc
#rm "templevel*.nc"
#mv TCo319_ssp585_uML50.nc TCo319_ssp585_uML50_regrid1x1.nc
#mv TCo319_ssp585_vML50.nc TCo319_ssp585_vML50_regrid1x1.nc
# Change the Name
#cdo chname,vML50,vvel_ML50 TCo319_ssp585_vML50_regrid1x1.nc TCo319_ssp585_vvel_ML50_regrid1x1.nc

# ---------------
# W, select at 50 meters
# rename to wvel50
# ---------------

# w, 31km Control
cdo sellevel,50 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ctl1950d/monthly/awicm3_tco319_ctl1950d_w_fesom_r360x180_1950-2100.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel.nc
cdo chname,w,wvel50 /home/niu4/gliu8/projects/scrap/regrid_1x1/templevel.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo319_ctl1950d_wvel50_regrid1x1.nc
rm templevel.nc

# ---------------
# SSH 
# ---------------

# Copy over SSH (TCo319 Control)
cp /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ctl1950d/monthly/awicm3_tco319_ctl1950d_ssh_fesom_r360x180_1950-2100.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo319_ctl1950d_ssh_regrid1x1.nc

# Copy over SSH (TCo319 SSP)
cp /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/TCo319_ssp585_ssh_1m_2015-2117_1x1regrid.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo319_ssp585_ssh_regrid1x1.nc

# ===================================================
# Now Preprocess for variables (anomalize and detrend)
# Copied from anom_detrend1_awiloop
# ===================================================
dpath="/home/niu4/gliu8/projects/scrap/regrid_1x1"
detrendpath="/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1"
scyclepath="/home/niu4/gliu8/projects/scrap/regrid_1x1/scycle"
vnames=("tx_sur" "uvel_ML50" "vvel_ML50" "temp75" "ssh")
expnames=("TCo319_ssp585" "TCo319_ctl1950d")
for exp in ${expnames[@]}; do # This loop format is for zsh. Use ${expes[@]} if you are using bash.
    for vname in ${vnames[@]}; do

        infile=${dpath}/${exp}_${vname}_regrid1x1.nc
        scyclefile=${scyclepath}/${exp}_${vname}_regrid1x1.nc
        outfile=${detrendpath}/${exp}_${vname}_regrid1x1.nc

        cdo ymonmean ${infile} ${scyclefile}
        cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
        cdo detrend ${detrendpath}/temp1.nc ${outfile}
        echo "Completed $vname for $exp"

    done
done




