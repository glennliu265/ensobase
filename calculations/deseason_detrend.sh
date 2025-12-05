




# Calculate Seasonal Cycle
cdo ymonmean /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo2559-DART-1950C_sst_regrid1x1.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/scycle/TCo2559-DART-1950C_sst.nc
# Subtract
cdo ymonsub /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo2559-DART-1950C_sst_regrid1x1.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/scycle/TCo2559-DART-1950C_sst.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/temp.nc
# Now Detrend (Linear)
cdo detrend /home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/temp.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/TCo2559-DART-1950C_sst.nc

expname='TCo2559-DART-1950C'
vname='sst'
infile='TCo2559-DART-1950C_sst_regrid1x1.nc'
rawpath='/home/niu4/gliu8/projects/scrap/regrid_1x1/'

scyclepath='/home/niu4/gliu8/projects/scrap/regrid_1x1/scycle/'
detrendpath='/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/'




# Anomalize and Detrend (linearly) the SST for original resolution data
# 31km Control ------------------
expname='TCo319_ctl1950d'
vname='sst'
infile='TCo319_ctl1950d_sst.nc'
rawpath='/home/niu4/gliu8/projects/scrap/processed_global'
scyclepath='/home/niu4/gliu8/projects/scrap/processed_global/scycle'
detrendpath='/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1'
cdo ymonmean ${rawpath}/${infile} ${scyclepath}/${expname}_${vname}.nc
cdo ymonsub ${rawpath}/${infile} ${scyclepath}/${expname}_${vname}.nc ${detrendpath}/temp.nc
cdo detrend ${detrendpath}/temp.nc ${detrendpath}/${expname}_${vname}.nc

# 31km Future ------------------
expname='TCo319_ssp585'
vname='sst'
infile='TCo319_ssp585_sst.nc'
rawpath='/home/niu4/gliu8/projects/scrap/processed_global/'
scyclepath='/home/niu4/gliu8/projects/scrap/processed_global/scycle/'
detrendpath='/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1/'
cdo ymonmean ${rawpath}/${infile} ${scyclepath}/${expname}_${vname}.nc
cdo ymonsub ${rawpath}/${infile} ${scyclepath}/${expname}_${vname}.nc ${detrendpath}/temp4.nc
cdo detrend ${detrendpath}/temp4.nc ${detrendpath}/${expname}_${vname}.nc

# 9km Control ---------------
expname='TCo1279-DART-1950'
vname='sst'
infile='TCo1279-DART-1950_atm_remapped_1m_sst_240months.nc'
rawpath='/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/'
scyclepath='/home/niu4/gliu8/projects/scrap/processed_global/scycle/'
detrendpath='/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1/'
cdo ymonmean ${rawpath}/${infile} ${scyclepath}/${expname}_${vname}.nc
cdo ymonsub ${rawpath}/${infile} ${scyclepath}/${expname}_${vname}.nc ${detrendpath}/temp3.nc
cdo detrend ${detrendpath}/temp3.nc ${detrendpath}/${expname}_${vname}.nc

# 9km Future ---------------
expname='TCo1279-DART-2090'
vname='sst'
infile='TCo1279_2090slice_aleph_sst_1m_2090-2099.nc'
rawpath='/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/'
scyclepath='/home/niu4/gliu8/projects/scrap/processed_global/scycle/'
detrendpath='/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1/'
cdo ymonmean ${rawpath}/${infile} ${scyclepath}/${expname}_${vname}.nc
cdo ymonsub ${rawpath}/${infile} ${scyclepath}/${expname}_${vname}.nc ${detrendpath}/temp2.nc
cdo detrend ${detrendpath}/temp2.nc ${detrendpath}/${expname}_${vname}.nc


# 5km Control --------------
expname='TCo2559-DART-1950C'
vname='sst'
infile='TCo2559-DART-1950C_sst.nc'
rawpath='/home/niu4/gliu8/projects/scrap/processed_global/'
scyclepath='/home/niu4/gliu8/projects/scrap/processed_global/scycle/'
detrendpath='/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1/'
cdo ymonmean ${rawpath}/${infile} ${scyclepath}/${expname}_${vname}.nc
cdo ymonsub ${rawpath}/${infile} ${scyclepath}/${expname}_${vname}.nc ${detrendpath}/temp1.nc
cdo detrend ${detrendpath}/temp1.nc ${detrendpath}/${expname}_${vname}.nc





