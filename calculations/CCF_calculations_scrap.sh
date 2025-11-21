





# ----------------------------------------------------------------------------------------------------
# Recalculate ucc (upper cloud content) for each simulation 
# ----------------------------------------------------------------------------------------------------

# 5km Control Simulation, calculate [TCC - LCC]
cdo -sub /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo2559-DART-1950C_5km/atm/TCo2559-DART-1950C_atm_10256x5120_1m_tcc_1950-1959.nc /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo2559-DART-1950C_5km/atm/TCo2559-DART-1950C_atm_10256x5120_1m_lcc_1950-1959.nc /home/niu4/gliu8/projects/scrap/processed_global/temp.nc 
cdo chname,tcc,ucc /home/niu4/gliu8/projects/scrap/processed_global/temp.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo2559-DART-1950C_ucc.nc

# 9km Control Simulation, calculate [TCC - LCC]
cdo -sub /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/TCo1279-DART-1950_atm_remapped_1m_tcc_240months.nc /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/TCo1279-DART-1950_atm_remapped_1m_lcc_240months.nc /home/niu4/gliu8/projects/scrap/processed_global/temp.nc 
cdo chname,tcc,ucc /home/niu4/gliu8/projects/scrap/processed_global/temp.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_ucc.nc

# Now, regrid each case
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo2559-DART-1950C_ucc.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo2559-DART-1950C_ucc_regrid1x1.nc
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_ucc.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo1279-DART-1950_ucc_regrid1x1.nc


# ----------------------------------------------------------------------------------------------------
# Wind Stress: Simply Regrid 
# ----------------------------------------------------------------------------------------------------

# Wind Stress 9km Control, 1950
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc  /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/ocn/mon/TCo1279_1950ctl_ty_surf_1m_1950-1969.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo1279-DART-1950_ty_sur_regrid1x1.nc
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc  /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/ocn/mon/TCo1279_1950ctl_tx_surf_1m_1950-1969.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo1279-DART-1950_tx_sur_regrid1x1.nc

# Wind Stress 5km Control, 1950
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc  /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo2559-DART-1950C_5km/ocn/TCo2559-DART-1950C_ocn_remapped_3600x1800_1m_tx_sur_1950-1959.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo2559-DART-1950C_tx_sur_regrid1x1.nc
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc  /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo2559-DART-1950C_5km/ocn/TCo2559-DART-1950C_ocn_remapped_3600x1800_1m_ty_sur_1950-1959.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo2559-DART-1950C_ty_sur_regrid1x1.nc


# ----------------------------------------------------------------------------------------------------
# Regrid SST
# ----------------------------------------------------------------------------------------------------


cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo2559-DART-1950C_5km/ocn/TCo2559-DART-1950C_ocn_remapped_3600x1800_1m_sst_1950-1958.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo2559-DART-1950C_sst_regrid1x1.nc

cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/TCo1279-DART-1950_atm_remapped_1m_sst_240months.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo1279-DART-1950_sst_regrid1x1.nc



# ----------------------------------------------------------------------------------------------------
# Regrid w700
# ----------------------------------------------------------------------------------------------------

# 9km 
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/TCo1279-DART-1950_atm_remapped_1m_w700_240months.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/temp.nc
cdo chname,w,w700 /home/niu4/gliu8/projects/scrap/regrid_1x1/temp.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo1279-DART-1950_w700_regrid1x1.nc

# 5km
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo2559-DART-1950C_5km/atm/TCo2559-DART-1950C_atm_10256x5120_1m_pl_w_700_1950-1959_jan1950missing.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/temp.nc
cdo chname,w,w700 /home/niu4/gliu8/projects/scrap/regrid_1x1/temp.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo2559-DART-1950C_w700_regrid1x1.nc


# ----------------------------------------------------------------------------------------------------
# Regrid r700
# ----------------------------------------------------------------------------------------------------

# 9km
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/TCo1279-DART-1950_atm_remapped_1m_r700_240months.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/temp.nc
cdo chname,r,r700 /home/niu4/gliu8/projects/scrap/regrid_1x1/temp.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo1279-DART-1950_r700_regrid1x1.nc

# 5km
cdo remapbil,U10/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo2559-DART-1950C_5km/atm/TCo2559-DART-1950C_atm_10256x5120_1m_pl_r_700_1950-1959_jan1950missing_feb1950insteadtwice.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/temp.nc
cdo chname,r,r700 /home/niu4/gliu8/projects/scrap/regrid_1x1/temp.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo2559-DART-1950C_r700_regrid1x1.nc



#cdo -chname,w,w700 /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo1279-DART-1950_w700_regrid1x1.nc remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/TCo1279-DART-1950_atm_remapped_1m_w700_240months.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/temp.nc
