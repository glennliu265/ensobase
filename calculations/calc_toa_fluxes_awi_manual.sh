

# Control 31km
expname="TCo319_ctl1950d"
rawpath="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_control"
outpath="/home/niu4/gliu8/projects/scrap/processed_global"
#ncstr="ctl1950d_atm_remapped_1m_tsr_1850-2134.nc"

# All Sky
cdo -add ${rawpath}/ctl1950d_atm_remapped_1m_tsr_1850-2134.nc  ${rawpath}/ctl1950d_atm_remapped_1m_ttr_1850-2134.nc ${outpath}/temp.nc
cdo chname,tsr,allsky ${outpath}/temp.nc ${outpath}/${expname}_allsky.nc

# Clear Sky
cdo -add ${rawpath}/ctl1950d_atm_remapped_1m_tsrc_1850-2134.nc  ${rawpath}/ctl1950d_atm_remapped_1m_ttrc_1850-2134.nc ${outpath}/temp.nc
cdo chname,tsrc,clearsky ${outpath}/temp.nc ${outpath}/${expname}_clearsky.nc

# CRE
cdo -sub ${outpath}/${expname}_allsky.nc  ${outpath}/${expname}_clearsky.nc ${outpath}/temp.nc
cdo chname,allsky,cre ${outpath}/temp.nc ${outpath}/${expname}_cre.nc

# TTCRE 
cdo -sub ${rawpath}/ctl1950d_atm_remapped_1m_ttr_1850-2134.nc  ${rawpath}/ctl1950d_atm_remapped_1m_ttrc_1850-2134.nc ${outpath}/temp.nc
cdo chname,ttr,ttcre ${outpath}/temp.nc ${outpath}/${expname}_ttcre.nc

# TSCRE
cdo -sub ${rawpath}/ctl1950d_atm_remapped_1m_tsr_1850-2134.nc  ${rawpath}/ctl1950d_atm_remapped_1m_tsrc_1850-2134.nc ${outpath}/temp.nc
cdo chname,tsr,tscre ${outpath}/temp.nc ${outpath}/${expname}_tscre.nc

# 1950 9km
expname="TCo1279-DART-1950"
rawpath="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon"
outpath="/home/niu4/gliu8/projects/scrap/processed_global"
cdo -sub ${rawpath}/TCo1279-DART-1950_atm_remapped_1m_ttr_240months.nc ${rawpath}/TCo1279-DART-1950_atm_remapped_1m_ttrc_240months.nc ${outpath}/temp.nc
cdo chname,ttr,ttcre ${outpath}/temp.nc ${outpath}/${expname}_ttcre.nc


# ----- Runnning


# SSP585 31km
expname="TCo319_ssp585"
rawpath="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585"
outpath="/home/niu4/gliu8/projects/scrap/processed_global"

# All Sky
cdo -add ${rawpath}/TCo319_ssp585_tsr_1m_2015-2114_atmogrid.nc  ${rawpath}/TCo319_ssp585_ttr_1m_2015-2114_atmogrid.nc ${outpath}/temp.nc
cdo chname,tsr,allsky ${outpath}/temp.nc ${outpath}/${expname}_allsky.nc

# Clear Sky
cdo -add ${rawpath}/TCo319_ssp585_tsrc_1m_2015-2114_atmogrid.nc  ${rawpath}/TCo319_ssp585_ttrc_1m_2015-2114_atmogrid.nc ${outpath}/temp.nc
cdo chname,tsrc,clearsky ${outpath}/temp.nc ${outpath}/${expname}_clearsky.nc

# CRE
cdo -sub ${outpath}/${expname}_allsky.nc  ${outpath}/${expname}_clearsky.nc ${outpath}/temp.nc
cdo chname,allsky,cre ${outpath}/temp.nc ${outpath}/${expname}_cre.nc

# TTCRE 
cdo -sub ${rawpath}/TCo319_ssp585_ttr_1m_2015-2114_atmogrid.nc  ${rawpath}/TCo319_ssp585_ttrc_1m_2015-2114_atmogrid.nc ${outpath}/temp.nc
cdo chname,ttr,ttcre ${outpath}/temp.nc ${outpath}/${expname}_ttcre.nc

# TSCRE
cdo -sub ${rawpath}/TCo319_ssp585_tsr_1m_2015-2114_atmogrid.nc  ${rawpath}/TCo319_ssp585_tsrc_1m_2015-2114_atmogrid.nc ${outpath}/temp.nc
cdo chname,tsr,tscre ${outpath}/temp.nc ${outpath}/${expname}_tscre.nc



# 1950 5km
expname="TCo2559-DART-1950C"
rawpath="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo2559-DART-1950C_5km/atm"
outpath="/home/niu4/gliu8/projects/scrap/processed_global"

# All Sky
cdo -add ${rawpath}/TCo2559-DART-1950C_atm_10256x5120_1m_tsr_1950-1959.nc  ${rawpath}/TCo2559-DART-1950C_atm_10256x5120_1m_ttr_1950-1959.nc ${outpath}/temp.nc
cdo chname,tsr,allsky ${outpath}/temp.nc ${outpath}/${expname}_allsky.nc

# Clear Sky
cdo -add ${rawpath}/TCo2559-DART-1950C_atm_10256x5120_1m_tsrc_1950-1959.nc  ${rawpath}/TCo2559-DART-1950C_atm_10256x5120_1m_ttrc_1950-1959.nc ${outpath}/temp.nc
cdo chname,tsrc,clearsky ${outpath}/temp.nc ${outpath}/${expname}_clearsky.nc

# CRE
cdo -sub ${outpath}/${expname}_allsky.nc  ${outpath}/${expname}_clearsky.nc ${outpath}/temp.nc
cdo chname,allsky,cre ${outpath}/temp.nc ${outpath}/${expname}_cre.nc

# TTCRE 
cdo -sub ${rawpath}/TCo2559-DART-1950C_atm_10256x5120_1m_ttr_1950-1959.nc  ${rawpath}/TCo2559-DART-1950C_atm_10256x5120_1m_ttrc_1950-1959.nc ${outpath}/temp.nc
cdo chname,ttr,ttcre ${outpath}/temp.nc ${outpath}/${expname}_ttcre.nc

# TSCRE
cdo -sub ${rawpath}/TCo2559-DART-1950C_atm_10256x5120_1m_tsr_1950-1959.nc  ${rawpath}/TCo2559-DART-1950C_atm_10256x5120_1m_tsrc_1950-1959.nc ${outpath}/temp.nc
cdo chname,tsr,tscre ${outpath}/temp.nc ${outpath}/${expname}_tscre.nc

# ---- Not Yet Run

# 2090 9km
expname="TCo1279-DART-2090"
rawpath="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km"
outpath="/home/niu4/gliu8/projects/scrap/processed_global"

# All Sky
cdo -add ${rawpath}/TCo1279_2090slice_aleph_tsr_1m_2090-2099.nc  ${rawpath}/TCo1279_2090slice_aleph_ttr_1m_2090-2099.nc ${outpath}/temp.nc
cdo chname,tsr,allsky ${outpath}/temp.nc ${outpath}/${expname}_allsky.nc

# Clear Sky
cdo -add ${rawpath}/TCo1279_2090slice_aleph_tsrc_1m_2090-2099.nc  ${rawpath}/TCo1279_2090slice_aleph_ttrc_1m_2090-2099.nc ${outpath}/temp.nc
cdo chname,tsrc,clearsky ${outpath}/temp.nc ${outpath}/${expname}_clearsky.nc

# CRE
cdo -sub ${outpath}/${expname}_allsky.nc  ${outpath}/${expname}_clearsky.nc ${outpath}/temp.nc
cdo chname,allsky,cre ${outpath}/temp.nc ${outpath}/${expname}_cre.nc

# TTCRE 
cdo -sub ${rawpath}/TCo1279_2090slice_aleph_ttr_1m_2090-2099.nc  ${rawpath}/TCo1279_2090slice_aleph_ttrc_1m_2090-2099.nc ${outpath}/temp.nc
cdo chname,ttr,ttcre ${outpath}/temp.nc ${outpath}/${expname}_ttcre.nc

# TSCRE
cdo -sub ${rawpath}/TCo1279_2090slice_aleph_tsr_1m_2090-2099.nc  ${rawpath}/TCo1279_2090slice_aleph_tsrc_1m_2090-2099.nc ${outpath}/temp.nc
cdo chname,tsr,tscre ${outpath}/temp.nc ${outpath}/${expname}_tscre.nc




#%% Rerun for 31 km

# Calculate All Sky
cdo -add /home/niu4/gliu8/projects/common_data/awi_cm3/Tco319_ctl1950d_tsr_1950-2100.nc /home/niu4/gliu8/projects/common_data/awi_cm3/Tco319_ctl1950d_ttr_1950-2100.nc /home/niu4/gliu8/projects/scrap/processed_global/temp.nc
cdo chname,tsr,allsky /home/niu4/gliu8/projects/scrap/processed_global/temp.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_allsky.nc

# Recalculate Clear Sky
cdo -add /home/niu4/gliu8/projects/common_data/awi_cm3/Tco319_ctl1950d_tsrc_1950-2100.nc /home/niu4/gliu8/projects/common_data/awi_cm3/Tco319_ctl1950d_ttr_1950-2100.nc /home/niu4/gliu8/projects/scrap/processed_global/temp.nc
cdo chname,tsrc,clearsky /home/niu4/gliu8/projects/scrap/processed_global/temp.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_clearsky.nc

# Recalculate CRE
cdo -sub /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_allsky.nc  /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_clearsky.nc /home/niu4/gliu8/projects/scrap/processed_global/temp.nc
cdo chname,allsky,cre /home/niu4/gliu8/projects/scrap/processed_global/temp.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_cre.nc


/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo2559-DART-1950C_5km/atm






# Calculate for ERA5 # ============================================================
# expname="TCo2559-DART-1950C"
rawpath='/home/niu4/gliu8/share/ERA5/processed'
outpath='/home/niu4/gliu8/share/ERA5/processed'

# All Sky
cdo -add ${rawpath}/tsr_1979_2024.nc  ${rawpath}/ttr_1979_2024.nc ${outpath}/temp.nc
cdo chname,tsr,allsky ${outpath}/temp.nc ${outpath}/allsky_1979_2024.nc

# Clear Sky
cdo -add ${rawpath}/tsrc_1979_2024.nc  ${rawpath}/ttrc_1979_2024.nc ${outpath}/temp.nc
cdo chname,tsrc,clearsky ${outpath}/temp.nc ${outpath}/clearsky_1979_2024.nc

# CRE
cdo -sub ${outpath}/allsky_1979_2024.nc  ${outpath}/clearsky_1979_2024.nc ${outpath}/temp.nc
cdo chname,allsky,cre ${outpath}/temp.nc ${outpath}/cre_1979_2024.nc

# TTCRE 
cdo -sub ${outpath}/ttr_1979_2024.nc  ${outpath}/ttrc_1979_2024.nc ${outpath}/temp.nc
cdo chname,ttr,ttcre ${outpath}/temp.nc ${outpath}/ttcre_1979_2024.nc

# TSCRE 
#rawpath='/home/niu4/gliu8/share/ERA5/processed'
#outpath='/home/niu4/gliu8/share/ERA5/processed'
cdo -sub ${outpath}/tsr_1979_2024.nc  ${outpath}/tsrc_1979_2024.nc ${outpath}/temp.nc
cdo chname,tsr,tscre ${outpath}/temp.nc ${outpath}/tscre_1979_2024.nc

# ==================================================================================

# Calculate for CERES-EBAF =========================================================
rawpath='/home/niu4/gliu8/share/CERES/processed'
outpath='/home/niu4/gliu8/share/CERES/processed'
# CRE
cdo -sub ${rawpath}/CERES_EBAF_allsky_2000-03_to_2025-08.nc  ${rawpath}/CERES_EBAF_clearsky_2000-03_to_2025-08.nc ${outpath}/temp.nc
cdo chname,allsky,cre ${outpath}/temp.nc ${outpath}/CERES_EBAF_cre_2000-03_to_2025-08.nc

# TTCRE
rawpath='/home/niu4/gliu8/share/CERES/processed'
outpath='/home/niu4/gliu8/share/CERES/processed'
cdo -sub ${rawpath}/CERES_EBAF_ttr_2000-03_to_2025-08.nc  ${rawpath}/CERES_EBAF_ttrc_2000-03_to_2025-08.nc ${outpath}/temp.nc
cdo chname,ttr,ttcre ${outpath}/temp.nc ${outpath}/CERES_EBAF_ttcre_2000-03_to_2025-08.nc

# TSCRE
rawpath='/home/niu4/gliu8/share/CERES/processed'
outpath='/home/niu4/gliu8/share/CERES/processed'
cdo -sub ${rawpath}/CERES_EBAF_tsr_2000-03_to_2025-08.nc  ${rawpath}/CERES_EBAF_tsrc_2000-03_to_2025-08.nc ${outpath}/temp.nc
cdo chname,tsr,tscre ${outpath}/temp.nc ${outpath}/CERES_EBAF_tscre_2000-03_to_2025-08.nc

# ==================================================================================
