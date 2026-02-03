
#  ===================================
# 9km (1950)
#  ===================================

# Rename \10ws --> \ws10
cdo chname,\10ws,ws10 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/TCo1279-DART-1950_atm_remapped_1m_10ws_240months.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_ws10.nc

# Rename \10u --> u10
cdo chname,\10u,u10 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/TCo1279-DART-1950_atm_remapped_1m_10u_240months.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_u10.nc

# Rename \10v --> v10
cdo chname,\10v,v10 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/TCo1279-DART-1950_atm_remapped_1m_10v_240months.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_v10.nc

# Copy over SST (9km)
cp /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/TCo1279-DART-1950_atm_remapped_1m_sst_240months.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_sst.nc


# Copy over Tropospheric Variables
cdo chname,t,t700 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/TCo1279-DART-1950_atm_remapped_1m_t700_240months.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_t700.nc
cdo chname,r,r700 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/TCo1279-DART-1950_atm_remapped_1m_r700_240months.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_r700.nc
cdo chname,w,w700 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/TCo1279-DART-1950_atm_remapped_1m_w700_240months.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_w700.nc



#  ===================================
# 9km (2090)
#  ===================================

# Tropospheric Variables (rename and move to processed folder)
cdo chname,t,t700 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/TCo1279-2090slice_aleph_pl_t_700_1m_2090-2099.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-2090_t700.nc
cdo chname,r,r700 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/TCo1279_2090slice_aleph_pl_r_700_1m_2090-2099.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-2090_r700.nc
cdo chname,w,w700 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/TCo1279-2090slice_aleph_pl_w_700_1m_2090-2099.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-2090_w700.nc

# Wind Variables
cdo chname,\10ws,ws10 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/TCo1279-2090slice_aleph_10ws_1m_2090-2099.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-2090_ws10.nc
cdo chname,\10u,u10 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/TCo1279-2090slice_aleph_10u_1m_2090-2099.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-2090_u10.nc
cdo chname,\10v,v10 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/TCo1279-2090slice_aleph_10v_1m_2090-2099.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-2090_v10.nc

# Copy over sst
cp /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/TCo1279_2090slice_aleph_sst_1m_2090-2099.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-2090_sst.nc

#  ===================================
# 31km (Control) 
#  ===================================

# Tropospheric Variables (rename and move to processed folder)
cdo chname,t,t700 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ctl1950d/monthly/awicm3_tco319_ctl1950d_pl_t_700_1950-2100.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_t700.nc
cdo chname,r,r700 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ctl1950d/monthly/awicm3_tco319_ctl1950d_pl_r_700_1950-2100.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_r700.nc
cdo chname,w,w700 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ctl1950d/monthly/awicm3_tco319_ctl1950d_pl_w_700_1950-2100.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_w700.nc

# 10u/10/10ws: Rename files merged the aleph_data_pull_consolidate.sh
cdo chname,10u,u10 /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_10u.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_u10.nc
rm /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_10u.nc

cdo chname,10v,v10 /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_10v.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_v10.nc
rm /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_10v.nc

cdo chname,10ws,ws10 /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_10ws.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_ws10.nc
rm /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_10ws.nc


