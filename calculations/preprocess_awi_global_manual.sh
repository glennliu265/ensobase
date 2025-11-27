
# 9km (1950)

# Rename \10ws --> \ws10
cdo chname,\10ws,ws10 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/TCo1279-DART-1950_atm_remapped_1m_10ws_240months.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_ws10.nc

# Rename \10u --> u10
cdo chname,\10u,u10 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/TCo1279-DART-1950_atm_remapped_1m_10u_240months.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_u10.nc

# Rename \10v --> v10
cdo chname,\10v,v10 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/TCo1279-DART-1950_atm_remapped_1m_10v_240months.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_v10.nc

# Copy over SST (9km)
cp /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/TCo1279-DART-1950_atm_remapped_1m_sst_240months.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_sst.nc

# 9km (2090)


