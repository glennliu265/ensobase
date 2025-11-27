

# 31km Future
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ssp585_sst.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo319_ssp585_sst_regrid1x1.nc

# 31km Control
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_sst.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo319_ctl1950d_sst_regrid1x1.nc

# 9km Future
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/TCo1279_2090slice_aleph_sst_1m_2090-2099.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo1279-DART-2090_sst_regrid1x1.nc

# Regrid u10, v10, ws10 (9km Control)
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_u10.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo1279-DART-1950_u10_regrid1x1.nc
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_v10.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo1279-DART-1950_v10_regrid1x1.nc
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_ws10.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo1279-DART-1950_ws10_regrid1x1.nc
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_Tadv.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo1279-DART-1950_Tadv_regrid1x1.nc

# Also Anomalize and Detrend

