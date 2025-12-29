# Crop File to NAO Region


# 31km Control
cdo sellonlatbox,270,40,20,80 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_control/TCo319_ctl1950d_msl_1m_1850-2134.nc /home/niu4/gliu8/projects/scrap/nao_crop/TCo319_ctl1950d_msl.nc

# 31km Future
cdo sellonlatbox,270,40,20,80 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/TCo319_ssp585_msl_1m_2015-2114.nc /home/niu4/gliu8/projects/scrap/nao_crop/TCo319_ssp585_msl.nc

# 9km 1950
cdo sellonlatbox,270,40,20,80 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/TCo1279_DART-1950_msl_1m_1950-1969.nc /home/niu4/gliu8/projects/scrap/nao_crop/TCo1279-DART-1950_msl.nc

# 9km 2090
cdo sellonlatbox,270,40,20,80 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/TCo1279_DART-2090_msl_1m_2090-2099.nc /home/niu4/gliu8/projects/scrap/nao_crop/TCo1279-DART-2090_msl.nc