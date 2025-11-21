
# Do LnCC calculation for 9km simulation
# Calculate 1-U
cdo -z zip_7 -expr,"nucc=1-ucc" TCo1279-DART-1950_ucc.nc TCo1279-DART-1950_nucc.nc
# Calculate L/(1-U)
cdo -div, /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/TCo1279-DART-1950_atm_remapped_1m_lcc_240months.nc TCo1279-DART-1950_nucc.nc temp.nc
cdo chname,lcc,lncc temp.nc TCo1279-DART-1950_lncc.nc


# Repeat for 5km imulation
# Calculate 1-U
cdo -expr,"nucc=1-ucc" TCo2559-DART-1950C_ucc.nc TCo2559-DART-1950C_nucc.nc
# Calculate L/(1-U)
cdo -div, /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo2559-DART-1950C_5km/atm/TCo2559-DART-1950C_atm_10256x5120_1m_lcc_1950-1959.nc TCo2559-DART-1950C_nucc.nc temp1.nc
cdo -z zip_7 chname,lcc,lncc temp1.nc TCo2559-DART-1950C_lncc.nc

# Do calculation for ERA5 (see Obsidian md cdo_scrap_20251114)
# calculate 1-U
cdo -expr,"nucc=1-ucc" /home/niu4/gliu8/share/ERA5/processed/ucc_1979_2024.nc /home/niu4/gliu8/share/ERA5/processed/nucc_1979_2024.nc
# Calculate L/(1-U)
cdo -div, /home/niu4/gliu8/share/ERA5/processed/lcc_1979_2024.nc /home/niu4/gliu8/share/ERA5/processed/nucc_1979_2024.nc temp1.nc
cdo chname,lcc,lncc temp1.nc lncc_1979_2024.nc



