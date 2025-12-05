

# FBCT -------

# Make Folders
mkdir /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_FBCT_ERA5
mkdir /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_FBCT_ERA5/raw/
mkdir /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_FBCT_ERA5/anom_detrend1/

# Copy Over Radiation
cp /home/niu4/gliu8/share/CERES/processed/CERES_FBCT_creln_anom_2002-07_to_2023-02.nc /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_FBCT_ERA5/anom_detrend1/CERES_FBCT_creln.nc

# Copy Over CCFs
vnames=("sst" "eis" "Tadv" "ws10" "w700" "r700")
for vname in ${vnames[@]}; do
    cp /home/niu4/gliu8/projects/common_data/ERA5/regrid_1x1/anom_detrend1/${vname}_1979_2024.nc /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_FBCT_ERA5/anom_detrend1/${vname}.nc
done

# EBAF -------

# Make Folders
mkdir /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_EBAF_ERA5
mkdir /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_EBAF_ERA5/raw/
mkdir /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_EBAF_ERA5/anom_detrend1/

# Anomalize and detrend cre
cdo ymonmean /home/niu4/gliu8/share/CERES/processed/CERES_EBAF_cre_2000-03_to_2025-08.nc  /home/niu4/gliu8/share/CERES/processed/CERES_EBAF_cre_2000-03_to_2025-08_scycle.nc
cdo ymonsub /home/niu4/gliu8/share/CERES/processed/CERES_EBAF_cre_2000-03_to_2025-08.nc /home/niu4/gliu8/share/CERES/processed/CERES_EBAF_cre_2000-03_to_2025-08_scycle.nc /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_EBAF_ERA5/anom_detrend1/temp1.nc
cdo detrend /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_EBAF_ERA5/anom_detrend1/temp1.nc /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_EBAF_ERA5/anom_detrend1/cre.nc

# Copy ERA5 data over
vnames=("sst" "eis" "Tadv" "ws10" "w700" "r700")
for vname in ${vnames[@]}; do
    cp /home/niu4/gliu8/projects/common_data/ERA5/regrid_1x1/anom_detrend1/${vname}_1979_2024.nc /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_EBAF_ERA5/anom_detrend1/${vname}.nc
done

# Calculate EBAF Unobscured Component
cdo -expr,"nucc=1-ucc" /home/niu4/gliu8/share/CERES/processed/CERES_FBCT_ucc_2002-07_to_2023-02.nc /home/niu4/gliu8/share/CERES/processed/CERES_FBCT_nucc_2002-07_to_2023-02.nc
cdo -div, /home/niu4/gliu8/share/CERES/processed/CERES_FBCT_lcc_2002-07_to_2023-02.nc /home/niu4/gliu8/share/CERES/processed/CERES_FBCT_nucc_2002-07_to_2023-02.nc temp1.nc
cdo chname,lcc,lncc temp1.nc /home/niu4/gliu8/share/CERES/processed/CERES_FBCT_lncc_2002-07_to_2023-02.nc

# Anomalize and Detrend EBAF 
datpath="/home/niu4/gliu8/share/CERES/processed"
outpath="/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_EBAF_ERA5/anom_detrend1"
cdo ymonmean ${datpath}/CERES_EBAF_creln_2000-03_to_2025-08.nc ${datpath}/CERES_EBAF_creln_2000-03_to_2025-08_scycle.nc
cdo ymonsub ${datpath}/CERES_EBAF_creln_2000-03_to_2025-08.nc ${datpath}/CERES_EBAF_creln_2000-03_to_2025-08_scycle.nc ${datpath}/temp1.nc
cdo detrend ${datpath}/temp1.nc ${outpath}/creln.nc


# ERA5_1979_2024 ------- 

# Make Folders
expname="ERA5_1979_2024"
mkdir /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/${expname}/
mkdir /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/${expname}/raw/
mkdir /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/${expname}/anom_detrend1/

# Copy ERA5 Data Over (regridded)
expname="ERA5_1979_2024"
vnames=("sst" "eis" "Tadv" "ws10" "w700" "r700" "cre" "creln")
for vname in ${vnames[@]}; do
    cp /home/niu4/gliu8/projects/common_data/ERA5/regrid_1x1/anom_detrend1/${vname}_1979_2024.nc /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/${expname}/anom_detrend1/${vname}.nc
done

# ========================================
# Clone Data for a EBAF Range
# ========================================
expname="ERA5_2000_2024"
mkdir /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/${expname}/
mkdir /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/${expname}/raw/
mkdir /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/${expname}/anom_detrend1/
vnames=("sst" "eis" "Tadv" "ws10" "w700" "r700" "cre")

# Subset, Anom/Detrend, Regrid
expname="ERA5_2000_2024"
vnames=("sst" "eis" "Tadv" "ws10" "w700" "r700" "cre")
for vname in ${vnames[@]}; do
    cdo selyear,2000/2024 /home/niu4/gliu8/share/ERA5/processed/${vname}_1979_2024.nc /home/niu4/gliu8/share/ERA5/processed/rawtemp.nc
    cdo ymonmean /home/niu4/gliu8/share/ERA5/processed/rawtemp.nc /home/niu4/gliu8/share/ERA5/processed/scycletemp.nc
    cdo ymonsub /home/niu4/gliu8/share/ERA5/processed/rawtemp.nc /home/niu4/gliu8/share/ERA5/processed/scycletemp.nc /home/niu4/gliu8/share/ERA5/processed/anomtemp.nc
    cdo detrend /home/niu4/gliu8/share/ERA5/processed/anomtemp.nc /home/niu4/gliu8/projects/ccfs/input_data/ERA5_2000_2024/anom_detrend1/${vname}.nc
    cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /home/niu4/gliu8/projects/ccfs/input_data/ERA5_2000_2024/anom_detrend1/${vname}.nc /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/${expname}/anom_detrend1/${vname}.nc
done


# Crop/Process msl in ERA5
cdo selyear,1979/2024 /home/niu4/gliu8/share/ERA5/downloaded/msl_1940_2024.nc /home/niu4/gliu8/share/ERA5/processed/msl_1979_2024.nc

# Scrap ----

# Calculate tcc in ERA5
cdo -add /home/niu4/gliu8/share/ERA5/processed/ucc_1979_2024.nc  /home/niu4/gliu8/share/ERA5/processed/lcc_1979_2024.nc /home/niu4/gliu8/share/ERA5/processed/temp.nc
cdo chname,ucc,tcc /home/niu4/gliu8/share/ERA5/processed/temp.nc /home/niu4/gliu8/share/ERA5/processed/tcc_1979_2024.nc

# Regrid ERA5 tcc
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /home/niu4/gliu8/share/ERA5/processed/tcc_1979_2024.nc /home/niu4/gliu8/projects/common_data/ERA5/regrid_1x1/tcc_1979_2024.nc

# ========================================
# Setup Input Data for AWI-CM3 simulations
# ========================================

# 9km 1950
expname="TCo1279-DART-1950"
mkdir /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/${expname}/
mkdir /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/${expname}/raw/
mkdir /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/${expname}/anom_detrend1/
vnames=("sst" "eis" "Tadv" "ws10" "w700" "r700")
for vname in ${vnames[@]}; do
    cp /home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/${expname}_${vname}_regrid1x1.nc /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/${expname}/anom_detrend1/${vname}.nc
done
cp /home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/TCo1279-DART-1950_cre_regrid1x1.nc /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/TCo1279-DART-1950/anom_detrend1/cre.nc


# ========================================================
# Anomalize and Detrend CERES for global mean calculations
# ========================================================
# See "global_mean_ERA5.sh" for next calculations
vnames=("cre" "allsky" "clearsky" "tscre" "ttcre")
rawpath='/home/niu4/gliu8/share/CERES/processed'
scyclepath='/home/niu4/gliu8/share/CERES/processed/scycle'
detrendpath='/home/niu4/gliu8/share/CERES/processed/anom_detrend1'
expname="CERES_EBAF"
for vname in ${vnames[@]}; do
    
    infile=${rawpath}/${expname}_${vname}_2000-03_to_2025-08.nc
    scyclefile=${scyclepath}/${expname}_${vname}.nc
    outfile=${detrendpath}/${expname}_${vname}.nc

    cdo ymonmean ${infile} ${scyclefile}
    cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
    cdo detrend ${detrendpath}/temp1.nc ${outfile}

done