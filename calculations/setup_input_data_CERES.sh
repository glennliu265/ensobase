
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





