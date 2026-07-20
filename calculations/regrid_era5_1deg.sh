
# # Regrid for anomalized and detrended
# vnames=("creln") #"sst" "eis" "ws10" "Tadv" "w700" "r700" "cre" "mask") 
# rawpath='/home/niu4/gliu8/projects/common_data/ERA5/anom_detrend1'
# outpath='/home/niu4/gliu8/projects/common_data/ERA5/regrid_1x1/anom_detrend1'

# Regrid for raw data
vnames=("skt") #"u10" "v10") #("allsky" "clearsky" "cre" "ttcre" "tscre" "ttr" "ttrc" "tsr" "tsrc") #"sst" "eis" "ws10" "Tadv" "w700" "r700" "cre" "mask") 
rawpath='/home/niu4/gliu8/share/ERA5/processed'
outpath='/home/niu4/gliu8/projects/common_data/ERA5/regrid_1x1'
yearstr="1979_2024"

for vname in ${vnames[@]}; do
    
    infile=${rawpath}/${vname}_${yearstr}.nc
    outfile=${outpath}/${vname}_${yearstr}.nc
    
    cdo selyear,1979/2024 ${infile} ${outpath}/temp.nc
    cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc $infile $outfile

done


# (2026.07.20) Special Regrid of ERA5 Land Mask
infile=/home/niu4/gliu8/projects/common_data/ERA5/era5_landmask_fromsst.nc
outfile=/home/niu4/gliu8/projects/common_data/ERA5/era5_landmask_fromsst_regrid1x1.nc
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc $infile $outfile