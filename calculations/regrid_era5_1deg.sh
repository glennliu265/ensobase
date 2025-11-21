
vnames=("creln") #"sst" "eis" "ws10" "Tadv" "w700" "r700" "cre" "mask") 
rawpath='/home/niu4/gliu8/projects/common_data/ERA5/anom_detrend1'
outpath='/home/niu4/gliu8/projects/common_data/ERA5/regrid_1x1/anom_detrend1'
yearstr="1979_2024"

for vname in ${vnames[@]}; do
    
    infile=${rawpath}/${vname}_${yearstr}.nc
    outfile=${outpath}/${vname}_${yearstr}.nc
    
    cdo selyear,1979/2024 ${infile} ${outpath}/temp.nc
    cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc $infile $outfile

done