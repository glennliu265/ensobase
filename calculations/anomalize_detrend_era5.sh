#! /bin/zsh

# Set the variables
#vnames=("skt" "tsr" "ttr" "tsrc" "ttrc" "u10" "v10" "t700" "w700" "r700")
vnames=("eis" "ws10" "Tadv") #("sst") #("allsky" "clearsky" "cre" "ttcre" "tscre")
yearstr="1979_2024"
rawpath='/home/niu4/gliu8/share/ERA5/processed'
scyclepath='/home/niu4/gliu8/projects/common_data/ERA5/scycle'
detrendpath='/home/niu4/gliu8/projects/common_data/ERA5/anom_detrend1'

for vname in ${vnames[@]}; do
  
    infile=${rawpath}/${vname}_${yearstr}.nc
    scyclefile=${scyclepath}/${vname}_${yearstr}.nc
    outfile=${detrendpath}/${vname}_${yearstr}.nc

    cdo ymonmean ${infile} ${scyclefile}
    cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
    cdo detrend ${detrendpath}/temp1.nc ${outfile}

done