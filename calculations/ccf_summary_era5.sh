
vnames=("sst" "eis" "ws10" "Tadv" "w700" "r700") 
rawpath='/home/niu4/gliu8/share/ERA5/processed'
outpath='/home/niu4/gliu8/projects/ccfs/time_mean/'
yearstr="1979_2024"

for vname in ${vnames[@]}; do
  
    infile=${rawpath}/${vname}_${yearstr}.nc
    outfile=${outpath}/${vname}_${yearstr}.nc
    
    cdo selyear,1979/2024 ${infile} ${outpath}/temp.nc
    cdo timmean ${outpath}/temp.nc ${outfile}
done
