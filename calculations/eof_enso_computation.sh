# Calculate EOF for different cases

inpath='/home/niu4/gliu8/projects/common_data/ERA5/anom_detrend1/sst_1979_2024.nc'
outpath='/home/niu4/gliu8/projects/ccfs/enso_eof'
expname='ERA5_1979_2024'
# Step 1 Crop to tropical pacific (following Takahashi et al. 2011)
cdo sellonlatbox,120,-70,-10,10 $inpath ${outpath}/tp_crop.nc
cdo setmissval,0 ${outpath}/tp_crop.nc ${outpath}/temp.nc
cdo selyear,1979/2024 ${outpath}/temp.nc ${outpath}/tp_crop.nc

# Step 2 Perform EOF (https://code.mpimet.mpg.de/boards/53/topics/5741)
cdo eof,10 ${outpath}/tp_crop.nc ${outpath}/${expname}_eigval.nc ${outpath}/${expname}_pattern.nc
cdo eofcoeff ${outpath}/${expname}_pattern.nc ${outpath}/tp_crop.nc ${outpath}/${expname}_pc.nc
cdo -div ${outpath}/${expname}_eigval.nc -timsum ${outpath}/${expname}_eigval.nc /${outpath}/${expname}_explvar.nc



