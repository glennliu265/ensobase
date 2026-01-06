# Take the Global mean of All Sky, Clear Sky, and other Radiative Fluxes
datpath="/home/niu4/gliu8/projects/common_data/ERA5/anom_detrend1"
#outpath="/home/niu4/gliu8/projects/scrap/global_mean"
outpath="/home/niu4/gliu8/projects/ccfs/global_mean"
vnames=("cre" "ttcre" "tscre" "allsky" "clearsky" "sst") #("sst") #
expname="ERA5_1979_2024"
for vname in ${vnames[@]}; do
    
    infile=${datpath}/${vname}_1979_2024.nc
    outfile=${outpath}/${expname}_${vname}_global_mean.nc
    cdo fldmean $infile $outfile
    echo "Completed $vname"
done

# Manual Change for ERA5
datpath="/home/niu4/gliu8/projects/common_data/ERA5/anom_detrend1"
#outpath="/home/niu4/gliu8/projects/scrap/global_mean"
outpath="/home/niu4/gliu8/projects/ccfs/global_mean"
expname="ERA5_1979_2024"
vname="tscre"
infile=${datpath}/${vname}_1979_2024.nc
outfile=${outpath}/${expname}_${vname}_global_mean.nc
cdo fldmean $infile $outfile
echo "Completed $vname"


# =====================
# Repeat for CERES EBAF
# =====================
datpath="/home/niu4/gliu8/share/CERES/processed/anom_detrend1"
outpath="/home/niu4/gliu8/projects/ccfs/global_mean"
vnames=("cre" "ttcre" "tscre" "allsky" "clearsky" "sst") #("sst") #
expname="CERES_EBAF"
for vname in ${vnames[@]}; do
    
    infile=${datpath}/${expname}_${vname}.nc
    outfile=${outpath}/${expname}_${vname}_global_mean.nc
    cdo fldmean $infile $outfile
    echo "Completed $vname"
done

# ==============================
# Repeat for AWI-CM3 Experiments
# ==============================
datpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
outpath="/home/niu4/gliu8/projects/ccfs/global_mean"
vnames=("cre" "ttcre" "tscre" "allsky" "clearsky") #"sst") #("sst") #
expnames=("TCo2559-DART-1950C" "TCo1279-DART-1950" "TCo319_ctl1950d" "TCo319_ssp585" "TCo1279-DART-2090")
for expname in ${expnames[@]}; do
    for vname in ${vnames[@]}; do
        
        infile=${datpath}/${expname}_${vname}.nc
        outfile=${outpath}/${expname}_${vname}_global_mean.nc
        cdo fldmean $infile $outfile
        echo "Completed $vname"
    done
done

# Manually Redo some calculations ~~~~~

vname="cre"
expname="TCo2559-DART-1950C"
infile=${datpath}/${expname}_${vname}.nc
outfile=${outpath}/${expname}_${vname}_global_mean.nc
cdo fldmean $infile $outfile
echo "Completed $vname"

datpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
outpath="/home/niu4/gliu8/projects/ccfs/global_mean"
vname="cre"
expname="TCo319_ssp585"
infile=${datpath}/${expname}_${vname}.nc
outfile=${outpath}/${expname}_${vname}_global_mean.nc
cdo fldmean $infile $outfile
echo "Completed $vname"

datpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
outpath="/home/niu4/gliu8/projects/ccfs/global_mean"
vname="ttcre"
expname="TCo1279-DART-1950"
infile=${datpath}/${expname}_${vname}.nc
outfile=${outpath}/${expname}_${vname}_global_mean.nc
cdo fldmean $infile $outfile
echo "Completed $vname"

