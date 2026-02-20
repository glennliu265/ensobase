# Regrid AWI files in a loop
expnames=("TCo319_ctl1950d") # "TCo319_ssp585" "TCo2559-DART-1950C" "TCo1279-DART-1950" "TCo1279-DART-2090")
vnames=("cre") #"ttcre" "tscre")
dpath=/home/niu4/gliu8/projects/scrap/regrid_1x1
for exp in ${expnames[@]}; do # This loop format is for zsh. Use ${expes[@]} if you are using bash.
    for vname in ${vnames[@]}; do 
        cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc /home/niu4/gliu8/projects/scrap/processed_global/${exp}_${vname}.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/${exp}_${vname}_regrid1x1.nc
        echo "Completed $vname for $exp"
    done
done

# Same as above, but regrid anomalized output (2026.01.06 Update)
expnames=("TCo1279-DART-2060" "TCo319-DART-ctl1950d-gibbs-charn" "TCo319-DART-ssp585d-gibbs-charn")
vnames=("msl") #"ttcre" "tscre")
rawpath=/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1
dpath=/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1
for exp in ${expnames[@]}; do # This loop format is for zsh. Use ${expes[@]} if you are using bash.
    for vname in ${vnames[@]}; do 
        cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc ${rawpath}/${exp}_${vname}.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/${exp}_${vname}_regrid1x1.nc
        echo "Completed $vname for $exp"
    done
done

# Same as above, but regrid anomalized output (2026.01.06 Update)
expnames=("TCo1279-DART-2060" "TCo319-DART-ctl1950d-gibbs-charn" "TCo319-DART-ssp585d-gibbs-charn")
vnames=("msl") #"ttcre" "tscre")
rawpath=/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1
dpath=/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1
for exp in ${expnames[@]}; do # This loop format is for zsh. Use ${expes[@]} if you are using bash.
    for vname in ${vnames[@]}; do 
        cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc ${rawpath}/${exp}_${vname}.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/${exp}_${vname}_regrid1x1.nc
        echo "Completed $vname for $exp"
    done
done

# Also Regrid FULL ERA5 variables
# And Crop between 2001 and 2024, and copy into folder
expnames=("ERA5")
vnames=("ttcre" "tscre") #("sst" "eis" "Tadv" "w700" "r700" "ws10")
rawpath=/home/niu4/gliu8/share/ERA5/processed
dpath=/home/niu4/gliu8/projects/common_data/ERA5/regrid_1x1
dpath2=/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_EBAF_ERA5_2001_2024/raw
for exp in ${expnames[@]}; do # This loop format is for zsh. Use ${expes[@]} if you are using bash.
    for vname in ${vnames[@]}; do 
        cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc ${rawpath}/${vname}_1979_2024.nc ${dpath}/${vname}_1979_2024.nc
        cdo selyear,2001/2024 ${dpath}/${vname}_1979_2024.nc ${dpath2}/${vname}.nc
        echo "Completed $vname for $exp"
    done
done

# Redo for EIS
rawpath=/home/niu4/gliu8/share/ERA5/processed
dpath=/home/niu4/gliu8/projects/common_data/ERA5/regrid_1x1
dpath2=/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_EBAF_ERA5_2001_2024/raw
vname="eis"
cdo remapbil,/home/niu4/gliu8/projects/common_data/regrid_re1x1.nc ${rawpath}/${vname}_1979_2024.nc ${dpath}/${vname}_1979_2024.nc
cdo selyear,2001/2024 ${dpath}/${vname}_1979_2024.nc ${dpath2}/${vname}.nc

# Copy over cre from other experiment
cdo selyear,2001/2024 /home/niu4/gliu8/share/CERES/processed/CERES_EBAF_cre_2000-03_to_2025-08.nc /home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/CERES_EBAF_ERA5_2001_2024/raw/cre.nc