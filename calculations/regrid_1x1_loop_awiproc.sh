# Regrid AWI files in a loop
expnames=("TCo319_ctl1950d") # "TCo319_ssp585" "TCo2559-DART-1950C" "TCo1279-DART-1950" "TCo1279-DART-2090")
vnames=("sst") #"ttcre" "tscre")
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
