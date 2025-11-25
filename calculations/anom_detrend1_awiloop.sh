# Detrend (Linear) and Anomalize (remove scycle) using CDO
# Note this doesnt work, I wonder if this is due to the time dimension

expnames=("TCo319_ctl1950d" "TCo319_ssp585" "TCo2559-DART-1950C" "TCo1279-DART-1950" "TCo1279-DART-2090")
vnames=("sst") #("allsky" "clearsky" "cre" "ttcre" "tscre")
dpath=/home/niu4/gliu8/projects/scrap/regrid_1x1
outpath=/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1
for exp in ${expnames[@]}; do # This loop format is for zsh. Use ${expes[@]} if you are using bash.
    for vname in ${vnames[@]}; do 
        cdo ymonsub ${dpath}/${exp}_${vname}_regrid1x1.nc -ymonmean ${dpath}/${exp}_${vname}_regrid1x1.nc ${outpath}/temp.nc
        cdo detrend ${outpath}/temp.nc ${outpath}/${exp}_${vname}_anom_regrid1x1.nc
        rm ${outpath}/temp.nc
        echo "Completed $vname for $exp"
    done
done



#cdo ymonsub /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo1279-DART-2090_cre_regrid1x1.nc -ymonmean /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo1279-DART-2090_cre_regrid1x1.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/temp.nc
#cdo detrend /home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/temp.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/TCo1279-DART-2090_cre_anom_regrid1x1.nc

