# Detrend (Linear) and Anomalize (remove scycle) using CDO
# Note this doesnt work, I wonder if this is due to the time dimension

expnames=("TCo2559-DART-1950C" "TCo1279-DART-1950") #("TCo319_ctl1950d" "TCo319_ssp585" "TCo2559-DART-1950C" "TCo1279-DART-1950" "TCo1279-DART-2090")
vnames=("Tadv" "ws10" "eis" "w700" "r700") #("sst") #("allsky" "clearsky" "cre" "ttcre" "tscre")
dpath=/home/niu4/gliu8/projects/scrap/regrid_1x1
detrendpath=/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1
scyclepath=/home/niu4/gliu8/projects/scrap/regrid_1x1/scycle
for exp in ${expnames[@]}; do # This loop format is for zsh. Use ${expes[@]} if you are using bash.
    for vname in ${vnames[@]}; do

        infile=${dpath}/${exp}_${vname}_regrid1x1.nc
        scyclefile=${scyclepath}/${exp}_${vname}_regrid1x1.nc
        outfile=${detrendpath}/${exp}_${vname}_regrid1x1.nc

        cdo ymonmean ${infile} ${scyclefile}
        cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
        cdo detrend ${detrendpath}/temp1.nc ${outfile}


        # cdo ymonsub ${dpath}/${exp}_${vname}_regrid1x1.nc -ymonmean ${dpath}/${exp}_${vname}_regrid1x1.nc ${outpath}/temp.nc
        # cdo detrend ${outpath}/temp.nc ${outpath}/${exp}_${vname}_anom_regrid1x1.nc
        # rm ${outpath}/temp.nc
        echo "Completed $vname for $exp"

    done
done




# vnames=("ws10") #"eis" "ws10" "Tadv") #("sst") #("allsky" "clearsky" "cre" "ttcre" "tscre")
# yearstr="1979_2024"
# rawpath='/home/niu4/gliu8/share/ERA5/processed'
# scyclepath='/home/niu4/gliu8/projects/common_data/ERA5/scycle'
# detrendpath='/home/niu4/gliu8/projects/common_data/ERA5/anom_detrend1'

# for vname in ${vnames[@]}; do
  
#     infile=${rawpath}/${vname}_${yearstr}.nc
#     scyclefile=${scyclepath}/${vname}_${yearstr}.nc
#     outfile=${detrendpath}/${vname}_${yearstr}.nc

#     cdo ymonmean ${infile} ${scyclefile}
#     cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
#     cdo detrend ${detrendpath}/temp1.nc ${outfile}

# done

#cdo ymonsub /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo1279-DART-2090_cre_regrid1x1.nc -ymonmean /home/niu4/gliu8/projects/scrap/regrid_1x1/TCo1279-DART-2090_cre_regrid1x1.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/temp.nc
#cdo detrend /home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/temp.nc /home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1/TCo1279-DART-2090_cre_anom_regrid1x1.nc

