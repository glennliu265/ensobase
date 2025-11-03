
#! /bin/zsh
# Works with regridded datasets, add things together
# Indicate Experiment Name
#expname="TCo1279-DART-1950"
##

expnames=("TCo319_ctl1950d" "TCo319_ssp585" "TCo2559-DART-1950C" "TCo1279-DART-1950" "TCo1279-DART-2090")
dpath=/home/niu4/gliu8/projects/scrap/regrid_1x1

for exp in $expnames; do # This loop format is for zsh. Use ${expes[@]} if you are using bash.

    # Calcualte Allsky
    cdo -add ${dpath}/${exp}_tsr_regrid1x1.nc ${dpath}/${exp}_ttr_regrid1x1.nc temp.nc
    cdo chname,tsr,allsky temp.nc ${dpath}/${exp}_allsky_regrid1x1.nc 

    # Calculate Clearsky
    cdo -add ${dpath}/${exp}_tsrc_regrid1x1.nc ${dpath}/${exp}_ttrc_regrid1x1.nc temp.nc
    cdo chname,tsrc,clearsky temp.nc ${dpath}/${exp}_clearsky_regrid1x1.nc

    # Calculate tsr cre
    cdo -sub ${dpath}/${exp}_tsr_regrid1x1.nc ${dpath}/${exp}_tsrc_regrid1x1.nc temp.nc #
    cdo chname,tsr,tscre temp.nc ${dpath}/${exp}_tscre_regrid1x1.nc


    # Calculate ttr cre
    cdo -sub ${dpath}/${exp}_ttr_regrid1x1.nc ${dpath}/${exp}_ttrc_regrid1x1.nc temp.nc #
    cdo chname,ttr,ttcre temp.nc ${dpath}/${exp}_ttcre_regrid1x1.nc

    # Calculare cre
    cdo -sub ${dpath}/${exp}_allsky_regrid1x1.nc ${dpath}/${exp}_clearsky_regrid1x1.nc temp.nc
    cdo chname,allsky,cre temp.nc ${dpath}/${exp}_cre_regrid1x1.nc

    echo "Completed $exp"
done

# # Manual Loop
# # Calculate allsky
# #cdo enssum TCo319_ctl1950d_tsr_regrid1x1.nc TCo319_ctl1950d_ttr_regrid1x1.nc TCo319_ctl1950d_allsky_regrid1x1.nc
# cdo -add TCo319_ctl1950d_tsr_regrid1x1.nc TCo319_ctl1950d_ttr_regrid1x1.nc temp.nc
# cdo chname,tsr,allsky temp.nc TCo319_ctl1950d_allsky_regrid1x1.nc

# # Calculate clearsky
# #cdo enssum TCo319_ctl1950d_tsrc_regrid1x1.nc TCo319_ctl1950d_ttrc_regrid1x1.nc TCo319_ctl1950d_clearsky_regrid1x1.nc
# cdo -add TCo319_ctl1950d_tsrc_regrid1x1.nc TCo319_ctl1950d_ttrc_regrid1x1.nc temp.nc
# cdo chname,tsrc,clearsky temp.nc TCo319_ctl1950d_clearsky_regrid1x1.nc

# # Calculate tsr cre
# cdo -sub TCo319_ctl1950d_tsr_regrid1x1.nc TCo319_ctl1950d_tsrc_regrid1x1.nc temp.nc #
# cdo chname,tsr,tscre temp.nc TCo319_ctl1950d_tscre_regrid1x1.nc

# # calculate ttr cre
# cdo -sub TCo319_ctl1950d_ttr_regrid1x1.nc TCo319_ctl1950d_ttrc_regrid1x1.nc temp.nc #
# cdo chname,ttr,ttcre temp.nc TCo319_ctl1950d_ttcre_regrid1x1.nc

# # Calculare cre
# cdo -sub TCo319_ctl1950d_allsky_regrid1x1.nc TCo319_ctl1950d_clearsky_regrid1x1.nc temp.nc
# cdo chname,allsky,cre temp.nc TCo319_ctl1950d_cre_regrid1x1.nc
# #cdo enssum TCo319_ctl1950d_tsrc_regrid1x1.nc TCo319_ctl1950d_ttrc_regrid1x1.nc TCo319_ctl1950d_clearsky_regrid1x1.nc