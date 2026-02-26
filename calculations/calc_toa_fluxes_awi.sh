
#! /bin/zsh
# Works with regridded datasets, add things together
# Indicate Experiment Name
#expname="TCo1279-DART-1950"
##

expnames=("TCo319_ctl1950d" "TCo319_ssp585" "TCo2559-DART-1950C" "TCo1279-DART-1950" "TCo1279-DART-2090")
dpath=/home/niu4/gliu8/projects/scrap/regrid_1x1

for exp in ${expnames[@]}; do # This loop format is for zsh. Use ${expes[@]} if you are using bash.

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

# ===================================================
# Update 2026.02.26: Work with non re-gridded outputs
# ===================================================
expnames=("TCo319-DART-ctl1950d-gibbs-charn" "TCo319-DART-ssp585d-gibbs-charn") # "TCo1279-DART-2060")
dpath=/home/niu4/gliu8/projects/scrap/processed_global
for exp in ${expnames[@]}; do # This loop format is for zsh. Use ${expes[@]} if you are using bash.

    # Calcualte Allsky
    cdo -L -chname,tsr,allsky -add ${dpath}/${exp}_tsr.nc ${dpath}/${exp}_ttr.nc ${dpath}/${exp}_allsky.nc 

    # Calculate Clearsky
    cdo -L -chname,tsrc,clearsky -add ${dpath}/${exp}_tsrc.nc ${dpath}/${exp}_ttrc.nc ${dpath}/${exp}_clearsky.nc

    # Calculate tsr cre
    cdo -L -chname,tsr,tscre -sub ${dpath}/${exp}_tsr.nc ${dpath}/${exp}_tsrc.nc ${dpath}/${exp}_tscre.nc
    
    # Calculate ttr cre
    cdo -L -chname,ttr,ttcre -sub ${dpath}/${exp}_ttr.nc ${dpath}/${exp}_ttrc.nc ${dpath}/${exp}_ttcre.nc

    # Calculare cre
    cdo -L -chname,allsky,cre -sub ${dpath}/${exp}_allsky.nc ${dpath}/${exp}_clearsky.nc ${dpath}/${exp}_cre.nc

    echo "Completed $exp"
done
