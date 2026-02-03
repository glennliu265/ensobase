# Detrend (Linear) and Anomalize (remove scycle) using CDO
# Note this doesnt work, I wonder if this is due to the time dimension

# expnames=("TCo2559-DART-1950C" "TCo1279-DART-1950" "TCo319_ctl1950d" "TCo319_ssp585" "TCo1279-DART-2090")
# vnames=("cre" "allsky" "clearsky") #("Tadv" "ws10" "eis" "w700" "r700" "sst") #("allsky" "clearsky" "cre" "ttcre" "tscre")
# dpath=/home/niu4/gliu8/projects/scrap/regrid_1x1
# detrendpath=/home/niu4/gliu8/projects/scrap/regrid_1x1/global_anom_detrend1
# scyclepath=/home/niu4/gliu8/projects/scrap/regrid_1x1/scycle
# for exp in ${expnames[@]}; do # This loop format is for zsh. Use ${expes[@]} if you are using bash.
#     for vname in ${vnames[@]}; do

#         infile=${dpath}/${exp}_${vname}_regrid1x1.nc
#         scyclefile=${scyclepath}/${exp}_${vname}_regrid1x1.nc
#         outfile=${detrendpath}/${exp}_${vname}_regrid1x1.nc

#         cdo ymonmean ${infile} ${scyclefile}
#         cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
#         cdo detrend ${detrendpath}/temp1.nc ${outfile}


#         # cdo ymonsub ${dpath}/${exp}_${vname}_regrid1x1.nc -ymonmean ${dpath}/${exp}_${vname}_regrid1x1.nc ${outpath}/temp.nc
#         # cdo detrend ${outpath}/temp.nc ${outpath}/${exp}_${vname}_anom_regrid1x1.nc
#         # rm ${outpath}/temp.nc
#         echo "Completed $vname for $exp"

#     done
#done

# Repeat for nonregridded case
expnames=("TCo2559-DART-1950C" "TCo1279-DART-1950" "TCo319_ctl1950d" "TCo319_ssp585" "TCo1279-DART-2090")
vnames=("cre" "allsky" "clearsky" "sst" "ttcre" "tscre") #("Tadv" "ws10" "eis" "w700" "r700" "sst") #("allsky" "clearsky" "cre" "ttcre" "tscre")
dpath="/home/niu4/gliu8/projects/scrap/processed_global"
detrendpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
scyclepath="/home/niu4/gliu8/projects/scrap/processed_global/scycle"
for exp in ${expnames[@]}; do # This loop format is for zsh. Use ${expes[@]} if you are using bash.
    for vname in ${vnames[@]}; do

        infile=${dpath}/${exp}_${vname}.nc
        scyclefile=${scyclepath}/${exp}_${vname}.nc
        outfile=${detrendpath}/${exp}_${vname}.nc

        cdo ymonmean ${infile} ${scyclefile}
        cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
        cdo detrend ${detrendpath}/temp1.nc ${outfile}



        echo "Completed $vname for $exp"

    done
done

# 2026.01.06: Run for Updated AWI-CM3 @ Original Resolution
expnames=("TCo1279-DART-2060" "TCo319-DART-ctl1950d-gibbs-charn" "TCo319-DART-ssp585d-gibbs-charn")
vnames=("sst" "msl" "lcc")
dpath="/home/niu4/gliu8/projects/scrap/processed_global"
detrendpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
scyclepath="/home/niu4/gliu8/projects/scrap/processed_global/scycle"
for exp in ${expnames[@]}; do # This loop format is for zsh. Use ${expes[@]} if you are using bash.
    for vname in ${vnames[@]}; do

        infile=${dpath}/${exp}_${vname}.nc
        scyclefile=${scyclepath}/${exp}_${vname}.nc
        outfile=${detrendpath}/${exp}_${vname}.nc

        cdo ymonmean ${infile} ${scyclefile}
        cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
        cdo detrend ${detrendpath}/temp1.nc ${outfile}
        echo "Completed $vname for $exp"

    done
done


# 2026.02.02: Run for AWI-CM3, Regridded CCF Raw Output
expnames=("TCo319_ctl1950d")
vnames=("sst" "cre" "w700" "r700" "Tadv" "eis" "ws10")

for exp in ${expnames[@]}; do # This loop format is for zsh. Use ${expes[@]} if you are using bash.
    dpath="/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/${exp}/raw"
    detrendpath="/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/${exp}/anom_detrend1"
    scyclepath="/home/niu4/gliu8/projects/ccfs/input_data/regrid_1x1/${exp}/raw/scycle"

    for vname in ${vnames[@]}; do

        infile=${dpath}/${vname}.nc
        scyclefile=${scyclepath}/${vname}.nc
        outfile=${detrendpath}/${vname}.nc

        cdo ymonmean ${infile} ${scyclefile}
        cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
        cdo detrend ${detrendpath}/temp1.nc ${outfile}
        echo "Completed $vname for $exp"

    done
done



# Manually Do this for 2 missing cases

# exp="TCo1279-DART-1950"
# vname="ttcre"
# dpath="/home/niu4/gliu8/projects/scrap/processed_global"
# detrendpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
# scyclepath="/home/niu4/gliu8/projects/scrap/processed_global/scycle"
# infile=${dpath}/${exp}_${vname}.nc
# scyclefile=${scyclepath}/${exp}_${vname}.nc
# outfile=${detrendpath}/${exp}_${vname}.nc
# cdo ymonmean ${infile} ${scyclefile}
# cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
# cdo detrend ${detrendpath}/temp1.nc ${outfile}

exp="TCo319_ssp585"
vname="cre"
dpath="/home/niu4/gliu8/projects/scrap/processed_global"
detrendpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
scyclepath="/home/niu4/gliu8/projects/scrap/processed_global/scycle"
infile=${dpath}/${exp}_${vname}.nc
scyclefile=${scyclepath}/${exp}_${vname}.nc
outfile=${detrendpath}/${exp}_${vname}.nc
cdo ymonmean ${infile} ${scyclefile}
cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
cdo detrend ${detrendpath}/temp1.nc ${outfile}

exp="TCo1279-DART-2090"
vname="cre"
dpath="/home/niu4/gliu8/projects/scrap/processed_global"
detrendpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
scyclepath="/home/niu4/gliu8/projects/scrap/processed_global/scycle"
infile=${dpath}/${exp}_${vname}.nc
scyclefile=${scyclepath}/${exp}_${vname}.nc
outfile=${detrendpath}/${exp}_${vname}.nc
cdo ymonmean ${infile} ${scyclefile}
cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
cdo detrend ${detrendpath}/temp1.nc ${outfile}


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


# ======================================================
# Manually Detrend some variables for hackathon analysis


## MSL -------------------------------------------------
# MSL 9km 1950
exp="TCo1279-DART-1950"
vname="msl"
ncname="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/TCo1279_DART-1950_msl_1m_1950-1969.nc"
detrendpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
scyclepath="/home/niu4/gliu8/projects/scrap/processed_global/scycle"
infile=${ncname}
scyclefile=${scyclepath}/${exp}_${vname}.nc
outfile=${detrendpath}/${exp}_${vname}.nc
cdo ymonmean ${infile} ${scyclefile}
cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
cdo detrend ${detrendpath}/temp1.nc ${outfile}

# MSL 9km 2090
exp="TCo1279-DART-2090"
vname="msl"
ncname="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/TCo1279_DART-2090_msl_1m_2090-2099.nc"
detrendpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
scyclepath="/home/niu4/gliu8/projects/scrap/processed_global/scycle"
infile=${ncname}
scyclefile=${scyclepath}/${exp}_${vname}.nc
outfile=${detrendpath}/${exp}_${vname}.nc
cdo ymonmean ${infile} ${scyclefile}
cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
cdo detrend ${detrendpath}/temp1.nc ${outfile}

# MSL 31km ssp585
exp="TCo319_ssp585"
vname="msl"
ncname="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/TCo319_ssp585_msl_1m_2015-2114.nc"
detrendpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
scyclepath="/home/niu4/gliu8/projects/scrap/processed_global/scycle"
infile=${ncname}
scyclefile=${scyclepath}/${exp}_${vname}.nc
outfile=${detrendpath}/${exp}_${vname}.nc
cdo ymonmean ${infile} ${scyclefile}
cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
cdo detrend ${detrendpath}/temp1.nc ${outfile}

# MSL 31km control (WARNING SeLECT YEAR Implemented)
exp="TCo319_ctl1950d"
vname="msl"
ncname="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_control/TCo319_ctl1950d_msl_1m_1850-2134.nc"
detrendpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
scyclepath="/home/niu4/gliu8/projects/scrap/processed_global/scycle"
infile=${ncname}
scyclefile=${scyclepath}/${exp}_${vname}.nc
outfile=${detrendpath}/${exp}_${vname}.nc
cdo selyear,1950/2100 ${ncname} ${detrendpath}/temp.nc
cdo ymonmean ${detrendpath}/temp.nc ${scyclefile}
cdo ymonsub ${detrendpath}/temp.nc ${scyclefile} ${detrendpath}/temp1.nc
cdo detrend ${detrendpath}/temp1.nc ${outfile}


# Remap pr -----------------------
# 9km 1950
exp="TCo1279-DART-1950"
vname="pr"
ncname="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/TCo1279-DART-1950_atm_remapped_1m_pr_240months.nc"
detrendpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
scyclepath="/home/niu4/gliu8/projects/scrap/processed_global/scycle"
infile=${ncname}
scyclefile=${scyclepath}/${exp}_${vname}.nc
outfile=${detrendpath}/${exp}_${vname}.nc
cdo ymonmean ${infile} ${scyclefile}
cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
cdo detrend ${detrendpath}/temp1.nc ${outfile}


# 9km 2090
exp="TCo1279-DART-2090"
vname="pr"
ncname="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/TCo1279-2090slice_aleph_pr_1m_2090-2099.nc"
detrendpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
scyclepath="/home/niu4/gliu8/projects/scrap/processed_global/scycle"
infile=${ncname}
scyclefile=${scyclepath}/${exp}_${vname}.nc
outfile=${detrendpath}/${exp}_${vname}.nc
cdo ymonmean ${infile} ${scyclefile}
cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
cdo detrend ${detrendpath}/temp1.nc ${outfile}


# lcc -------------------------------------------------
exp="TCo1279-DART-2090"
vname="lcc"
ncname="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/TCo1279_2090slice_aleph_lcc_1m_2090-2099.nc"
detrendpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
scyclepath="/home/niu4/gliu8/projects/scrap/processed_global/scycle"
infile=${ncname}
scyclefile=${scyclepath}/${exp}_${vname}.nc
outfile=${detrendpath}/${exp}_${vname}.nc
cdo ymonmean ${infile} ${scyclefile}
cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
cdo detrend ${detrendpath}/temp1.nc ${outfile}

exp="TCo1279-DART-1950"
vname="lcc"
ncname="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/atm/mon/TCo1279-DART-1950_atm_remapped_1m_lcc_240months.nc"
detrendpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
scyclepath="/home/niu4/gliu8/projects/scrap/processed_global/scycle"
infile=${ncname}
scyclefile=${scyclepath}/${exp}_${vname}.nc
outfile=${detrendpath}/${exp}_${vname}.nc
cdo ymonmean ${infile} ${scyclefile}
cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
cdo detrend ${detrendpath}/temp1.nc ${outfile}

# 31km Control
exp="TCo319_ctl1950d"
vname="lcc"
ncname="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_control/ctl1950d_atm_remapped_1m_lcc_1850-2134.nc"
detrendpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
scyclepath="/home/niu4/gliu8/projects/scrap/processed_global/scycle"
infile=${ncname}
scyclefile=${scyclepath}/${exp}_${vname}.nc
outfile=${detrendpath}/${exp}_${vname}.nc
cdo selyear,1950/2100 ${ncname} ${detrendpath}/temp.nc
cdo ymonmean ${detrendpath}/temp.nc ${scyclefile}
cdo ymonsub ${detrendpath}/temp.nc ${scyclefile} ${detrendpath}/temp1.nc
cdo detrend ${detrendpath}/temp1.nc ${outfile}

# 31km SSP5.85
exp="TCo319_ssp585"
vname="lcc"
ncname="/export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/TCo319_ssp585_lcc_1m_2015-2114_atmogrid.nc"
detrendpath="/home/niu4/gliu8/projects/scrap/processed_global/global_anom_detrend1"
scyclepath="/home/niu4/gliu8/projects/scrap/processed_global/scycle"
infile=${ncname}
scyclefile=${scyclepath}/${exp}_${vname}.nc
outfile=${detrendpath}/${exp}_${vname}.nc
cdo ymonmean ${infile} ${scyclefile}
cdo ymonsub ${infile} ${scyclefile} ${detrendpath}/temp1.nc
cdo detrend ${detrendpath}/temp1.nc ${outfile}




# ======================================================


