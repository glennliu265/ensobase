# Consolidate netcdfs downloaded from aleph
# Move to project folder for analysis
cdo mergetime /home/niu4/gliu8/share/awicm3/TCo319-DART-ctl1950d-gibbs-charn/oifs/sst/*.nc /home/niu4/gliu8/share/awicm3/TCo319-DART-ctl1950d-gibbs-charn/merged/atm_remapped_1m_sst_1m_1950-2004.nc

# Test Version (with 31km Control sst)
# Merge Files
cdo -O -mergetime -apply,-selname,sst [ /home/niu4/gliu8/share/awicm3/TCo319-DART-ctl1950d-gibbs-charn/oifs/sst/*.nc ] /home/niu4/gliu8/share/awicm3/TCo319-DART-ctl1950d-gibbs-charn/merged/atm_remapped_1m_sst_1m_1950-2004.nc 
# Move the Processed Folder
cp /home/niu4/gliu8/share/awicm3/TCo319-DART-ctl1950d-gibbs-charn/merged/atm_remapped_1m_sst_1m_1950-2004.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo319-DART-ctl1950d-gibbs-charn_sst.nc

# Generalized Version (note this includes selecting sst, dropping time bnds, will this cause issues?)
vname="sst"
expname="TCo1279-DART-2060" # "TCo319-DART-ssp585d-gibbs-charn", "TCo319-DART-ctl1950d-gibbs-charn", 
compname="oifs"
cdo -O -mergetime -apply,-selname,sst [ /home/niu4/gliu8/share/awicm3/${expname}/${compname}/1m/${vname}/*.nc ] /home/niu4/gliu8/projects/scrap/processed_global/${expname}_${vname}.nc 

# Finalized Loops Below!!

# Loop for "TCo319-DART-ctl1950d-gibbs-charn"
expname="TCo319-DART-ctl1950d-gibbs-charn"
compname="oifs"
vnames=("cp" "lsp" "lcc" "tcc" "slhf" "sshf" "slhf" "ssr" "str" "tsr" "ttr" "tsrc" "ttrc" "msl")
for vname in ${vnames[@]}; do
    cdo -O -mergetime [ /home/niu4/gliu8/share/awicm3/${expname}/${compname}/1m/${vname}/*.nc ] /home/niu4/gliu8/projects/scrap/processed_global/${expname}_${vname}.nc 
done

# Loop for "TCo319-DART-ssp585d-gibbs-charn"
expname="TCo319-DART-ssp585d-gibbs-charn"
compname="oifs"
vnames=("cp" "lsp" "lcc" "tcc" "slhf" "sshf" "slhf" "ssr" "str" "tsr" "ttr" "tsrc" "ttrc" "msl" "ci")
for vname in ${vnames[@]}; do
    cdo -O -mergetime [ /home/niu4/gliu8/share/awicm3/${expname}/${compname}/1m/${vname}/*.nc ] /home/niu4/gliu8/projects/scrap/processed_global/${expname}_${vname}.nc 
done

# Loop for "TCo1279-DART-2060"
expname="TCo1279-DART-2060"
compname="oifs"
vnames=("cp" "lsp" "lcc" "tcc" "slhf" "sshf" "slhf" "ssr" "str" "tsr" "ttr" "tsrc" "ttrc" "msl" "ci")
for vname in ${vnames[@]}; do
    cdo -O -mergetime [ /home/niu4/gliu8/share/awicm3/${expname}/${compname}/1m/${vname}/*.nc ] /home/niu4/gliu8/projects/scrap/processed_global/${expname}_${vname}.nc 
done


# Loop for "TCo95-hi1950d"
expname="TCo95-hi1950d"
compname="oifs"
vnames=("sst" "msl")
for vname in ${vnames[@]}; do
    cdo -O -mergetime [ /home/niu4/gliu8/share/awicm3/${expname}/${compname}/1m/${vname}/*.nc ] /home/niu4/gliu8/projects/scrap/processed_global/${expname}_${vname}.nc 
done

# Loop for "TCo95-ssp585d"
expname="TCo95-ssp585d"
compname="oifs"
vnames=("sst" "msl")
for vname in ${vnames[@]}; do
    cdo -O -mergetime [ /home/niu4/gliu8/share/awicm3/${expname}/${compname}/1m/${vname}/*.nc ] /home/niu4/gliu8/projects/scrap/processed_global/${expname}_${vname}.nc 
done


# Loop for monthly data, TCo319_ctl1950d
expname="TCo319_ctl1950d"
compname="oifs"
vnames=("10u" "10v" "10ws")
for vname in ${vnames[@]}; do
    cdo -O -mergetime [ /home/niu4/gliu8/share/awicm3/TCo319-ctl1950d/${compname}/1m/${vname}/*.nc ] /home/niu4/gliu8/projects/scrap/processed_global/${expname}_${vname}.nc 
done


# Consolidate for another experiment
expnames=("TCo319-DART-hi1950d-gibbs-charn" "TCo319-DART-ctl1950d-gibbs-charn" "TCo1279-DART-2060") #("TCo319-DART-ssp585d-gibbs-charn")
vnames=( "msl" )
compname="oifs"
freq="1m"
outpath="/home/niu4/gliu8/projects/scrap/processed_global"
for expname in ${expnames[@]}; do
    for vname in ${vnames[@]}; do
        cdo -O -mergetime [ /home/niu4/gliu8/share/awicm3/${expname}/${compname}/${freq}/${vname}/*.nc ] ${outpath}/${expname}_${vname}.nc 
    done
done

# Also Copy MSL over from other runs
cp /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_ssp585/TCo319_ssp585_msl_1m_2015-2114.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ssp585_msl.nc
cdo selyear,1950/2134 /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo319_control/TCo319_ctl1950d_msl_1m_1850-2134.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo319_ctl1950d_msl.nc
cp /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-1950_9km/TCo1279_DART-1950_msl_1m_1950-1969.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-1950_msl.nc
cp /export/niu2/stuecker/MODELOUTPUT/awicm3_highres/TCo1279-DART-2090_9km/TCo1279_DART-2090_msl_1m_2090-2099.nc /home/niu4/gliu8/projects/scrap/processed_global/TCo1279-DART-2090_msl.nc
cp 

/home/niu4/gliu8/share/awicm3/TCo319-DART-ssp585d-gibbs-charn/oifs/1m/msl


# Rename Files


# Loop it and move to project folder
#vnames = ("")
#cdo mergetime /home/niu4/gliu8/share/awicm3/TCo319-DART-ctl1950d-gibbs-charn/oifs/1m/sst/*.nc /home/niu4/gliu8/share/awicm3/TCo319-DART-ctl1950d-gibbs-charn/merged/atm_remapped_1m_sst_1m_1950-2004.nc



