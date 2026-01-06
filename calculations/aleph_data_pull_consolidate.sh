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


# Loop it and move to project folder
#vnames = ("")
#cdo mergetime /home/niu4/gliu8/share/awicm3/TCo319-DART-ctl1950d-gibbs-charn/oifs/1m/sst/*.nc /home/niu4/gliu8/share/awicm3/TCo319-DART-ctl1950d-gibbs-charn/merged/atm_remapped_1m_sst_1m_1950-2004.nc



