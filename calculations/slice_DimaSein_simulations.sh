
# N39 NAO Slice (msl)
cdo sellonlatbox,270,40,20,80 /work/uo0119/a270067/runtime/awicm3-v3.1/N39/outdata/oifs/N39_1m_msl_1950_2014.nc /home/k/k206179/scrap/nao_crop/temp.nc
cdo selyear,1950/1979 /home/k/k206179/scrap/nao_crop/temp.nc DimaSein_N39_msl.nc
rm /home/k/k206179/scrap/nao_crop/temp.nc

# N44 NAO Slice (msl)
cdo sellonlatbox,270,40,20,80 /work/uo0119/a270067/runtime/awicm3-v3.1/N44/outdata/oifs/N44_1m_msl_1980_2014.nc /home/k/k206179/scrap/nao_crop/DimaSein_N44_msl.nc

# N43 NAO Slice (msl)
cdo sellonlatbox,270,40,20,80 /work/uo0119/a270067/runtime/awicm3-v3.1/N43/outdata/oifs/N43_1m_sst_2015_2100.nc /home/k/k206179/scrap/nao_crop/DimaSein_N43_msl.nc

# N39 SST slice
cdo sellonlatbox,120,-70,-20,20 /work/uo0119/a270067/runtime/awicm3-v3.1/N39/outdata/oifs/N39_1m_sst_1950_2014.nc /home/k/k206179/scrap/tropical_crop/temp.nc
cdo selyear,1950/1979 /home/k/k206179/scrap/tropical_crop/temp.nc DimaSein_N39_sst.nc
rm /home/k/k206179/scrap/tropical_crop/temp.nc

# N44 SST slice
cdo sellonlatbox,120,-70,-20,20 /work/uo0119/a270067/runtime/awicm3-v3.1/N44/outdata/oifs/N44_1m_sst_1980_2014.nc /home/k/k206179/scrap/tropical_crop/DimaSein_N44_sst.nc

# N43 SST slice
cdo sellonlatbox,120,-70,-20,20 /work/uo0119/a270067/runtime/awicm3-v3.1/N43/outdata/oifs/N43_1m_sst_2015_2100.nc /home/k/k206179/scrap/tropical_crop/DimaSein_N43_sst.nc

# ------------------------------------------------------



# ------------------------------------------------------

# Process MSL
expnames=("DimaSein_N44" "DimaSein_N39" "DimaSein_N43")
dpath="/home/k/k206179/scrap/nao_crop"
vname="msl"
for exp in ${expnames[@]}; do # This loop format is for zsh. Use ${expes[@]} if you are using bash.
    infile=${dpath}/${exp}_${vname}.nc
    scyclefile=${dpath}/${exp}_${vname}_scycle.nc
    outfile=${dpath}/${exp}_${vname}_anom.nc
    cdo ymonmean ${infile} ${scyclefile}
    cdo ymonsub ${infile} ${scyclefile} ${dpath}/temp1.nc
    cdo detrend ${dpath}/temp1.nc ${outfile}
    echo "Completed $vname for $exp"
done

# ------------------------------------------------------

# Process SST
expnames=("DimaSein_N44" "DimaSein_N39" "DimaSein_N43")
dpath="/home/k/k206179/scrap/tropical_crop"
vname="sst"
for exp in ${expnames[@]}; do # This loop format is for zsh. Use ${expes[@]} if you are using bash.
    infile=${dpath}/${exp}_${vname}.nc
    scyclefile=${dpath}/${exp}_${vname}_scycle.nc
    outfile=${dpath}/${exp}_${vname}_anom.nc
    cdo ymonmean ${infile} ${scyclefile}
    cdo ymonsub ${infile} ${scyclefile} ${dpath}/temp1.nc
    cdo detrend ${dpath}/temp1.nc ${outfile}
    echo "Completed $vname for $exp"
done