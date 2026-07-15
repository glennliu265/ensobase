
dpath=/home/niu4/gliu8/projects/mesaclip/region_crops
rawnc=HiRes_BHISTC5_TS_Pacific_ens001_anom_detrended.nc
meanfile=${dpath}/mu.nc
outnc=${dpath}/skewtest2.nc

echo Processing file : ${dpath}/${rawnc}
echo Output : ${outnc}

# Calculate the Mean (Perhaps can skip if already anomalized...)
cdo -O -L -timmean ${dpath}/${rawnc} ${meanfile}

# Calculate the Numerator : timsum [ (x-mu)^3 ]  --> 366.79s
cdo -O -L -timsum -pow,3 -sub ${dpath}/${rawnc} ${meanfile} ${dpath}/skew_numer.nc

# Calculate the Denominator : stdev^3
cdo -O -L -pow,3 -timstd ${dpath}/${rawnc} ${dpath}/skew_denom.nc

# Calculate the Skewness (1/ntime) * (numer/denom)
ntime=$(cdo -ntime ${dpath}/${rawnc}) # Get Length of Time Dimension
cdo -O -L -divc,${ntime} -div ${dpath}/skew_numer.nc ${dpath}/skew_denom.nc ${outnc}

# Clean Files
rm ${dpath}/skew_denom.nc
rm ${dpath}/skew_numer.nc
rm ${meanfile}

#echo ${ntime}
#dpath=/home/niu4/gliu8/projects/mesaclip/region_crops
#rawnc=HiRes_BHISTC5_TS_Pacific_ens001_anom_detrended.nc
#cdo -O -timsum -pow,3 -sub ${dpath}/${rawnc} -timmean ${dpath}/${rawnc} ${dpath}/skew_numer.nc