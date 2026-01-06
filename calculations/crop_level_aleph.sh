

module load cdo
module load netcdf

# Use CDO to crop single levels for data on cdo

cdo /scratch/awicm3-TCo319/TCo319-DART-ssp585d-gibbs-charn/outdata/oifs/atm_remapped_1m_pl_t_1m_pl_2015-2015.nc





cdo sellevel,70000 /scratch/awicm3-TCo319/TCo319-DART-ssp585d-gibbs-charn/outdata/oifs/atm_remapped_1m_pl_t_1m_pl_2015-2015.nc /home/gliu/processed_data/test700.nc



for i in $(ls /scratch/awicm3-TCo319/TCo319-DART-ssp585d-gibbs-charn/outdata/oifs/*pl_t_1m*.nc); do echo ${i}; done



# Crop 70000 hPa for (pl_t_700, 31km Future)
cd /scratch/awicm3-TCo319/TCo319-DART-ssp585d-gibbs-charn/outdata/oifs/ # Move to Directory
# Get list of files and operate
for i in $(ls *pl_t_1m*.nc); do cdo sellevel,70000 ${i} /home/gliu/processed_data/TCo319-DART-ssp585d-gibbs-charn/pl_t_700/${i}_700.nc; done



# Redo Loop (but for pl_w_700, 31km Future)
# Extract Filename based on based on https://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash
vname="pl_w_1m"
dpath="/scratch/awicm3-TCo319/TCo319-DART-ssp585d-gibbs-charn/outdata/oifs"
outpath="/home/gliu/processed_data/TCo319-DART-ssp585d-gibbs-charn/pl_w_700"
#mkdir outpath
flist=$(ls ${dpath}/*${vname}*.nc)
echo flist
for f in ${flist[@]}; do
    filename="${f##*/}" # Remove the full path (just get file name)
    filename="${filename%.*}" # Remove the extension
    outname=${outpath}/${filename}_700.nc
    cdo sellevel,70000 ${f} ${outname}
    echo $outname
done

# (pl_r_700, 31km Future)
vname="pl_r_1m"
outpath="/home/gliu/processed_data/TCo319-DART-ssp585d-gibbs-charn/pl_r_700"
dpath="/scratch/awicm3-TCo319/TCo319-DART-ssp585d-gibbs-charn/outdata/oifs"
flist=$(ls ${dpath}/*${vname}*.nc)
echo flist
for f in ${flist[@]}; do
    filename="${f##*/}" # Remove the full path (just get file name)
    filename="${filename%.*}" # Remove the extension
    outname=${outpath}/${filename}_700.nc
    cdo sellevel,70000 ${f} ${outname}
    echo $outname
done

# --- Above has been completed --- X

# (pl_t_700, 31km Control)
vname="pl_t_1m"
outpath="/home/gliu/processed_data/TCo319-DART-ctl1950d-gibbs-charn/pl_t_700"
dpath="/scratch/awicm3-TCo319/TCo319-DART-ctl1950d-gibbs-charn/outdata/oifs/"
flist=$(ls ${dpath}/*${vname}*.nc)
echo flist
for f in ${flist[@]}; do
    filename="${f##*/}" # Remove the full path (just get file name)
    filename="${filename%.*}" # Remove the extension
    outname=${outpath}/${filename}_700.nc
    cdo sellevel,70000 ${f} ${outname}
    echo $outname
done

# (pl_w_700, 31km Control)
module load cdo
vname="pl_w_1m"
outpath="/home/gliu/processed_data/TCo319-DART-ctl1950d-gibbs-charn/pl_w_700"
dpath="/scratch/awicm3-TCo319/TCo319-DART-ctl1950d-gibbs-charn/outdata/oifs"
flist=$(ls ${dpath}/*${vname}*.nc)
echo flist
for f in ${flist[@]}; do
    filename="${f##*/}" # Remove the full path (just get file name)
    filename="${filename%.*}" # Remove the extension
    outname=${outpath}/${filename}_700.nc
    cdo sellevel,70000 ${f} ${outname}
    echo $outname
done

# (pl_r_700, 31km Control)
module load cdo
vname="pl_r_1m"
outpath="/home/gliu/processed_data/TCo319-DART-ctl1950d-gibbs-charn/pl_r_700"
dpath="/scratch/awicm3-TCo319/TCo319-DART-ctl1950d-gibbs-charn/outdata/oifs"
flist=$(ls ${dpath}/*${vname}*.nc)
echo flist
for f in ${flist[@]}; do
    filename="${f##*/}" # Remove the full path (just get file name)
    filename="${filename%.*}" # Remove the extension
    outname=${outpath}/${filename}_700.nc
    cdo sellevel,70000 ${f} ${outname}
    echo $outname
done

#

# --- Below has not been Run --- X


