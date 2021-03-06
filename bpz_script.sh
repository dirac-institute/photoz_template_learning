
# settings
export CONDALOC=/home/jfcrenshaw/miniconda3/etc/profile.d/conda.sh
export PY2ENV=py2
export BPZPATH=$HOME/documents/bpz-1.99.3
export OUTDIR=bpz_files
export IBANDS=('i' 'i2' 'Icfh12k' 'i+')
export NINTERP=0 # number of templates to interpolate for BPZ

# Number of each galaxy type in the CWW+SB4 trained templates
export NEl_8cw=1 # number of eliptical galaxy templates 
export NSp_8cw=2 # number of spiral galaxy templates
export NIS_8cw=5 # number of irregular/star burst galaxy templates

# Number of each galaxy type in the naive set with 8 trained templates
export NEl_8=1 # number of eliptical galaxy templates 
export NSp_8=4 # number of spiral galaxy templates
export NIS_8=3 # number of irregular/star burst galaxy templates

# Number of each galaxy type in the naive set with 16 trained templates
export NEl_16=2 # number of eliptical galaxy templates
export NSp_16=8 # number of spiral galaxy templates
export NIS_16=6 # number of irregular/star burst galaxy templates

STARTTIME=$(date +%s)

# print reminders
echo " "
echo "NOTE:" 
echo "Make sure the settings at the top of the script are correct."
echo "This includes:"
echo "  -location of the conda bash function"
echo "  -name of the python 2 environment"
echo "  -location of BPZ"
echo "  -where to save BPZ results"
echo "  -list of i-bands used for the magnitude prior"
echo "  -number of templates to interpolate in BPZ"
echo "  -number of each spectral type in the trained naive template set (determined"
echo "      in color_classify.ipynb)"
echo "  "

# switch to environment with python 2 
source $CONDALOC
conda activate $PY2ENV

# copy templates and filters to the BPZ folders
cp templates/* $BPZPATH/SED/
cp filters/*res $BPZPATH/FILTER/
cp prior_calibrated.py $BPZPATH/


echo "Running BPZ on CWW+SB4 templates..."
export OUTFILE=$OUTDIR/cwwsb4_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA cwwsb4.list -INTERP $NINTERP -PRIOR calibrated -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/cwwsb4_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/cwwsb4_$BAND\_photoz.bpz
done
echo " "

echo "Running BPZ on trained CWW+SB4 templates..."
export OUTFILE=$OUTDIR/cwwsb4_trained_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA cwwsb4_trained.list -INTERP $NINTERP -PRIOR calibrated -NTYPES $NEl_8cw $NSp_8cw $NIS_8cw -NEW_AB no -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/cwwsb4_trained_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/cwwsb4_trained_$BAND\_photoz.bpz
done
echo " "

#echo "Running BPZ on the 8 trained naive templates..."
#export OUTFILE=$OUTDIR/N8_trained_output.txt
#rm $OUTFILE 2> /dev/null
#echo "Saving output to" $OUTFILE"..."
#for BAND in "${IBANDS[@]}"; do
#    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA N8_trained.list -INTERP $NINTERP -PRIOR calibrated -NTYPES $NEl_8 $NSp_8 $NIS_8 -NEW_AB no -VERBOSE no &>> $OUTFILE
#    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
#    echo "Saving" $OUTDIR"/N8_trained_"$BAND"_photoz.bpz..."
#    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/N8_trained_$BAND\_photoz.bpz
#done
#echo " "

echo "Running BPZ on the 16 trained naive templates..."
export OUTFILE=$OUTDIR/N16_trained_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA N16_trained.list -INTERP $NINTERP -PRIOR calibrated -NTYPES $NEl_16 $NSp_16 $NIS_16 -NEW_AB no -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/N16_trained_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/N16_trained_$BAND\_photoz.bpz
done
echo " "


rm $OUTDIR/*catalog*bpz 2> /dev/null
rm $OUTDIR/*bak 2> /dev/null
rm $OUTDIR/*flux_comparison 2> /dev/null
rm $OUTDIR/*probs 2> /dev/null


ENDTIME=$(date +%s)
echo "Duration $(($ENDTIME - $STARTTIME)) seconds"
echo " "