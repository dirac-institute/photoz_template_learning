
# settings
export CONDALOC=/opt/miniconda3/etc/profile.d/conda.sh
export PY2ENV=py2
export BPZPATH=$HOME/documents/dirac/bpz-1.99.3
export OUTDIR=bpz_files
export IBANDS=('i')
export NINTERP=2 # number of templates to interpolate for BPZ

# Number of each galaxy type in the naive set with 8 templates
export NEl_8=1 # number of eliptical galaxy templates 
export NSp_8=5 # number of spiral galaxy templates
export NIS_8=2 # number of irregular/star burst galaxy templates

# Number of each galaxy type in the naive set with 20 templates
export NEl_20=2 # number of eliptical galaxy templates
export NSp_20=11 # number of spiral galaxy templates
export NIS_20=7 # number of irregular/star burst galaxy templates

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


echo "Running BPZ on CWW+SB4 templates..."
export OUTFILE=$OUTDIR/cwwsb4_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA cwwsb4.list -INTERP $NINTERP -VERBOSE no &>> $OUTFILE
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
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA cwwsb4_trained.list -INTERP $NINTERP -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/cwwsb4_trained_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/cwwsb4_trained_$BAND\_photoz.bpz
done
echo " "

echo "Running BPZ on the 8 trained naive templates..."
export OUTFILE=$OUTDIR/naive_8_trained_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA naive_8_trained.list -INTERP $NINTERP -NTYPES $NEl_8 $NSp_8 $NIS_8 -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/naive_8_trained_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/naive_8_trained_$BAND\_photoz.bpz
done
echo " "

echo "Running BPZ on the 20 trained naive templates..."
export OUTFILE=$OUTDIR/naive_20_trained_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA naive_20_trained.list -INTERP $NINTERP -NTYPES $NEl_20 $NSp_20 $NIS_20 -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/naive_20_trained_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/naive_20_trained_$BAND\_photoz.bpz
done
echo " "


rm $OUTDIR/*catalog*bpz 2> /dev/null
rm $OUTDIR/*bak 2> /dev/null
rm $OUTDIR/*flux_comparison 2> /dev/null
rm $OUTDIR/*probs 2> /dev/null


ENDTIME=$(date +%s)
echo "Duration $(($ENDTIME - $STARTTIME)) seconds"
echo " "