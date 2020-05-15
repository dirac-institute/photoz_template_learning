# settings
export CONDALOC=/home/jfcrenshaw/miniconda3/etc/profile.d/conda.sh
export PY2ENV=py2
export BPZPATH=$HOME/documents/bpz-1.99.3
export OUTDIR=bpz_files
export IBANDS=('i' 'i2' 'Icfh12k' 'i+')
export NINTERP=0 # number of templates to interpolate for BPZ

# Number of each galaxy type in various sets

export NEl_6=1 # number of eliptical galaxy templates 
export NSp_6=2 # number of spiral galaxy templates
export NIS_6=3 # number of irregular/star burst galaxy templates

export NEl_7=1 # number of eliptical galaxy templates 
export NSp_7=3 # number of spiral galaxy templates
export NIS_7=3 # number of irregular/star burst galaxy templates

export NEl_8=1 # number of eliptical galaxy templates 
export NSp_8=4 # number of spiral galaxy templates
export NIS_8=3 # number of irregular/star burst galaxy templates

export NEl_9=1 # number of eliptical galaxy templates 
export NSp_9=4 # number of spiral galaxy templates
export NIS_9=4 # number of irregular/star burst galaxy templates

export NEl_10=1 # number of eliptical galaxy templates 
export NSp_10=5 # number of spiral galaxy templates
export NIS_10=4 # number of irregular/star burst galaxy templates

export NEl_11=1 # number of eliptical galaxy templates 
export NSp_11=5 # number of spiral galaxy templates
export NIS_11=5 # number of irregular/star burst galaxy templates

export NEl_12=1 # number of eliptical galaxy templates 
export NSp_12=5 # number of spiral galaxy templates
export NIS_12=6 # number of irregular/star burst galaxy templates

export NEl_13=1 # number of eliptical galaxy templates 
export NSp_13=6 # number of spiral galaxy templates
export NIS_13=6 # number of irregular/star burst galaxy templates

export NEl_14=1 # number of eliptical galaxy templates 
export NSp_14=7 # number of spiral galaxy templates
export NIS_14=6 # number of irregular/star burst galaxy templates

export NEl_16=2 # number of eliptical galaxy templates 
export NSp_16=8 # number of spiral galaxy templates
export NIS_16=6 # number of irregular/star burst galaxy templates

export NEl_18=2 # number of eliptical galaxy templates 
export NSp_18=8 # number of spiral galaxy templates
export NIS_18=8 # number of irregular/star burst galaxy templates

export NEl_20=2 # number of eliptical galaxy templates 
export NSp_20=9 # number of spiral galaxy templates
export NIS_20=9 # number of irregular/star burst galaxy templates

export NEl_22=2 # number of eliptical galaxy templates 
export NSp_22=10 # number of spiral galaxy templates
export NIS_22=10 # number of irregular/star burst galaxy templates

export NEl_24=3 # number of eliptical galaxy templates 
export NSp_24=10 # number of spiral galaxy templates
export NIS_24=11 # number of irregular/star burst galaxy templates



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


echo "Running BPZ on the 6 trained templates..."
export OUTFILE=$OUTDIR/N6_Ntemplate_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA N6_trained.list -INTERP $NINTERP -NTYPES $NEl_6 $NSp_6 $NIS_6 -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/N6_Ntemplates_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/N6_Ntemplates_$BAND\_photoz.bpz
done
echo " "

echo "Running BPZ on the 7 trained templates..."
export OUTFILE=$OUTDIR/N7_Ntemplate_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA N7_trained.list -INTERP $NINTERP -NTYPES $NEl_7 $NSp_7 $NIS_7 -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/N7_Ntemplates_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/N7_Ntemplates_$BAND\_photoz.bpz
done
echo " "

echo "Running BPZ on the 8 trained templates..."
export OUTFILE=$OUTDIR/N8_Ntemplate_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA N8_trained.list -INTERP $NINTERP -NTYPES $NEl_8 $NSp_8 $NIS_8 -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/N8_Ntemplates_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/N8_Ntemplates_$BAND\_photoz.bpz
done
echo " "

echo "Running BPZ on the 9 trained templates..."
export OUTFILE=$OUTDIR/N9_Ntemplate_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA N9_trained.list -INTERP $NINTERP -NTYPES $NEl_9 $NSp_9 $NIS_9 -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/N9_Ntemplates_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/N9_Ntemplates_$BAND\_photoz.bpz
done
echo " "

echo "Running BPZ on the 10 trained templates..."
export OUTFILE=$OUTDIR/N10_Ntemplate_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA N10_trained.list -INTERP $NINTERP -NTYPES $NEl_10 $NSp_10 $NIS_10 -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/N10_Ntemplates_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/N10_Ntemplates_$BAND\_photoz.bpz
done
echo " "

echo "Running BPZ on the 11 trained templates..."
export OUTFILE=$OUTDIR/N11_Ntemplate_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA N11_trained.list -INTERP $NINTERP -NTYPES $NEl_11 $NSp_11 $NIS_11 -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/N11_Ntemplates_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/N11_Ntemplates_$BAND\_photoz.bpz
done
echo " "

echo "Running BPZ on the 12 trained templates..."
export OUTFILE=$OUTDIR/N12_Ntemplate_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA N12_trained.list -INTERP $NINTERP -NTYPES $NEl_12 $NSp_12 $NIS_12 -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/N12_Ntemplates_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/N12_Ntemplates_$BAND\_photoz.bpz
done
echo " "

echo "Running BPZ on the 13 trained templates..."
export OUTFILE=$OUTDIR/N13_Ntemplate_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA N13_trained.list -INTERP $NINTERP -NTYPES $NEl_13 $NSp_13 $NIS_13 -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/N13_Ntemplates_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/N13_Ntemplates_$BAND\_photoz.bpz
done
echo " "

echo "Running BPZ on the 14 trained templates..."
export OUTFILE=$OUTDIR/N14_Ntemplate_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA N14_trained.list -INTERP $NINTERP -NTYPES $NEl_14 $NSp_14 $NIS_14 -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/N14_Ntemplates_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/N14_Ntemplates_$BAND\_photoz.bpz
done
echo " "

echo "Running BPZ on the 16 trained templates..."
export OUTFILE=$OUTDIR/N16_Ntemplate_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA N16_trained.list -INTERP $NINTERP -NTYPES $NEl_16 $NSp_16 $NIS_16 -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/N16_Ntemplates_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/N16_Ntemplates_$BAND\_photoz.bpz
done
echo " "

echo "Running BPZ on the 18 trained templates..."
export OUTFILE=$OUTDIR/N18_Ntemplate_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA N18_trained.list -INTERP $NINTERP -NTYPES $NEl_18 $NSp_18 $NIS_18 -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/N18_Ntemplates_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/N18_Ntemplates_$BAND\_photoz.bpz
done
echo " "

echo "Running BPZ on the 20 trained templates..."
export OUTFILE=$OUTDIR/N20_Ntemplate_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA N20_trained.list -INTERP $NINTERP -NTYPES $NEl_20 $NSp_20 $NIS_20 -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/N20_Ntemplates_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/N20_Ntemplates_$BAND\_photoz.bpz
done
echo " "

echo "Running BPZ on the 22 trained templates..."
export OUTFILE=$OUTDIR/N22_Ntemplate_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA N22_trained.list -INTERP $NINTERP -NTYPES $NEl_22 $NSp_22 $NIS_22 -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/N22_Ntemplates_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/N22_Ntemplates_$BAND\_photoz.bpz
done
echo " "

echo "Running BPZ on the 24 trained templates..."
export OUTFILE=$OUTDIR/N24_Ntemplate_output.txt
rm $OUTFILE 2> /dev/null
echo "Saving output to" $OUTFILE"..."
for BAND in "${IBANDS[@]}"; do
    python $BPZPATH/bpz.py $OUTDIR/bpz_catalog_$BAND.cat -SPECTRA N24_trained.list -INTERP $NINTERP -NTYPES $NEl_24 $NSp_24 $NIS_24 -VERBOSE no &>> $OUTFILE
    python $BPZPATH/bpzfinalize.py $OUTDIR/bpz_catalog_$BAND &>> $OUTFILE
    echo "Saving" $OUTDIR"/N24_Ntemplates_"$BAND"_photoz.bpz..."
    mv $OUTDIR/bpz_catalog_$BAND\_bpz.cat $OUTDIR/N24_Ntemplates_$BAND\_photoz.bpz
done
echo " "


rm $OUTDIR/*catalog*bpz 2> /dev/null
rm $OUTDIR/*bak 2> /dev/null
rm $OUTDIR/*flux_comparison 2> /dev/null
rm $OUTDIR/*probs 2> /dev/null


ENDTIME=$(date +%s)
echo "Duration $(($ENDTIME - $STARTTIME)) seconds"
echo " "