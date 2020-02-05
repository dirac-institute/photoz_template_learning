
STARTTIME=$(date +%s)

echo " "
echo "WARNING:" 
echo "Don't forget to switch to python 2. On my laptop: conda activate py2"
echo "You also might need to change the BPZPATH in the script"
echo " "

export BPZPATH=$HOME/documents/dirac/bpz-1.99.3

cp templates/* $BPZPATH/SED/
cp filters/*res $BPZPATH/FILTER/

echo "Running BPZ on CWW+SB4 templates..."
python $BPZPATH/bpz.py data/bpz_catalog.cat -SPECTRA cwwsb4.list -INTERP 2 > /dev/null
python $BPZPATH/bpzfinalize.py data/bpz_catalog > /dev/null
echo "Saving data/cwwsb4_photoz.bpz ..."
mv data/bpz_catalog_bpz.cat data/cwwsb4_photoz.bpz

echo "Running BPZ on trained CWW+SB4 templates..."
python $HOME/documents/dirac/bpz-1.99.3/bpz.py data/bpz_catalog.cat -SPECTRA cwwsb4_trained.list -INTERP 2 > /dev/null
python $HOME/documents/dirac/bpz-1.99.3/bpzfinalize.py data/bpz_catalog > /dev/null
echo "Saving data/cwwsb4_trained_photoz.bpz ..."
mv data/bpz_catalog_bpz.cat data/cwwsb4_trained_photoz.bpz

echo "Running BPZ on trained naive templates..."
python $HOME/documents/dirac/bpz-1.99.3/bpz.py data/bpz_catalog.cat -SPECTRA naive_trained.list -INTERP 2 -NTYPES 1 3 5 > /dev/null
python $HOME/documents/dirac/bpz-1.99.3/bpzfinalize.py data/bpz_catalog > /dev/null
echo "Saving data/naive_trained_photoz.bpz ..."
mv data/bpz_catalog_bpz.cat data/naive_trained_photoz.bpz

rm data/bpz_catalog.bpz
rm data/bpz_catalog.bpz.bak
rm data/bpz_catalog.flux_comparison
rm data/bpz_catalog.probs

echo " "

ENDTIME=$(date +%s)
echo "Duration $(($ENDTIME - $STARTTIME)) seconds"
echo " "