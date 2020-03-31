
NOTE: This is out of date and needs to be cleaned up!

Repo that contains all the code and plots for the template learning paper about my algorithm to learn photo-z templates from broadband galaxy photometry.


Workflow:

1. <code>catalog.ipynb</code> - this notebook opens the DEEP2/3 and HST spec-z catalogs that are in the data folder. It applies quality cuts, combines the qualifying galaxies into a single catalog for my use, and also makes plots of the data set.

2. <code>training_example.ipynb</code> - this notebook takes the 5Myr starburst template that comes with BPZ, simulates observations, and then trains a flat template until it matches the data. The classes for filters and seds, and the all the functions used in the learning algorithm are found in <code>modules.py</code>.

3. <code>training_naive.ipynb</code> - this takes the combined galaxy catalog, establishes naive templates, and trains them to match the data. Matching and training are iterated until templates approximately converge. The notebook also contains plots analyzing the final trained templates.

4. <code>training_cwwsb4.ipynb</code> - same as number 3, except using the CWW+SB4 templates as the starting point.

5. <code>create_bpz_catalog.py</code> - this script takes the combined catalog from number 1, and creates a catalog file that is readable by BPZ.

6. <code>bpz_script.py</code> - this script runs BPZ to estimate photo-z's using the 3 template sets: CWW+SB4 original, CWW+SB4 trained, and naive trained. It runs the original bpz script, then the bpz finalize script, and then cleans up after BPZ. Note that all BPZ output is supressed. The script requires that you have BPZ installed on your system at <code>BPZPATH</code>, which is defined at the top of the script, and that you are in an environment with python 2. BPZ can be installed [here](http://www.stsci.edu/~dcoe/BPZ/).

7. <code>photoz_analysis.ipynb</code> - this notebook loads all of the BPZ photo-z results, and makes plots/computes stats to analyze them.

Note all of the dependencies can be installed via <code>pip install --user --requirement requirements.txt</code>. Python version is 3.7.3, except BPZ, which was run in version 2.7.17. Jupyter notebooks can be opened via running <code>jupyter lab</code> in the root directory, and then double clicking on notebooks in the sidebar.
