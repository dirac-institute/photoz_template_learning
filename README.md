
Repo that contains all the code and plots for the template learning paper about my algorithm to learn photo-z templates from broadband galaxy photometry.


The catalog notebook loads the spec-z catalogs that are in the data folder. It applies quality cuts, creates training and test sets, and creates a table and plots characterizing the data sets.

The training_example notebook simulates observations for the 5Myr starburst CWW+SB4 template, then applies the training algorithm to a flat template until it matches the original SED.

The training_N8, training_N16, and training_cwwsb4 apply the training algorithm to these sets repectively, using the training set created in the catalog notebook.

The N8_spectral_lines notebook performs post-processing on the N8 templates to reconstruct emission lines in the bluest templates.

The pca_analysis notebook does some basic pca analysis of the N16 and N8 templates. This wasn't used in the paper.

The N_templates notebook applies the training algorithm to sets of 6-24 templates. The goal is to see how the template number affects photo-z estimation.

create_bpz_catalog.py is a script that converts the test set into a format readable by BPZ.

The calibrate_prior notebook finds the parameters for the prior by calibrating to the training set.

prior_calibrated.py is the calibrated prior for BPZ. The parameters from the calibrate_prior notebook are hardcoded inside. This file is automatically copied to the BPZ folder by bpz_script.sh before running BPZ on the test set of galaxies.

bpz_script.sh is a shell script that runs BPZ on the CWW+SB4, trained CWW+SB4, N8, and N16 template sets. The resulting files are all saved in the bpz_files folder. This script runs the original bpz script, then the bpz finalize script, then cleans up the files I don't need. All of the bpz stdout and stderr are saved in an output file in the bpz_files folder. The script requires that you have BPZ installed on your system at BPZPATH, which is defined at the top of the script. BPZ can be installed [here](http://www.stsci.edu/~dcoe/BPZ/). You also need a python 2 environment. Check the top of this script for all the settings and variables to set.

Ntemplate_bpz_script.sh is a shell script that runs BPZ on the template sets created by the N_templates notebook. See the descripton of bpz_script.sh above for more details.

The photoz_analysis notebook analyzes the results of BPZ.

Note all of the dependencies can be installed via <code>pip install --user --requirement requirements.txt</code>. Python version is 3.7.3, except BPZ, which was run in version 2.7.17. Jupyter notebooks can be opened via running <code>jupyter lab</code> in the root directory, and then double clicking on notebooks in the sidebar.
