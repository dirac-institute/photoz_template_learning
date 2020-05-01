#!/usr/bin/env python3

import numpy as np
import pickle

# This script creates the catalog and column files for BPZ 

# load the validation catalog
with open('data/test_catalog.pkl', 'rb') as input:
    allgalaxies = pickle.load(input)

# I need to make catalog and column files for each of the various i-bands I am
# using for the magnitude prior
ibands = ['i','i2','Icfh12k','i+']
for i,band in enumerate(ibands):
    
    # get the galaxies that will be using this i band
    galaxies = np.array([galaxy for galaxy in allgalaxies
                            if band in galaxy.filters and
                            set(galaxy.filters).isdisjoint(ibands[:i])])
    if len(galaxies) == 0:
        continue

    # create the catalog file
    f = open("bpz_files/bpz_catalog_"+band+".cat",'w')

    # load the names of the filters
    names, files = np.loadtxt('filters/filters.list', unpack=True, dtype=str)

    # write the header
    header = "#  1 id\n"
    colnames = " id"
    for j,name in enumerate(names):
        j = 2*j+2
        header += "# {0:>2} {1}\n# {2:>2} d{1}\n".format(j,name,j+1)
        colnames += "{0:>9}{1:>9}".format(name,"d"+name)
    header += "# {0:>2} zspec\n#\n#".format(j+2)
    colnames += "{0:>8}".format("zspec")
    header += colnames
    header += '\n#\n'

    f.write(header)

    # write the rows for the galaxy data
    for j,galaxy in enumerate(galaxies):

        # save the mags we have, set the rest to -99 (err = 0)
        idx = np.where( np.in1d(names,galaxy.filters) )
        mags = -99*np.ones(len(names))
        mags[idx] = galaxy.mags
        errs = np.zeros(len(names))
        errs[idx] = galaxy.mag_err

        # interweave all the mags and mag errs
        output = np.empty((mags.size + errs.size,), dtype=mags.dtype)
        output[0::2] = mags
        output[1::2] = errs

        # save the line which is id, mags and errs, specz
        line = ("{:>4}"+"{:>9.3f}"*len(output)+"{:>8.3f}\n").format(j,*output,galaxy.redshift)
        f.write(line)

    f.close()
    print("Saving","bpz_files/bpz_catalog_"+band+".cat")


    # create the columns file
    # these are the same for each, except for the M0 column number
    f = open("bpz_files/bpz_catalog_"+band+".columns",'w')
    f.write("# Filter    columns  AB/Vega  zp_error  zp_offset\n")
    for j,name in enumerate(files):
        name = name[:-4]
        col1 = 2*j+2
        col2 = 2*j+3
        f.write("{:<14}{:>2},{:>2}{:>8}{:>11.2f}{:>11.2f}\n".format(name,col1,col2,'AB',0,0))
    icol = 2*np.where(names == band)[0][0] + 2
    f.write("{:<14}{:>5}\n".format("M_0",icol))
    f.write("{:<14}{:>5}\n".format("Z_S",2*len(names)+2))
    f.write("{:<14}{:>5}".format("ID",1))
    f.close()
    print("Saving","bpz_files/bpz_catalog_"+band+".columns")