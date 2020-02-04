#!/usr/bin/env python3

import numpy as np
import pickle

# FIRST, I will create the catalog file
f = open("data/bpz_catalog.cat",'w')

names, files = np.loadtxt('filters/filters.list', unpack=True, dtype=str)

# first, write the header
header = "#  1 id\n"
colnames = " id"
for i,name in enumerate(names):
    i = 2*i+2
    header += "# {0:>2} {1}\n# {2:>2} d{1}\n".format(i,name,i+1)
    colnames += "{0:>9}{1:>9}".format(name,"d"+name)
header += "# {0:>2} zspec\n#\n#".format(i+2)
colnames += "{0:>8}".format("zspec")
header += colnames
header += '\n#\n'

f.write(header)

# now write the galaxy catalog
with open('data/galaxy_catalog.pkl', 'rb') as input:
    galaxies = pickle.load(input)

for i,galaxy in enumerate(galaxies):

    # save the mags we have, set the rest to -99 (err = 0)
    idx = np.where(np.in1d(names,galaxy.filters))
    mags = -99*np.ones(len(names))
    mags[idx] = galaxy.mags
    errs = np.zeros(len(names))
    errs[idx] = galaxy.mag_err

    # interweave all the mags and mag errs
    output = np.empty((mags.size + errs.size,), dtype=mags.dtype)
    output[0::2] = mags
    output[1::2] = errs

    # save the line which is id, mags and errs, specz
    line = ("{:>4}"+"{:>9.3f}"*len(output)+"{:>8.3f}\n").format(i,*output,galaxy.redshift)
    f.write(line)

f.close()


# SECOND, I will create the columns file
f = open("data/bpz_catalog.columns",'w')
f.write("# Filter    columns  AB/Vega  zp_error  zp_offset\n")
for i,name in enumerate(files):
    name = name[:-4]
    col1 = 2*i+2
    col2 = 2*i+3
    f.write("{:<14}{:>2},{:>2}{:>8}{:>11.2f}{:>11.2f}\n".format(name,col1,col2,'AB',0,0))
icol = 2*np.where(names == 'i')[0][0] + 2
f.write("{:<14}{:>5}\n".format("M_0",icol))
f.write("{:<14}{:>5}\n".format("Z_S",2*len(names)+2))
f.write("{:<14}{:>5}".format("ID",1))
f.close()