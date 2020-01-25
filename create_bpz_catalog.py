#!/usr/bin/env python3

import numpy as np
from astropy.table import Table

data = Table.read('data/combined_catalog.fits')

data = np.array([np.arange(len(data)),
                 data['u'],  data['uerr'],
                 data['g'],  data['gerr'],
                 data['r'],  data['rerr'],
                 data['i2'], data['i2err'],
                 data['i'],  data['ierr'],
                 data['z'],  data['zerr'],
                 data['y'],  data['yerr'],
                 data['redshift']]).T

header = """#  1 id   
#  2 u    
#  3 du   
#  4 g   
#  5 dg   
#  6 r   
#  7 dr   
#  8 i2   
#  9 di2   
# 10 i    
# 11 di   
# 12 z    
# 13 dz
# 14 y    
# 15 dy
# 16 zspec
#
#\n"""


f = open("data/bpz_catalog.cat",'w')
f.write(header)
cols = ['id','u','du','g','dg','r','dr','i2','di2','i','di','z','dz','y','dy','zspec']
# id        u      du       g      dg       r      dr      i2     di2       i      di       z      dz       y      dy   zspec\n"""

f.write("# {0:>2}{1:>8}{2:>8}{3:>8}{4:>8}{5:>8}{6:>8}{7:>8}{8:>8}{9:>8}{10:>8}{11:>8}{12:>8}{13:>8}{14:>8}{15:>8}\n#\n".format(*cols))
for row in data:
    row = list(row)
    row[0] = int(row[0])
    for j,k in enumerate(row[1:]):
        row[j+1] = round(k,3)
    f.write("{0:> 4}{1:>8}{2:>8}{3:>8}{4:>8}{5:>8}{6:>8}{7:>8}{8:>8}{9:>8}{10:>8}{11:>8}{12:>8}{13:>8}{14:>8}{15:>8}\n".format(*row))
f.close()