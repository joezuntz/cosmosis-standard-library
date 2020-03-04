"""
This module loads data from the powtable files which summarize the simulation
results made by Hunjing Huang et al.

This is justified from OWLS module in cosmosis standard library.

I interpolate into that data using a bivariate spline to get an estimate of the
effect on the matter power from baryons at a given z and k.

I could have made mistakes here easily in understanding what is happening! 

"""
from astropy.table import Table
from builtins import range
from builtins import object
import numpy as np
import re
try:
    import scipy.interpolate
except ImportError:
    sys.stderr.write(
        "The eagle baryon power code requires the scipy python module\n")
    sys.stderr.write("but it could not be found.  This code will fail.\n")
    raise ImportError(
        "Require scipy module to run baryon code.  Could not find it.")

import os
dirname = os.path.split(__file__)[0]

# File names
        
class RatioDatPowerModulator():
    """ using hunjing's logpk ratio fits file. redshift list is specific to her file"""

    def __init__(self, ratiofile):
        #print(ratiofile)
        f = open(ratiofile, 'r')
        line = f.readline()
        f.close()
        colname_list = line.replace('\n','').split(' ')[1:]
        # assuming read-in redshift list is in decending order
        colname_list.reverse()
        #print(colname_list)
        redshift_list = np.zeros(len(colname_list))
        for i,col in enumerate(colname_list):
            zintstr = re.findall('\d+',col)[0]
            redshift_list[i]  = int(zintstr)/10.0**(len(zintstr)-1)
        #print("hydro sim redshifts:",redshift_list)
        data = np.loadtxt(ratiofile,skiprows=1)
        logk = data[:,0]
        #print((len(logk),len(redshift_list)))
        logpkratio = np.zeros((len(logk),len(redshift_list)))
        for i in range(len(redshift_list)):
            logpkratio[:,i] = data[:,-(i+1)]
        self.ratio = scipy.interpolate.RectBivariateSpline(logk, redshift_list, logpkratio)
        self.ratio_k = logk

    def modulate(self, k, z, P, r=None):
        logk = np.log10(k)
        s = self.ratio(logk,z)
        #s[k < self.ratio_k[0]] = 0.0
        modulation = 10**s
        print(modulation.shape,P.shape)
        P_out = P*modulation
        return P_out

