The FITS files in this directory contain the DES SV shear-shear data in the 
format described in github.com/joezuntz/2point/ though it should be mostly obvious.

There are two tomographic ones, for the im3shape and ngmix analyses with 3 tomographic bins each.


I have not yet converted the two non-tomographic data sets to the new format; they are in directories in this old format:

xis.txt - data measurement. First column is theta, then xi+(0,0),xi+(0,1),...,xi-(0,0),xi-(0,1)...
Cov.npy - covariance matrix (use np.load(<path>) to load as numpy array). Assumes the data vector is the columns (excluding the
	  first) of xis.txt, concatenated.
n_of_zs.hist - first column is z, subsequent columns are n(z) for each redshift bin.
