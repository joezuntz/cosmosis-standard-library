Each folder contains:

xis.txt - data measurement. First column is theta, then xi+(0,0),xi+(0,1),...,xi-(0,0),xi-(0,1)...
Cov.npy - covariance matrix (use np.load(<path>) to load as numpy array). Assumes the data vector is the columns (excluding the
	  first) of xis.txt, concatenated.
n_of_zs.hist - first column is z, subsequent columns are n(z) for each redshift bin.
