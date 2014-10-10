from cosmosis import names



def setup(options):
	return 1


def execute(block, config):
	sums = {}
	nbin = None

	#Load the GG
	nbin = block[names.shear_cl_gg, "nbin"]
	ell = block[names.shear_cl_gg, "ell"]
	for b1 in xrange(1,nbin+1):
		for b2 in xrange(1,b1+1):
			col = "bin_%d_%d"%(b1,b2)
			sums[(b1,b2)] = block[names.shear_cl_gg, col]

	#Add in the II term
	for b1 in xrange(1,nbin+1):
		for b2 in xrange(1,b1+1):
			col = "bin_%d_%d"%(b1,b2)
			sums[(b1,b2)] += block[names.shear_cl_ii, col]

	#Add in the GI term.  This one is different as there is
	#a contribution from the bin pair in both directions
	for b1 in xrange(1,nbin+1):
		for b2 in xrange(1,b1+1):
			col = "bin_%d_%d"%(b1,b2)
			sums[(b1,b2)] += block[names.shear_cl_gi, col]
			col = "bin_%d_%d"%(b2,b1)
			sums[(b1,b2)] += block[names.shear_cl_gi, col]

	block[names.shear_cl, "ell"] = ell
	block[names.shear_cl, "nbin"] = nbin

	for b1 in xrange(1,nbin+1):
		for b2 in xrange(1,b1+1):
			col = "bin_%d_%d"%(b1,b2)
			block[names.shear_cl, col] = sums[(b1,b2)]

	return 0

def cleanup(config):
	pass