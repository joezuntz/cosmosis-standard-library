from cosmosis.datablock import option_section, names

def setup(options):

	shear_shear=options.get_bool(option_section,"shear_shear",True)
	pos_shear=options.get_bool(option_section,"pos_shear",False)
	sec_names={}
	sec_names['GG']=options.get_string(option_section,"cl_shear_shear",names.shear_cl_gg)
	sec_names['GI']=options.get_string(option_section,"cl_shear_int",names.shear_cl_gi)
	sec_names['II']=options.get_string(option_section,"cl_int_int",names.shear_cl_ii)
	sec_names['gI']=options.get_string(option_section,"cl_gal_int","galaxy_intrinsic_cl")
	sec_names['gG']=options.get_string(option_section,"cl_gal_shear","galaxy_shear_cl")
	sec_names['GG_out']=options.get_string(option_section,"cl_shear_shear_out",names.shear_cl)
	sec_names['gG_out']=options.get_string(option_section,"cl_gal_shear_out","galaxy_shear_cl")
	return shear_shear,pos_shear,sec_names

def execute(block, config):

	shear_shear,pos_shear,sec_names=config
	
	if pos_shear:
		pos_shear_total={}
		ell = block[sec_names['gG'], "ell"]
		nbin_A,nbin_B=block[sec_names['gG'],"nbin_A"],block[sec_names['gG'],"nbin_B"]
		for b1 in xrange(1,nbin_A+1):
			for b2 in xrange(1,nbin_B+1):
				col = "bin_%d_%d"%(b1,b2)
				cl = block[sec_names['gG'], col]
				pos_shear_total[(b1,b2)] = cl

		for b1 in xrange(1,nbin_A+1):
			for b2 in xrange(1,nbin_B+1):
				col = "bin_%d_%d"%(b1,b2)
				pos_shear_total[(b1,b2)] += block[sec_names['gI'], col]

		block[sec_names['gG_out'], "ell"] = ell
		block[sec_names['gG_out'], "nbin_A"] = nbin_A
		block[sec_names['gG_out'], "nbin_B"] = nbin_A

		for b1 in xrange(1,nbin_A+1):
			for b2 in xrange(b1,nbin_B+1):
				col = "bin_%d_%d"%(b1,b2)
				block[sec_names['gG_out'], col] = pos_shear_total[(b1,b2)]

	if shear_shear:
		shear_shear_total = {}
		nbin = None

		#Load the GG
		try:
			nbin = block[sec_names['GG'], "nbin"]
		except:
			nbin = block[sec_names['GG'], "nbin_A"]
		ell = block[sec_names['GG'], "ell"]
		for b1 in xrange(1,nbin+1):
			for b2 in xrange(b1,nbin+1):
				col = "bin_%d_%d"%(b1,b2)
				cl = block[sec_names['GG'], col]
				shear_shear_total[(b1,b2)] = cl

		#Add in the II term
		for b1 in xrange(1,nbin+1):
			for b2 in xrange(b1,nbin+1):
				col = "bin_%d_%d"%(b1,b2)
				shear_shear_total[(b1,b2)] += block[sec_names['II'], col]

		#Add in the GI term.  This one is different as there is
		#a contribution from the bin pair in both directions
		for b1 in xrange(1,nbin+1):
			for b2 in xrange(b1,nbin+1):
				col = "bin_%d_%d"%(b1,b2)
				shear_shear_total[(b1,b2)] += block[sec_names['GI'], col]
				col = "bin_%d_%d"%(b2,b1)
				shear_shear_total[(b1,b2)] += block[sec_names['GI'], col]

		block[sec_names['GG_out'], "ell"] = ell
		block[sec_names['GG_out'], "nbin"] = nbin

		for b1 in xrange(1,nbin+1):
			for b2 in xrange(b1,nbin+1):
				col = "bin_%d_%d"%(b1,b2)
				col_rev = "bin_%d_%d"%(b2,b1)
				block[sec_names['GG_out'], col] = shear_shear_total[(b1,b2)]
				block[sec_names['GG_out'], col_rev] = shear_shear_total[(b1,b2)]
		block[sec_names['GG_out'],'nbin_A'],block[sec_names['GG_out'],'nbin_B']=nbin,nbin
	return 0

def cleanup(config):
	pass
