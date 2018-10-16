"""

Adapted from boss_rsd 
BOSS DR12 data are presented in Alam et al 2016, 1607.03155

"""

from numpy import log, pi, interp, where, loadtxt,dot, append, linalg
import os
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section

cosmo = section_names.cosmological_parameters
likes = section_names.likelihoods
growthparams  = section_names.growth_parameters 
dist = section_names.distances

ROOT_dir = os.path.split(os.path.abspath(__file__))[0]
RED_file = os.path.join(ROOT_dir,'final_consensus_redshift.txt')

c_km_per_s =  299792.458
default_rs_fiducial = 148.6


def setup(options):
	section = option_section
	mode = options.get_int(section, "mode", default=0)
	feedback = options.get_int(section, "feedback", default=0)
	if not mode:
		print("BAO only, Dm(z) and H(z)")
		#data and cov files
		DATA_file = os.path.join(ROOT_dir,'BAO_consensus_results_dM_Hz.txt')
		COV_file = os.path.join(ROOT_dir,'BAO_consensus_covtot_dM_Hz.txt')
	else:
		print("BAO+FS, 	Dm(z), H(z) and fsigma8")
		#data and cov files
		DATA_file = os.path.join(ROOT_dir,'final_consensus_results_dM_Hz_fsig.txt')
		COV_file = os.path.join(ROOT_dir,'final_consensus_covtot_dM_Hz_fsig.txt')

	data = loadtxt(DATA_file)
	cov = loadtxt(COV_file)
	redshift = loadtxt(RED_file)
	rs_fiducial = options.get_double(option_section, "rs_fiducial", default_rs_fiducial)

	return (mode,data,cov,redshift,rs_fiducial,feedback)

def execute(block, config):

	mode,data,cov,redshift,rs_fiducial,feedback = config

	#inversion of the covariance matrix
	precision = linalg.inv(cov)

	#theoretical computation
	dist_z = block[dist, 'z']
	rs = block[dist, "RS_ZDRAG"]
	dist_z = block[dist, 'z']
	d_a = block[dist, 'd_a']
	d_m = block[dist, 'd_m']
	h = c_km_per_s*block[dist, 'h']
	if dist_z[1] < dist_z[0]:
		dist_z = dist_z[::-1]
		d_a = d_a[::-1]
		d_m = d_m[::-1]
		h = h[::-1]

	h0 = block[cosmo, 'h0']
	
	# D_m and H interpolated at the 3 redshifts
	dm_z = interp(redshift,dist_z,d_m)
	h_z = interp(redshift,dist_z,h)
	
	#mult by r_s,fid/r_s
	dm_z_rs = dm_z * (rs_fiducial/rs)
	h_z_rs = h_z * (rs/rs_fiducial)

	#add f*sigma8 when mode == 1
	if mode:
		#redshift
		z = block[growthparams, 'z']
		#growth parameters
		d_z = block[growthparams, 'd_z']
		f_z = block[growthparams, 'f_z']
		sig = block[cosmo, 'sigma_8']
		try:
			z0 = where(z==0)[0][0]
		except IndexError:
			raise ValueError("You need to calculate f(z) and d(z) down to z=0 to use the BOSS f*sigma8 likelihood")
		
		# fsigma8 interpolated at the 3 redshifts
		fsigma = (sig*(d_z/d_z[z0]))*f_z
		fsig = interp(redshift, z, fsigma)

	
		if feedback:
			print("Growth parameters: z = ",redshift, "fsigma_8  = ",fsig, " z0 = ", z0)

	if not mode:
		#reordering the parameters
		params = [dm_z_rs,h_z_rs]
		params = append((params[0])[:],(params[1])[:])
		params = [params[0::3],params[1::3],params[2::3]]
		params = append(append((params[0])[:],(params[1])[:]),(params[2])[:])

	else:
		#reordering the parameters
		params = [dm_z_rs,h_z_rs,fsig]
		params = append(append((params[0])[:],(params[1])[:]),(params[2])[:])
		params = [params[0::3],params[1::3],params[2::3]]
		params = append(append((params[0])[:],(params[1])[:]),(params[2])[:])
	
	#computation of chi square
	d = params - data
	chi2 = dot(d,dot(precision,d))

	if not mode:
		if feedback:
			print ("[H*rs/rs_fid,D_m*rs_fid/rs] at redshift %lf = [%lf,%lf]"%(redshift[0],h_z_rs[0],dm_z_rs[0]))
			print ("[H*rs/rs_fid,D_m*rs_fid/rs] at redshift %lf = [%lf,%lf]"%(redshift[1],h_z_rs[1],dm_z_rs[1]))
			print ("[H*rs/rs_fid,D_m*rs_fid/rs] at redshift %lf = [%lf,%lf]"%(redshift[2],h_z_rs[2],dm_z_rs[2]))
	
	else: 
		if feedback:
			print ("[H*rs/rs_fid,D_m*rs_fid/rs,fsigma8] at redshift %lf = [%lf,%lf,%lf]"%(redshift[0],h_z_rs[0],dm_z_rs[0],fsig[0]))
			print ("[H*rs/rs_fid,D_m*rs_fid/rs,fsigma8] at redshift %lf = [%lf,%lf,%lf]"%(redshift[1],h_z_rs[1],dm_z_rs[1],fsig[1]))
			print ("[H*rs/rs_fid,D_m*rs_fid/rs,fsigma8] at redshift %lf = [%lf,%lf,%lf]"%(redshift[2],h_z_rs[2],dm_z_rs[2],fsig[2]))
	
	#ln(likelihood)
	block[likes, 'BOSS_DR12_LIKE'] = chi2*(-0.5) 
	#signal that everything went fine
	return 0

def cleanup(config):
    #nothing to do here!  We just include this 
    # for completeness.  The joy of python.
    return 0

